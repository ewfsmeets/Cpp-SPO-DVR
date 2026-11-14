// Reading the SCM data files (params.dat and the slab position file) with the list-directed reader,
// statement by statement as the Fortran's READs. Fortran: DistortedSurfacePES.dir/SCM_IO.f90
// (Read_SCM_Parameters, GetSurfaceParameters, ReadSlabFileParams, DummyReadSlabFile, ReadSCMSlab)
#include "potential/scm/scm_files.hpp"

#include <cstdint>
#include <format>
#include <optional>
#include <string>
#include <system_error>
#include <vector>

#include "core/error.hpp"
#include "fortran_text/list_directed_reader.hpp"

namespace spodvr {

namespace {

using fortran_text::ListDirectedReader;
using fortran_text::ListStatement;
using fortran_text::ReadError;

// The input keys that errors name (the Fortran's DVR.inp names).
constexpr const char* parameter_file_key = "SCMParFilepath";
constexpr const char* position_file_key = "PosFilepath";
constexpr const char* configuration_key = "PosFileSlabNum";
constexpr const char* cutoff_key = "SCMCutoff";

// The surface layers the correction uses. Fortran: SCM_Data%UsedLayers (= 2, never changed)
constexpr int used_layers = 2;
// The values of one layer in params.dat. Fortran: RydbergParams(i)%P(18)
constexpr int values_per_layer = 18;
// Lines before the first configuration of a position file. Fortran: DummyReadSlabFile
constexpr int header_lines = 3;

[[noreturn]] void fail(const char* key, const std::filesystem::path& file, const std::string& message) {
  throw InputError(key, file.string() + ": " + message);
}

// A reading error (its message names the file and line) with what was being read appended.
[[noreturn]] void fail_reading(const char* key, const ReadError& error, const std::string& what) {
  throw InputError(key, std::string(error.what()) + " (" + what + ")");
}

// Opens a data file for reading; a missing or unreadable file is an input error naming `key`.
ListDirectedReader open_reader(const char* key, const std::filesystem::path& file) {
  std::error_code ignored;
  if (!std::filesystem::is_regular_file(file, ignored)) fail(key, file, "no such file");
  try {
    return ListDirectedReader(file);
  } catch (const ReadError& e) {
    throw InputError(key, e.what());
  }
}

// ---------------------------------------------------------------------------------------------
// params.dat
// ---------------------------------------------------------------------------------------------

// Line 1 is "VERSION n": the tag in columns 1–7 (case-sensitive), n read list-directed from column
// 8 on (the first value; the rest of the line is ignored). Without the tag the Fortran assumes
// version 1 (SCM D8). Fortran: Read_SCM_Parameters, READ(SCMVersion(8:100), *)
int read_version_tag(const std::filesystem::path& file, const std::string& line) {
  const std::string tag = "VERSION";
  if (line.compare(0, tag.size(), tag) != 0) {
    const std::string shown = line.size() > 40 ? line.substr(0, 40) + "..." : line;
    fail(parameter_file_key, file,
         std::format("line 1 is '{}'; it must be 'VERSION 2' (upper case, from column 1). Without the tag the "
                     "Fortran assumes version 1, which it does not support",
                     shown));
  }
  const std::string rest = line.substr(tag.size());
  const std::size_t start = rest.find_first_not_of(" \t");
  if (start == std::string::npos) fail(parameter_file_key, file, "line 1: 'VERSION' without a number");
  const std::size_t end = rest.find_first_of(" \t,/", start);
  const ListStatement version(file.string(), {{rest.substr(start, end - start), 1}});
  return version.integer(0);
}

// One list-directed READ of a layer's 18 values, in file order P1 ... P18.
// Fortran: READ(ParFile, *) RydbergParams(i)%P
ScmLayerParameters read_layer_parameters(ListDirectedReader& reader) {
  const ListStatement values = reader.read(values_per_layer);
  ScmLayerParameters layer;
  int next = 0;
  for (double& value : layer.intercept) value = values.real(next++);
  for (double& value : layer.slope) value = values.real(next++);
  layer.r_clamp_min = values.real(next++);
  layer.r_clamp_max = values.real(next);
  return layer;
}

// P17 < P18. The Fortran does not check it: Setup_3Body_Rydberg_Par (SCM_Types.f90) clamps r with
// IF/ELSE IF, so with P17 ≥ P18 every r ≥ P17 becomes P18 (SCM D6).
void require_clamp_range(const std::filesystem::path& file, const ScmLayerParameters& p, int layer_number) {
  if (p.r_clamp_min < p.r_clamp_max) return;
  fail(parameter_file_key, file,
       std::format("layer {}: P17 = {} must be smaller than P18 = {} (the clamp range of r)", layer_number,
                   p.r_clamp_min, p.r_clamp_max));
}

// ---------------------------------------------------------------------------------------------
// The position file
// ---------------------------------------------------------------------------------------------

// The three header statements (each may span lines). Fortran: ReadSlabFileParams
ScmSlabHeader read_header(const std::filesystem::path& file) {
  ListDirectedReader reader = open_reader(position_file_key, file);
  ScmSlabHeader header;
  try {
    const ListStatement names = reader.read(2);
    header.tag = names.text(0);
    header.structure = names.text(1);
    const ListStatement sizes = reader.read(6);
    header.nx = sizes.integer(0);
    header.ny = sizes.integer(1);
    header.nz = sizes.integer(2);
    header.configuration_count = sizes.integer(3);
    header.version = sizes.integer(4);
    header.max_surface_atoms = sizes.integer(5);
    header.lattice_constant = reader.read(1).real(0);
  } catch (const ReadError& e) {
    fail_reading(position_file_key, e,
                 "the header: SlabFileTag SlabSurfStruct / nx ny nz NumSlabConfigs PosFileVersion MaxNumSurfAtoms / "
                 "DistLatticeConstant");
  }
  return header;
}

// The header checks, in the Fortran's order (GetSurfaceParameters, then ReadSCMSlab), plus the
// ones the C++ adds where the Fortran continues with garbage or silently gives ΔV = 0 (SCM D5, D6,
// D14 in docs/fortran_deviations.md).
void require_header(const std::filesystem::path& file, const ScmSlabHeader& header, int configuration) {
  if (header.structure != "FCC111") {
    fail(position_file_key, file,
         std::format("SlabSurfStruct '{}': only FCC111 is supported (case-sensitive)", header.structure));
  }
  if (configuration > header.configuration_count) {
    fail(configuration_key, file,
         std::format("configuration {} asked for, but the file has NumSlabConfigs = {}", configuration,
                     header.configuration_count));
  }
  if (header.nz < used_layers) {
    fail(position_file_key, file,
         std::format("nz = {}: the correction uses the top {} layers, so nz must be at least {}", header.nz,
                     used_layers, used_layers));
  }
  if (header.version == 1) {
    fail(position_file_key, file,
         "PosFileVersion 1 (ideal positions generated from the order of the atoms) is not supported; "
         "write the file in version 3 (distorted and ideal positions and the layer of every atom)");
  }
  if (header.version != 2 && header.version != 3) {
    fail(position_file_key, file, std::format("unknown PosFileVersion {} (supported: 3 and 2)", header.version));
  }
  if (header.version == 3 && (header.nx < 1 || header.ny < 1)) {
    fail(position_file_key, file, std::format("nx = {}, ny = {}: both must be at least 1", header.nx, header.ny));
  }
  if (header.version == 2 && header.max_surface_atoms < 1) {
    fail(position_file_key, file,
         std::format("MaxNumSurfAtoms = {} in a version 2 file: it must be at least 1 (the Fortran silently "
                     "uses no atom, ΔV = 0)",
                     header.max_surface_atoms));
  }
}

// The records of one configuration: one per atom, nx·ny·nz in version 3, MaxNumSurfAtoms in version 2.
// Fortran: the skip loops of ReadSCMSlab (numAtoms)
std::int64_t lines_per_configuration(const ScmSlabHeader& header) {
  if (header.version == 2) return header.max_surface_atoms;
  return std::int64_t{header.nx} * header.ny * header.nz;
}

// Skips the `count` records before configuration `configuration` without looking at them (empty
// READ statements). Fortran: DummyReadSlabFile and the skip loops of ReadSCMSlab, READ(loadFileUnit, *)
void skip_records(ListDirectedReader& reader, std::int64_t count, const std::filesystem::path& file,
                  int configuration) {
  for (std::int64_t i = 0; i < count; ++i) {
    if (!reader.try_read(0)) {
      fail(position_file_key, file,
           std::format("the file ends after line {}; configuration {} would begin on line {}", i, configuration,
                       count + 1));
    }
  }
}

// One atom record. Fortran: READ(loadFileUnit, *) currAtLoc(1:3), idAtLoc(1:3), zIdCor (version 3)
// and VectorsDisplaced, VectorsIdeal, AtomType (version 2): the same seven values.
ScmAtom read_atom(ListDirectedReader& reader) {
  const ListStatement values = reader.read(7);
  ScmAtom atom;
  atom.distorted = {values.real(0), values.real(1), values.real(2)};
  atom.ideal = {values.real(3), values.real(4), values.real(5)};
  atom.layer = values.integer(6);
  return atom;
}

// Version 3 keeps an atom whose distorted position lies strictly inside the cut-off sphere around
// the origin, compared squared: ((x² + y²) + z²) < c². The order must not change: the Fortran's
// expression decides the atoms on the boundary. (The Fortran first subtracts the molecule centre
// h2Cent = 0, which is exact.) Fortran: ReadSCMSlab, CASE(3)
bool inside_cutoff(const Vec3& position, double cutoff) {
  const double distance_squared = (position.x * position.x + position.y * position.y) + position.z * position.z;
  return distance_squared < cutoff * cutoff;
}

// The atoms of one configuration, read from where `reader` stands, as the correction uses them:
// every layer index checked, layers 1–2 kept, in version 3 only those inside the cut-off.
std::vector<ScmAtom> read_configuration(ListDirectedReader& reader, const std::filesystem::path& file,
                                        const ScmSlabHeader& header, int configuration, double cutoff) {
  std::vector<ScmAtom> kept;
  const std::int64_t atom_count = lines_per_configuration(header);
  for (std::int64_t i = 0; i < atom_count; ++i) {
    ScmAtom atom;
    try {
      atom = read_atom(reader);
    } catch (const ReadError& e) {
      fail_reading(position_file_key, e,
                   std::format("atom {} of configuration {}: xd yd zd xi yi zi layer", i + 1, configuration));
    }
    if (atom.layer < 1) {  // the Fortran keeps it and indexes Pars(0) (SCM D6)
      fail(position_file_key, file,
           std::format("atom {} of configuration {}: layer {}; layers are numbered from 1 (the top layer)", i + 1,
                       configuration, atom.layer));
    }
    // Layers below the second never contribute: version 3 drops them here, version 2 in the
    // Fortran's potential loop (SCM_Potential.f90 SCM_SurfInteract_3Body, before any arithmetic).
    if (atom.layer > used_layers) continue;
    if (header.version == 3 && !inside_cutoff(atom.distorted, cutoff)) continue;
    kept.push_back(atom);
  }
  return kept;
}

// At least one atom kept. Version 3: as the Fortran ("No atoms in SCM range", ReadSCMSlab).
// Version 2: the Fortran has no such check and silently gives ΔV = 0 (SCM D6).
void require_atoms(const std::filesystem::path& file, const ScmSlab& slab) {
  if (!slab.atoms.empty()) return;
  if (slab.header.version == 3) {
    fail(position_file_key, file,
         std::format("configuration {} has no atom of layers 1–{} within the cut-off of {} bohr", slab.configuration,
                     used_layers, slab.cutoff));
  }
  fail(position_file_key, file,
       std::format("configuration {} has no atom of layers 1–{} (the Fortran silently uses none, ΔV = 0)",
                   slab.configuration, used_layers));
}

}  // namespace

ScmParameters read_scm_parameters(const std::filesystem::path& file) {
  ListDirectedReader reader = open_reader(parameter_file_key, file);
  ScmParameters parameters;
  try {
    const std::optional<std::string> first_line = reader.read_line();
    if (!first_line) fail(parameter_file_key, file, "the file is empty (line 1 must be 'VERSION 2')");
    parameters.version = read_version_tag(file, *first_line);
    if (parameters.version != 2) {
      fail(parameter_file_key, file,
           std::format("line 1: VERSION {}: only version 2 is supported (the Fortran also rejects every other)",
                       parameters.version));
    }
    int layer_number = 1;
    for (ScmLayerParameters& layer : parameters.layers) {
      try {
        layer = read_layer_parameters(reader);
      } catch (const ReadError& e) {
        fail_reading(parameter_file_key, e,
                     std::format("the {} values P1–P18 of layer {}", values_per_layer, layer_number));
      }
      ++layer_number;
    }
  } catch (const ReadError& e) {
    throw InputError(parameter_file_key, e.what());
  }
  int layer_number = 1;
  for (const ScmLayerParameters& layer : parameters.layers) require_clamp_range(file, layer, layer_number++);
  return parameters;
}

ScmSlab read_scm_slab(const std::filesystem::path& file, int configuration, double cutoff) {
  if (!(cutoff > 0.0)) throw InputError(cutoff_key, std::format("the cut-off must be positive, got {} bohr", cutoff));
  if (configuration < 1) {  // the Fortran reads configuration 1 (SCM D6)
    throw InputError(configuration_key,
                     std::format("the configuration number must be at least 1, got {}", configuration));
  }

  ScmSlab slab;
  slab.configuration = configuration;
  slab.cutoff = cutoff;
  slab.header = read_header(file);
  require_header(file, slab.header, configuration);

  // The atoms start on line 4 whatever the header statements read: the Fortran reads the header
  // in GetSurfaceParameters, reopens the file and skips three records (DummyReadSlabFile).
  ListDirectedReader reader = open_reader(position_file_key, file);
  skip_records(reader, header_lines + (configuration - 1) * lines_per_configuration(slab.header), file, configuration);
  slab.atoms = read_configuration(reader, file, slab.header, configuration, cutoff);
  require_atoms(file, slab);
  return slab;
}

}  // namespace spodvr

// Reading, writing and validating a CRP data set. Fortran: CRP_Wrapper.f90 SetUp_CRP_*_PES and the
// Read_PES_Info / Read_*_Data routines of CRP_C*v_AtomSurface.f90 and CRP_C*v_MoleculeSurface.f90
// (file part; the unit conversion and the corrugation reduction are done by the PES constructors).
#include "crplib/dataset.hpp"

#include <format>
#include <fstream>
#include <string>

#include "crplib/errors.hpp"
#include "crplib/list_directed_reader.hpp"

namespace crplib {

using detail::ListDirectedReader;
using detail::ListStatement;

const char* symmetry_name(Symmetry symmetry) {
  switch (symmetry) {
    case Symmetry::c3v: return "C3v";
    case Symmetry::c4v: return "C4v";
    case Symmetry::c6v: return "C6v";
  }
  return "?";
}

// Fortran: SIZE(fouri3D) (C3v/C6v, ATII), Cuts(1..6) (C4v, Read_AtomSurface_Data)
int atom_cut_count(Symmetry symmetry) {
  switch (symmetry) {
    case Symmetry::c3v: return 10;
    case Symmetry::c4v: return 6;
    case Symmetry::c6v: return 7;
  }
  return 0;
}

// Fortran: the cut indices used by RRINTII (C3v: 1..29, C6v: 1..18) and RRINT (C4v: 1..16)
int molecule_cut_count(Symmetry symmetry) {
  switch (symmetry) {
    case Symmetry::c3v: return 29;
    case Symmetry::c4v: return 16;
    case Symmetry::c6v: return 18;
  }
  return 0;
}

namespace {

// Fortran: PESCut3D/PESCut6D FileName is CHARACTER(50); longer names are truncated silently.
constexpr std::size_t max_file_name_length = 50;

[[noreturn]] void fail(const std::filesystem::path& file, const std::string& message) {
  throw DataError(file.string() + ": " + message);
}

Symmetry parse_symmetry(const std::string& text, const std::filesystem::path& file) {
  // Fortran: Symmetry .EQ. "C3v" (exact, case-sensitive)
  if (text == "C3v") return Symmetry::c3v;
  if (text == "C4v") return Symmetry::c4v;
  if (text == "C6v") return Symmetry::c6v;
  fail(file, "unknown symmetry '" + text + "' (C3v, C4v or C6v)");
}

// count lines "x f" after the count has been read. Fortran: DO i = 1, NX; READ(10,*) X(i), F(i)
Table1D read_table(ListDirectedReader& in, int count, const std::filesystem::path& file) {
  if (count < 0) fail(file, "negative number of points " + std::to_string(count));
  Table1D table;
  table.x.resize(count);
  table.f.resize(count);
  for (int i = 0; i < count; ++i) {
    const ListStatement s = in.read(2);
    table.x[i] = s.real(0);
    table.f[i] = s.real(1);
  }
  return table;
}

// Fortran: SetUp_CRP_*_PES reading CRPData (and CRPData.info)
void read_crpdata(const std::filesystem::path& directory, DataSet& d) {
  const std::filesystem::path file = directory / "CRPData";
  if (!std::filesystem::exists(file)) {
    throw DataError("CRP data set " + directory.string() + " does not exist (no file CRPData)");
  }
  ListDirectedReader in(file);
  d.symmetry = parse_symmetry(in.read(1).text(0), file);
  const ListStatement units = in.read(2);
  try {
    d.units = Units{length_unit_from_code(units.integer(0)), energy_unit_from_code(units.integer(1))};
  } catch (const DataError& e) {
    fail(file, e.what());
  }
  d.lattice_constant = in.read(1).real(0);
  const int num_layers = in.read(1).integer(0);
  if (num_layers < 1) fail(file, "the number of layers must be at least 1");
  const ListStatement spacings = in.read(num_layers - 1);  // an empty READ still skips a record
  for (int i = 0; i < num_layers - 1; ++i) d.interlayer_spacings.push_back(spacings.real(i));
  // Fortran: READ(200, *, END=100, ERR=100) Disp1, Disp2 (atom PES only): absent or unreadable -> 0
  if (const std::optional<ListStatement> disp = in.try_read(2)) {
    try {
      d.displacement_1 = disp->real(0);
      d.displacement_2 = disp->real(1);
    } catch (const DataError&) {
      d.displacement_1 = 0.0;
      d.displacement_2 = 0.0;
    }
  }

  const std::filesystem::path info_file = directory / "CRPData.info";
  if (std::filesystem::exists(info_file)) {
    ListDirectedReader info(info_file);
    DataSetInfo text;
    if (auto line = info.read_line()) text.short_info = *line;
    if (auto line = info.read_line()) text.long_info = *line;
    if (auto line = info.read_line()) text.reference = *line;
    d.info = text;
  }
}

// Fortran: Read_PES_Info and Read_AtomSurface_Data (file part)
void read_atom_surface(const std::filesystem::path& directory, DataSet& d) {
  const std::filesystem::path list_file = directory / "CRPData.3D";
  ListDirectedReader list(list_file);
  const int n = list.read(1).integer(0);
  if (n < 0) fail(list_file, "negative number of cuts");
  for (int i = 0; i < n; ++i) {
    const ListStatement s = list.read(3);
    AtomCut cut;
    cut.x = s.real(0);
    cut.y = s.real(1);
    cut.file_name = s.text(2);
    d.atom_cuts.push_back(cut);
  }

  const std::filesystem::path corrugation_file = directory / "Atom" / "intrep.dat";
  ListDirectedReader corrugation(corrugation_file);
  d.corrugation_slope_first = corrugation.read(1).real(0);
  d.corrugation_slope_last = corrugation.read(1).real(0);
  d.corrugation = read_table(corrugation, corrugation.read(1).integer(0), corrugation_file);

  for (AtomCut& cut : d.atom_cuts) {
    if (cut.file_name.size() > max_file_name_length) {
      fail(list_file, "cut file name '" + cut.file_name + "' longer than 50 characters (the Fortran truncates it)");
    }
    const std::filesystem::path file = directory / "Atom" / cut.file_name;
    if (!std::filesystem::exists(file)) fail(file, "atom-surface cut file not found");
    ListDirectedReader in(file);
    const int count = in.read(1).integer(0);
    cut.slope_first = in.read(1).real(0);
    cut.slope_last = in.read(1).real(0);
    cut.potential = read_table(in, count, file);
  }
}

// Fortran: Read_PES_Info and Read_MoleculeSurface_Data (file part)
void read_molecule_surface(const std::filesystem::path& directory, DataSet& d) {
  const std::filesystem::path list_file = directory / "CRPData.6D";
  ListDirectedReader list(list_file);
  const ListStatement header = list.read(2);
  const int n = header.integer(0);
  if (n < 0) fail(list_file, "negative number of cuts");
  d.cuts_interpolation = header.logical(1);
  for (int i = 0; i < n; ++i) {
    const ListStatement s = list.read(5);
    MoleculeCut cut;
    cut.x = s.real(0);
    cut.y = s.real(1);
    cut.theta = s.real(2);
    cut.phi = s.real(3);
    cut.file_name = s.text(4);
    d.molecule_cuts.push_back(cut);
  }

  const std::filesystem::path extrapolation_file = directory / "Molecule" / "extrapol.dat";
  ListDirectedReader extrapolation(extrapolation_file);
  d.zso_min = extrapolation.read(1).real(0);
  d.zso_max = extrapolation.read(1).real(0);
  d.z_dependence = read_table(extrapolation, extrapolation.read(1).integer(0), extrapolation_file);

  const std::filesystem::path h2_file = directory / "Molecule" / "VH2.dat";
  ListDirectedReader h2(h2_file);
  d.h2_potential = read_table(h2, h2.read(1).integer(0), h2_file);

  for (MoleculeCut& cut : d.molecule_cuts) {
    if (cut.file_name.size() > max_file_name_length) {
      fail(list_file, "cut file name '" + cut.file_name + "' longer than 50 characters (the Fortran truncates it)");
    }
    const std::filesystem::path file = directory / "Molecule" / cut.file_name;
    // The Fortran skips a missing file (INQUIRE) and later evaluates an unallocated spline.
    if (!std::filesystem::exists(file)) fail(file, "molecule-surface cut file not found");
    ListDirectedReader in(file);
    const int nr = in.read(1).integer(0);
    const int nz = in.read(1).integer(0);
    if (nr < 0 || nz < 0) fail(file, "negative grid size");
    cut.r.resize(nr);
    cut.z.resize(nz);
    cut.v.resize(static_cast<std::size_t>(nr) * nz);
    // Rows "Z r V" with r fastest; r and Z are re-read in every row: the last row wins.
    for (int j = 0; j < nz; ++j) {
      for (int i = 0; i < nr; ++i) {
        const ListStatement s = in.read(3);
        cut.z[j] = s.real(0);
        cut.r[i] = s.real(1);
        cut.v[i + nr * j] = s.real(2);
      }
    }
  }
}

bool strictly_increasing(const std::vector<double>& x) {
  for (std::size_t i = 1; i < x.size(); ++i)
    if (!(x[i] > x[i - 1])) return false;
  return true;
}

// A spline table: at least 4 points (DSPLIN), strictly increasing knots.
void check_spline_grid(const std::vector<double>& x, const std::string& what) {
  if (x.size() < 4) throw DataError(what + ": a cubic spline needs at least 4 points, got " + std::to_string(x.size()));
  if (!strictly_increasing(x)) throw DataError(what + ": the grid is not strictly increasing");
}

}  // namespace

void validate_data_set(const DataSet& d) {
  const char* name = symmetry_name(d.symmetry);
  if (!(d.lattice_constant > 0.0)) throw DataError("CRPData: the lattice constant must be positive");
  if (d.symmetry == Symmetry::c3v && d.num_layers() < 4) {
    throw DataError("CRPData: C3v needs at least 4 layers (MURII sums over 4 layers; Fortran: 'not enough layers')");
  }
  if (d.symmetry == Symmetry::c4v && d.num_layers() < 2) {
    throw DataError("CRPData: C4v needs at least 2 layers (MURO uses the first interlayer spacing)");
  }
  if (d.symmetry == Symmetry::c3v && d.displacement_1 != 0.0) {
    throw DataError("CRPData: the C3v atom potential does not support displacements (Disp1 must be 0)");
  }

  check_spline_grid(d.corrugation.x, "Atom/intrep.dat");
  const int atom_cuts = static_cast<int>(d.atom_cuts.size());
  const int required_atom_cuts = atom_cut_count(d.symmetry);
  if (d.symmetry == Symmetry::c4v ? atom_cuts < required_atom_cuts : atom_cuts != required_atom_cuts) {
    throw DataError("CRPData.3D: " + std::string(name) + " needs " +
                    (d.symmetry == Symmetry::c4v ? "at least " : "exactly ") + std::to_string(required_atom_cuts) +
                    " atom-surface cuts, the data set has " + std::to_string(atom_cuts));
  }
  for (const AtomCut& cut : d.atom_cuts) {
    check_spline_grid(cut.potential.x, "Atom/" + cut.file_name);
    if (d.symmetry == Symmetry::c4v && cut.potential.x != d.atom_cuts.front().potential.x) {
      throw DataError("Atom/" + cut.file_name +
                      ": C4v needs the same z values for all sites (Fortran: 'Error in data for site')");
    }
  }

  if (d.zso_min > d.zso_max) throw DataError("Molecule/extrapol.dat: ZsoMin is larger than ZsoMax");
  check_spline_grid(d.z_dependence.x, "Molecule/extrapol.dat");
  check_spline_grid(d.h2_potential.x, "Molecule/VH2.dat");
  const int required_molecule_cuts = molecule_cut_count(d.symmetry);
  if (static_cast<int>(d.molecule_cuts.size()) < required_molecule_cuts) {
    throw DataError("CRPData.6D: " + std::string(name) + " needs at least " + std::to_string(required_molecule_cuts) +
                    " molecule-surface cuts, the data set has " + std::to_string(d.molecule_cuts.size()));
  }
  for (const MoleculeCut& cut : d.molecule_cuts) {
    const std::string what = "Molecule/" + cut.file_name;
    check_spline_grid(cut.r, what + " (r)");
    check_spline_grid(cut.z, what + " (Z)");
    if (cut.r.size() > 999 || cut.z.size() > 999)
      throw DataError(what + ": at most 999 points per direction (SPLIN2D)");
    if (cut.v.size() != cut.r.size() * cut.z.size()) throw DataError(what + ": wrong number of values");
  }
}

std::vector<std::string> data_set_warnings(const DataSet& d) {
  std::vector<std::string> warnings;
  const int required = molecule_cut_count(d.symmetry);
  if (static_cast<int>(d.molecule_cuts.size()) > required) {
    warnings.push_back("CRPData.6D has " + std::to_string(d.molecule_cuts.size()) + " cuts; only the first " +
                       std::to_string(required) + " are used (the others still limit Rmin, Rmax and Zmin)");
  }
  if (!d.z_dependence.x.empty() && d.z_dependence.x.back() < d.zso_max) {
    warnings.push_back("Molecule/extrapol.dat: the last Z (Zasy, the asymptotic Z of the PES) is below ZsoMax");
  }
  return warnings;
}

DataSet read_data_set(const std::filesystem::path& directory) {
  // Fortran (CRP_QD_Wrapper Read_PES_From_Disk): a file <data set>.dat selects a multi-PES set,
  // a branch that is empty in the Fortran.
  const std::filesystem::path header = directory.string() + ".dat";
  if (std::filesystem::exists(header)) {
    throw DataError(header.string() + ": multi-PES data sets (a <data set>.dat file) are not supported");
  }
  DataSet d;
  read_crpdata(directory, d);
  read_atom_surface(directory, d);
  read_molecule_surface(directory, d);
  validate_data_set(d);
  return d;
}

namespace {

// Shortest representation that reads back exactly (Fortran list-directed input accepts it).
std::string number(double x) { return std::format("{}", x); }

std::string quoted(const std::string& name) {
  if (name.find_first_of(" ,/'\"\t") == std::string::npos) return name;
  std::string q = "'";
  for (char c : name) q += (c == '\'' ? std::string("''") : std::string(1, c));
  return q + "'";
}

std::ofstream open_for_writing(const std::filesystem::path& file) {
  std::filesystem::create_directories(file.parent_path());
  std::ofstream out(file);
  if (!out) throw DataError(file.string() + ": cannot write the file");
  return out;
}

void write_table(std::ofstream& out, const Table1D& table) {
  out << table.x.size() << '\n';
  for (std::size_t i = 0; i < table.x.size(); ++i) out << number(table.x[i]) << ' ' << number(table.f[i]) << '\n';
}

}  // namespace

void write_data_set(const DataSet& d, const std::filesystem::path& directory) {
  {
    std::ofstream out = open_for_writing(directory / "CRPData");
    out << symmetry_name(d.symmetry) << '\n'
        << static_cast<int>(d.units.length) << ' ' << static_cast<int>(d.units.energy) << '\n'
        << number(d.lattice_constant) << '\n'
        << d.num_layers() << '\n';
    for (double spacing : d.interlayer_spacings) out << number(spacing) << ' ';
    out << '\n' << number(d.displacement_1) << ' ' << number(d.displacement_2) << '\n';
  }
  if (d.info) {
    std::ofstream out = open_for_writing(directory / "CRPData.info");
    out << d.info->short_info << '\n' << d.info->long_info << '\n' << d.info->reference << '\n';
  }
  {
    std::ofstream out = open_for_writing(directory / "CRPData.3D");
    out << d.atom_cuts.size() << '\n';
    for (const AtomCut& cut : d.atom_cuts)
      out << number(cut.x) << ' ' << number(cut.y) << ' ' << quoted(cut.file_name) << '\n';
  }
  {
    std::ofstream out = open_for_writing(directory / "Atom" / "intrep.dat");
    out << number(d.corrugation_slope_first) << '\n' << number(d.corrugation_slope_last) << '\n';
    write_table(out, d.corrugation);
  }
  for (const AtomCut& cut : d.atom_cuts) {
    std::ofstream out = open_for_writing(directory / "Atom" / cut.file_name);
    out << cut.potential.x.size() << '\n' << number(cut.slope_first) << '\n' << number(cut.slope_last) << '\n';
    for (std::size_t i = 0; i < cut.potential.x.size(); ++i)
      out << number(cut.potential.x[i]) << ' ' << number(cut.potential.f[i]) << '\n';
  }
  {
    std::ofstream out = open_for_writing(directory / "CRPData.6D");
    out << d.molecule_cuts.size() << ' ' << (d.cuts_interpolation ? 'T' : 'F') << '\n';
    for (const MoleculeCut& cut : d.molecule_cuts) {
      out << number(cut.x) << ' ' << number(cut.y) << ' ' << number(cut.theta) << ' ' << number(cut.phi) << ' '
          << quoted(cut.file_name) << '\n';
    }
  }
  {
    std::ofstream out = open_for_writing(directory / "Molecule" / "extrapol.dat");
    out << number(d.zso_min) << '\n' << number(d.zso_max) << '\n';
    write_table(out, d.z_dependence);
  }
  {
    std::ofstream out = open_for_writing(directory / "Molecule" / "VH2.dat");
    write_table(out, d.h2_potential);
  }
  for (const MoleculeCut& cut : d.molecule_cuts) {
    std::ofstream out = open_for_writing(directory / "Molecule" / cut.file_name);
    out << cut.r.size() << '\n' << cut.z.size() << '\n';
    for (std::size_t j = 0; j < cut.z.size(); ++j)
      for (std::size_t i = 0; i < cut.r.size(); ++i)
        out << number(cut.z[j]) << ' ' << number(cut.r[i]) << ' ' << number(cut.v[i + cut.r.size() * j]) << '\n';
  }
}

}  // namespace crplib

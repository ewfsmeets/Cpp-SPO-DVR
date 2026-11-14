// Tests: the SCM reference data (tests/scm/reference/, tools/fortran_reference/run_scm_reference.sh)
// are complete and cover what they should, and the C++ readers read their files bit for bit as the
// Fortran did (SCM_IO.f90 Read_SCM_Parameters, ReadSlabFileParams, ReadSCMSlab).
#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "potential/scm/scm_files.hpp"
#include "scm_driver_output.hpp"

using namespace spodvr;
using scm_test::ScmDriverOutput;

namespace {

// tests/scm/reference: the driver outputs and the synthetic inputs.
const std::filesystem::path reference_directory = SPODVR_SCM_REFERENCE_DIR;
// potentials/H2Cu111SRP_C3v.dir: the example params.dat and examplePosfile.dat.
const std::filesystem::path example_directory = SPODVR_SCM_EXAMPLE_DIR;

// One driver run of run_scm_reference.sh: its output and the files it read.
struct ReferenceRun {
  std::filesystem::path output;      // relative to tests/scm/reference
  std::filesystem::path parameters;  // params.dat as run
  std::filesystem::path slab;        // the position file as run
};

const std::vector<ReferenceRun>& reference_runs() {
  static const std::vector<ReferenceRun> runs = {
      {"example/rescaling_on.out", example_directory / "params.dat", example_directory / "examplePosfile.dat"},
      {"example/rescaling_off.out", example_directory / "params.dat", example_directory / "examplePosfile.dat"},
      {"synthetic/v3_slab3_rescaling_on.out", reference_directory / "synthetic/params.in",
       reference_directory / "synthetic/slab_v3.in"},
      {"synthetic/v2_slab2_rescaling_off.out", reference_directory / "synthetic/params.in",
       reference_directory / "synthetic/slab_v2.in"},
  };
  return runs;
}

// The bit pattern of a double: bitwise comparisons with EXPECT_EQ.
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

// The number on the first line of a point list.
int point_count(const std::filesystem::path& file) {
  std::ifstream in(file);
  int count = -1;
  in >> count;
  return count;
}

// The driver output of a run; a missing file fails (it would otherwise hide an incomplete tree).
ScmDriverOutput read_run(const ReferenceRun& run) {
  const std::filesystem::path file = reference_directory / run.output;
  EXPECT_TRUE(std::filesystem::exists(file)) << file << " missing: tools/fortran_reference/run_scm_reference.sh";
  return std::filesystem::exists(file) ? scm_test::read_scm_driver_output(file) : ScmDriverOutput{};
}

}  // namespace

// Every run finished and has every value of every point and atom, so tests can index them freely.
TEST(ReferenceFiles, DriverRunsAreComplete) {
  for (const ReferenceRun& run : reference_runs()) {
    SCOPED_TRACE(run.output.string());
    const ScmDriverOutput output = read_run(run);
    EXPECT_TRUE(output.done);
    const int points = point_count(reference_directory / run.output.parent_path() / "points.in");
    EXPECT_GT(points, 0);
    EXPECT_EQ(output.text.at("points"), std::to_string(points));
    for (const char* tag : {"x", "y", "delta_v", "v_crp", "v_total"})
      EXPECT_EQ(std::ssize(output.values.at(tag)), points) << tag;
    const int atoms = std::stoi(output.text.at("atoms"));
    EXPECT_GT(atoms, 0);
    EXPECT_EQ(std::ssize(output.atom_layers), atoms);
    for (const char* tag : {"distorted_x", "distorted_y", "distorted_z", "ideal_x", "ideal_y", "ideal_z"})
      EXPECT_EQ(std::ssize(output.values.at(tag)), atoms) << tag;
    EXPECT_EQ(output.values.at("parameter_layer1").size(), 18U);
    EXPECT_EQ(output.values.at("parameter_layer2").size(), 18U);
  }
}

// Every value of every output is finite: a NaN reference would pass a bitwise EXPECT_EQ against a
// NaN result.
TEST(ReferenceFiles, EveryValueIsFinite) {
  for (const ReferenceRun& run : reference_runs()) {
    SCOPED_TRACE(run.output.string());
    const ScmDriverOutput output = read_run(run);
    std::size_t count = 0;
    for (const auto& [tag, values] : output.values) {
      for (std::size_t i = 0; i < values.size(); ++i) EXPECT_TRUE(std::isfinite(values[i])) << tag << " " << i + 1;
      count += values.size();
    }
    EXPECT_GT(count, 0U);
  }
}

// The cases cover what the example files cannot: rescaling on and off, a configuration other than
// the first in both slab versions, layer-2 parameters unlike layer 1's (they are equal in the
// example, which would hide a mix-up of the layers), version-2 atoms below layer 2 (read by the
// Fortran, skipped by its correction) and version-2 atoms beyond the run's cut-off.
TEST(ReferenceFiles, CasesCoverRescalingConfigurationsLayersAndDeepAtoms) {
  const ScmDriverOutput on = read_run(reference_runs()[0]);
  const ScmDriverOutput off = read_run(reference_runs()[1]);
  const ScmDriverOutput v3 = read_run(reference_runs()[2]);
  const ScmDriverOutput v2 = read_run(reference_runs()[3]);
  EXPECT_EQ(on.text.at("rescaling"), "T");
  EXPECT_EQ(off.text.at("rescaling"), "F");
  EXPECT_NE(bits(on.value("scale")), bits(1.0));
  EXPECT_EQ(v3.text.at("slab_version"), "3");
  EXPECT_EQ(v2.text.at("slab_version"), "2");
  EXPECT_GT(std::stoi(v3.text.at("slab_number")), 1);
  EXPECT_GT(std::stoi(v2.text.at("slab_number")), 1);
  EXPECT_NE(v3.values.at("parameter_layer1"), v3.values.at("parameter_layer2"));
  EXPECT_TRUE(std::ranges::any_of(v2.atom_layers, [](int layer) { return layer > 2; }));
  // Version 2 has no cut-off: atoms of layers 1-2 lie outside the run's cut-off sphere, so a reader
  // that applied it would fail the slab test below.
  const double cutoff = v2.value("cutoff");
  int outside = 0;
  for (int i = 0; i < std::ssize(v2.atom_layers); ++i) {
    const double x = v2.values.at("distorted_x")[i];
    const double y = v2.values.at("distorted_y")[i];
    const double z = v2.values.at("distorted_z")[i];
    if (v2.atom_layers[i] <= 2 && (x * x + y * y) + z * z >= cutoff * cutoff) ++outside;
  }
  EXPECT_GT(outside, 0);
}

// read_scm_parameters gives the 36 values the Fortran read, bit for bit.
TEST(ReferenceFiles, ParameterFilesAreReadAsTheFortranReadsThem) {
  for (const ReferenceRun& run : reference_runs()) {
    SCOPED_TRACE(run.output.string());
    const ScmDriverOutput output = read_run(run);
    const ScmParameters parameters = read_scm_parameters(run.parameters);
    EXPECT_EQ(std::to_string(parameters.version), output.text.at("parameter_version"));
    for (int layer = 0; layer < 2; ++layer) {
      const std::vector<double>& fortran = output.values.at("parameter_layer" + std::to_string(layer + 1));
      const ScmLayerParameters& cpp = parameters.layers[layer];
      for (int k = 0; k < 8; ++k) {
        EXPECT_EQ(bits(cpp.intercept[k]), bits(fortran[k])) << "layer " << layer + 1 << " P" << k + 1;
        EXPECT_EQ(bits(cpp.slope[k]), bits(fortran[8 + k])) << "layer " << layer + 1 << " P" << k + 9;
      }
      EXPECT_EQ(bits(cpp.r_clamp_min), bits(fortran[16])) << "layer " << layer + 1 << " P17";
      EXPECT_EQ(bits(cpp.r_clamp_max), bits(fortran[17])) << "layer " << layer + 1 << " P18";
    }
  }
}

// read_scm_slab, given the run's configuration and cut-off, gives the header and the atoms the
// Fortran kept, bit for bit and in file order. The Fortran's list also holds the version-2 atoms
// below layer 2 (its correction skips them); the C++ reader drops them.
TEST(ReferenceFiles, SlabFilesAreReadAsTheFortranReadsThem) {
  for (const ReferenceRun& run : reference_runs()) {
    SCOPED_TRACE(run.output.string());
    const ScmDriverOutput output = read_run(run);
    const ScmSlab slab = read_scm_slab(run.slab, std::stoi(output.text.at("slab_number")), output.value("cutoff"));
    const ScmSlabHeader& header = slab.header;
    EXPECT_EQ(header.tag, output.text.at("slab_tag"));
    EXPECT_EQ(header.structure, output.text.at("slab_structure"));
    EXPECT_EQ(std::format("{} {} {}", header.nx, header.ny, header.nz), output.text.at("slab_size"));
    EXPECT_EQ(std::to_string(header.configuration_count), output.text.at("slab_configurations"));
    EXPECT_EQ(std::to_string(header.version), output.text.at("slab_version"));
    EXPECT_EQ(std::to_string(header.max_surface_atoms), output.text.at("slab_max_atoms"));
    EXPECT_EQ(bits(header.lattice_constant), bits(output.value("slab_lattice_constant")));

    std::vector<int> fortran_atoms;  // the indices of the Fortran's atoms of layers 1 and 2
    for (int i = 0; i < std::ssize(output.atom_layers); ++i)
      if (output.atom_layers[i] <= 2) fortran_atoms.push_back(i);
    ASSERT_EQ(slab.atoms.size(), fortran_atoms.size());
    for (int k = 0; k < std::ssize(fortran_atoms); ++k) {
      const int i = fortran_atoms[k];
      const ScmAtom& atom = slab.atoms[k];
      EXPECT_EQ(atom.layer, output.atom_layers[i]) << "atom " << i + 1;
      EXPECT_EQ(bits(atom.distorted.x), bits(output.values.at("distorted_x")[i])) << "atom " << i + 1;
      EXPECT_EQ(bits(atom.distorted.y), bits(output.values.at("distorted_y")[i])) << "atom " << i + 1;
      EXPECT_EQ(bits(atom.distorted.z), bits(output.values.at("distorted_z")[i])) << "atom " << i + 1;
      EXPECT_EQ(bits(atom.ideal.x), bits(output.values.at("ideal_x")[i])) << "atom " << i + 1;
      EXPECT_EQ(bits(atom.ideal.y), bits(output.values.at("ideal_y")[i])) << "atom " << i + 1;
      EXPECT_EQ(bits(atom.ideal.z), bits(output.values.at("ideal_z")[i])) << "atom " << i + 1;
    }
  }
}

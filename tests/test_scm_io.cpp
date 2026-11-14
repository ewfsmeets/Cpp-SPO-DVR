// Tests: the SCM data-file readers (src/potential/scm/scm_files.*): the example files of the SRP48
// data set, multi-configuration files written here (block skipping), every error, and the
// list-directed reading (records spanning lines, D exponents, skipped records).
#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <format>
#include <functional>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "potential/scm/scm_files.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using test_support::TemporaryDirectory;

namespace {

const std::filesystem::path example_directory = SPODVR_SCM_EXAMPLE_DIR;

// The 18 values of one layer in the example params.dat (both layers are equal there).
const std::string example_layer =
    "2.28970941935 1.29096748614 -0.0703557026734 0.0234889163443 -0.0632740237351 0.0272448940344 "
    "7.44024869632 10.4600719899 0.0236010392861 -0.00311521325544 0.0165908842195 0.0286500845806 "
    "-0.00722325023675 -0.00644192925251 0.000769448712731 4.63492299557e-05 0.725 1.847";

// A params.dat with the given first line and the example's two layers.
std::string parameter_file(const std::string& first_line) {
  return first_line + "\n" + example_layer + "\n" + example_layer + "\n";
}

// Expects `read` to throw an InputError naming `key`, whose message contains each of `words`.
void expect_input_error(const std::function<void()>& read, const std::string& key,
                        const std::vector<std::string>& words) {
  try {
    read();
    ADD_FAILURE() << "expected an InputError naming " << key;
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), key) << e.what();
    for (const std::string& word : words) {
      EXPECT_NE(std::string(e.what()).find(word), std::string::npos) << "'" << word << "' not in: " << e.what();
    }
  }
}

// params.dat with `text` must be rejected; the message names the file and each of `words`.
void expect_parameter_file_error(const std::string& text, const std::vector<std::string>& words) {
  const TemporaryDirectory directory("scm_io");
  const std::filesystem::path file = directory.write("params.dat", text);
  std::vector<std::string> expected = words;
  expected.push_back(file.string());
  expect_input_error([&] { (void)read_scm_parameters(file); }, "SCMParFilepath", expected);
}

// A position file with `text` must be rejected for (configuration, cutoff), naming `key`, the file
// and each of `words`.
void expect_position_file_error(const std::string& text, int configuration, double cutoff, const std::string& key,
                                const std::vector<std::string>& words) {
  const TemporaryDirectory directory("scm_io");
  const std::filesystem::path file = directory.write("slab.dat", text);
  std::vector<std::string> expected = words;
  expected.push_back(file.string());
  expect_input_error([&] { (void)read_scm_slab(file, configuration, cutoff); }, key, expected);
}

// read_scm_slab on a position file holding `text`.
ScmSlab read_slab_text(const std::string& text, int configuration, double cutoff) {
  const TemporaryDirectory directory("scm_io");
  return read_scm_slab(directory.write("slab.dat", text), configuration, cutoff);
}

// A version 3 file of one configuration, 1 × 1 × 2 atoms: layer 1 at the origin, layer 2 below it.
const std::string small_v3_file =
    "small FCC111\n 1 1 2 1 3 2\n 4.8\n"
    " 0.0 0.0 0.1  0.0 0.0 0.0 1\n"
    " 1.0 1.0 -3.9 1.0 1.0 -4.0 2\n";

}  // namespace

// ---------------------------------------------------------------------------------------------
// The example files (potentials/H2Cu111SRP_C3v.dir). The expected values are the file's text. The
// Fortran reads the same values bit for bit and keeps the same atoms at 10, 16 and 20 bohr (checked
// 2026-11-14 with a one-off gfortran driver against the r1713-dce0ecc objects).
// ---------------------------------------------------------------------------------------------

TEST(ScmParameterFile, ReadsTheExample) {
  const ScmParameters p = read_scm_parameters(example_directory / "params.dat");
  EXPECT_EQ(p.version, 2);
  const ScmLayerParameters& top = p.layers[0];
  EXPECT_EQ(top.intercept[0], 2.28970941935);
  EXPECT_EQ(top.intercept[1], 1.29096748614);
  EXPECT_EQ(top.intercept[2], -0.0703557026734);
  EXPECT_EQ(top.intercept[7], 10.4600719899);
  EXPECT_EQ(top.slope[0], 0.0236010392861);
  EXPECT_EQ(top.slope[7], 4.63492299557e-05);
  EXPECT_EQ(top.r_clamp_min, 0.725);
  EXPECT_EQ(top.r_clamp_max, 1.847);
  EXPECT_EQ(p.layers[1].intercept, top.intercept);
  EXPECT_EQ(p.layers[1].slope, top.slope);
  EXPECT_EQ(p.layers[1].r_clamp_min, top.r_clamp_min);
  EXPECT_EQ(p.layers[1].r_clamp_max, top.r_clamp_max);
}

TEST(ScmPositionFile, ReadsTheExampleHeader) {
  const ScmSlab slab = read_scm_slab(example_directory / "examplePosfile.dat", 1, 16.0);
  const ScmSlabHeader& h = slab.header;
  EXPECT_EQ(h.tag, "Multiseedv3_New");
  EXPECT_EQ(h.structure, "FCC111");
  EXPECT_EQ(h.nx, 7);
  EXPECT_EQ(h.ny, 7);
  EXPECT_EQ(h.nz, 7);
  EXPECT_EQ(h.configuration_count, 1);
  EXPECT_EQ(h.version, 3);
  EXPECT_EQ(h.max_surface_atoms, 343);
  EXPECT_EQ(h.lattice_constant, 4.8896663475);
  EXPECT_EQ(slab.configuration, 1);
  EXPECT_EQ(slab.cutoff, 16.0);
}

// Version 3 keeps the atoms of layers 1 and 2 inside the cut-off sphere, in file order.
TEST(ScmPositionFile, KeepsTheExampleAtomsInsideTheCutoff) {
  struct Case {
    double cutoff;
    int top_layer;
    int second_layer;
  };
  for (const Case c : {Case{10.0, 18, 12}, Case{16.0, 37, 36}, Case{20.0, 43, 42}}) {
    const ScmSlab slab = read_scm_slab(example_directory / "examplePosfile.dat", 1, c.cutoff);
    std::array<int, 2> counts{};
    int previous_layer = 1;
    for (const ScmAtom& atom : slab.atoms) {
      ASSERT_TRUE(atom.layer == 1 || atom.layer == 2);
      ++counts[atom.layer - 1];
      EXPECT_GE(atom.layer, previous_layer) << "file order (the file lists the atoms layer by layer)";
      previous_layer = atom.layer;
    }
    EXPECT_EQ(counts[0], c.top_layer) << c.cutoff;
    EXPECT_EQ(counts[1], c.second_layer) << c.cutoff;
  }
  // At 16 bohr: 73 atoms, the first from line 7 of the file, the last from line 97.
  const ScmSlab slab = read_scm_slab(example_directory / "examplePosfile.dat", 1, 16.0);
  ASSERT_EQ(slab.atoms.size(), 73U);
  const ScmAtom& first = slab.atoms.front();
  EXPECT_EQ(first.distorted.x, -7.3095338452999998);
  EXPECT_EQ(first.distorted.y, -12.565286989100001);
  EXPECT_EQ(first.distorted.z, -0.12005619700000000);
  EXPECT_EQ(first.ideal.x, -7.3742307077832585);
  EXPECT_EQ(first.ideal.y, -12.772542252615207);
  EXPECT_EQ(first.ideal.z, 0.0);
  EXPECT_EQ(first.layer, 1);
  const ScmAtom& last = slab.atoms.back();
  EXPECT_EQ(last.distorted.x, 5.0496853563000004);
  EXPECT_EQ(last.distorted.z, -4.1750847927999999);
  EXPECT_EQ(last.ideal.z, -3.9630013477364412);
  EXPECT_EQ(last.layer, 2);
}

// ---------------------------------------------------------------------------------------------
// Selection of the atoms and of the configuration
// ---------------------------------------------------------------------------------------------

// The cut-off test is strict and on the squared distance of the distorted position:
// (3, 4, 0) lies exactly on a sphere of 5 bohr, (9 + 16) + 0 = 25, and is not kept.
TEST(ScmPositionFile, CutoffIsStrictAndUsesTheDistortedPosition) {
  const std::string file =
      "cut FCC111\n 2 2 2 1 3 8\n 4.8\n"
      " 0.0 0.0 0.5   9.0 9.0 9.0 1\n"   // kept: distorted inside, ideal outside
      " 3.0 4.0 0.0   3.0 4.0 0.0 1\n"   // on the sphere
      " 9.0 0.0 0.0   0.0 0.0 0.0 1\n"   // dropped: distorted outside, ideal inside
      " 0.0 1.0 0.0   0.0 1.0 0.0 3\n"   // dropped: layer 3
      " 1.0 0.0 -1.0  1.0 0.0 -1.0 2\n"  // kept
      " 0.0 0.0 -9.0  0.0 0.0 -9.0 2\n"  // dropped: outside
      " 0.0 2.0 0.0   0.0 2.0 0.0 4\n"   // dropped: layer 4
      " 2.0 0.0 0.0   2.0 0.0 0.0 1\n";  // kept
  const ScmSlab on_sphere = read_slab_text(file, 1, 5.0);
  ASSERT_EQ(on_sphere.atoms.size(), 3U);
  EXPECT_EQ(on_sphere.atoms[0].distorted.z, 0.5);
  EXPECT_EQ(on_sphere.atoms[0].ideal.x, 9.0);
  EXPECT_EQ(on_sphere.atoms[1].layer, 2);
  EXPECT_EQ(on_sphere.atoms[2].distorted.x, 2.0);
  EXPECT_EQ(read_slab_text(file, 1, 5.000001).atoms.size(), 4U);
}

// The Fortran's expression ((x² + y²) + z²) < c² decides atoms on the boundary of a 16 bohr sphere
// (points found by search; the other forms evaluated in IEEE double):
// - dropped, but kept by x² + (y² + z²);
// - kept, but dropped by hypot(x, y, z) < c;
// - dropped, but kept by the FMA contractions fma(z, z, fma(x, x, y²)), fma(z, z, fma(y, y, x²))
//   and fma(z, z, x² + y²), which -ffast-math may produce (hence the strict flags of this code).
TEST(ScmPositionFile, CutoffUsesTheFortranExpression) {
  const std::string file =
      "edge FCC111\n 1 2 2 1 3 4\n 4.8\n"
      " 8.245996358163023 -9.64247267344831 -9.748141607658598    0.0 0.0 0.0 1\n"
      " 4.656727750282793 8.869184482886643 -12.47607523495546    0.0 0.0 0.0 1\n"
      " 2.765846690377188 -9.558365096724112 -12.529475199027923  0.0 0.0 0.0 2\n"
      " 0.0 0.0 -4.0  0.0 0.0 -4.0 2\n";
  const ScmSlab slab = read_slab_text(file, 1, 16.0);
  ASSERT_EQ(slab.atoms.size(), 2U);
  EXPECT_EQ(slab.atoms[0].distorted.x, 4.656727750282793);
  EXPECT_EQ(slab.atoms[1].distorted.z, -4.0);
}

// A version 3 file with three configurations of nx·ny·nz = 4 atoms; MaxNumSurfAtoms (1) is not
// what version 3 skips by. Atom a of configuration c lies at distorted (c, a, 0).
TEST(ScmPositionFile, Version3SkipsNxNyNzRecordsPerConfiguration) {
  std::string file = "multi FCC111\n 2 1 2 3 3 1\n 4.8\n";
  for (int c = 1; c <= 3; ++c) {
    for (int a = 1; a <= 4; ++a) {
      file += std::format("{}.0 {}.0 0.0  0.0 0.0 0.0 {}\n", c, a, a <= 2 ? 1 : 2);
    }
  }
  for (int c = 1; c <= 3; ++c) {
    const ScmSlab slab = read_slab_text(file, c, 100.0);
    ASSERT_EQ(slab.atoms.size(), 4U) << c;
    for (int a = 0; a < 4; ++a) {
      EXPECT_EQ(slab.atoms[a].distorted.x, c);
      EXPECT_EQ(slab.atoms[a].distorted.y, a + 1);
      EXPECT_EQ(slab.atoms[a].layer, a < 2 ? 1 : 2);
    }
  }
}

// Version 2 skips MaxNumSurfAtoms records per configuration (not nx·ny·nz = 343), applies no
// cut-off, and leaves out layers below the second (the Fortran's potential loop skips them).
TEST(ScmPositionFile, Version2SkipsMaxNumSurfAtomsRecordsAndHasNoCutoff) {
  std::string file = "v2 FCC111\n 7 7 7 3 2 3\n 4.8\n";
  for (int c = 1; c <= 3; ++c) {
    for (int a = 1; a <= 3; ++a) {
      file += std::format("{}.0 {}.0 0.0  1.0 2.0 3.0 {}\n", 10 * c, a, a);
    }
  }
  for (int c = 1; c <= 3; ++c) {
    const ScmSlab slab = read_slab_text(file, c, 1.0);
    ASSERT_EQ(slab.atoms.size(), 2U) << c;
    EXPECT_EQ(slab.atoms[0].distorted.x, 10.0 * c);
    EXPECT_EQ(slab.atoms[0].layer, 1);
    EXPECT_EQ(slab.atoms[1].distorted.y, 2.0);
    EXPECT_EQ(slab.atoms[1].ideal.z, 3.0);
    EXPECT_EQ(slab.atoms[1].layer, 2);
  }
}

// ---------------------------------------------------------------------------------------------
// List-directed reading as gfortran's READ(unit,*)
// ---------------------------------------------------------------------------------------------

// A layer's 18 values may span lines (blank lines included) with D exponents and commas; the
// rest of the statement's last line is ignored, and the next layer starts on a new line.
TEST(ScmParameterFile, RecordsSpanLinesWithDExponents) {
  const std::string text =
      "VERSION 2 the rest of the line is ignored\n"
      "2.28970941935D0 1.29096748614d+00, -7.03557026734D-02\n"
      "\n"
      "0.0234889163443 -0.0632740237351 0.0272448940344 7.44024869632 10.4600719899\n"
      "0.0236010392861 -0.00311521325544 0.0165908842195 0.0286500845806 -0.00722325023675\n"
      "-0.00644192925251 0.000769448712731 4.63492299557D-05 0.725, 1.847 99.0 not-a-number\n" +
      example_layer + "\n";
  const TemporaryDirectory directory("scm_io");
  const ScmParameters p = read_scm_parameters(directory.write("params.dat", text));
  const ScmParameters example = read_scm_parameters(example_directory / "params.dat");
  for (int layer = 0; layer < 2; ++layer) {
    EXPECT_EQ(p.layers[layer].intercept, example.layers[layer].intercept) << layer;
    EXPECT_EQ(p.layers[layer].slope, example.layers[layer].slope) << layer;
    EXPECT_EQ(p.layers[layer].r_clamp_min, example.layers[layer].r_clamp_min) << layer;
    EXPECT_EQ(p.layers[layer].r_clamp_max, example.layers[layer].r_clamp_max) << layer;
  }
}

// An atom record may span lines and use D exponents; skipped configurations are not parsed (an
// empty READ skips a record whatever it holds).
TEST(ScmPositionFile, RecordsSpanLinesAndSkippedRecordsAreNotParsed) {
  const std::string file =
      "span FCC111\n 1 1 2 2 3 2\n 4.8896663475D0\n"
      "not an atom\n"
      "also not an atom\n"
      " 0.0 0.0 1.5D-1\n"
      "    0.0 0.0 0.0 1 trailing text\n"
      " 1.0d0 1.0 -3.9 1.0 1.0 -4.0 2\n";
  const ScmSlab slab = read_slab_text(file, 2, 16.0);
  EXPECT_EQ(slab.header.lattice_constant, 4.8896663475);
  ASSERT_EQ(slab.atoms.size(), 2U);
  EXPECT_EQ(slab.atoms[0].distorted.z, 0.15);
  EXPECT_EQ(slab.atoms[0].layer, 1);
  EXPECT_EQ(slab.atoms[1].distorted.x, 1.0);
  EXPECT_EQ(slab.atoms[1].layer, 2);
}

// The atoms start on line 4 even when a header statement spans lines: the Fortran reads the
// header in GetSurfaceParameters and then skips three records (DummyReadSlabFile). Here the
// structure is on line 2, so the first atom record starts with the a_T line.
TEST(ScmPositionFile, AtomsStartOnLine4LikeTheFortran) {
  const std::string file =
      "tag\nFCC111\n 1 1 2 1 2 2\n"
      " 4.8\n"
      " 0.0 0.1 1.0 2.0 3.0 1 the rest is ignored\n"
      " 1.0 1.0 -3.9 1.0 1.0 -4.0 2\n";
  const ScmSlab slab = read_slab_text(file, 1, 16.0);
  EXPECT_EQ(slab.header.structure, "FCC111");
  EXPECT_EQ(slab.header.lattice_constant, 4.8);
  ASSERT_EQ(slab.atoms.size(), 2U);
  EXPECT_EQ(slab.atoms[0].distorted.x, 4.8);
  EXPECT_EQ(slab.atoms[0].ideal.z, 3.0);
  EXPECT_EQ(slab.atoms[1].distorted.z, -3.9);
}

// SCM D6: the Fortran passes the paths through CHARACTER(100) (cut there); the C++ does not.
TEST(ScmFiles, ScmD6LongPathsAreNotTruncated) {
  const TemporaryDirectory directory("scm_io");
  const std::string long_name(120, 'd');
  const std::filesystem::path parameters = directory.write(long_name + "/params.dat", parameter_file("VERSION 2"));
  const std::filesystem::path positions = directory.write(long_name + "/slab.dat", small_v3_file);
  ASSERT_GT(positions.string().size(), 100U);
  EXPECT_EQ(read_scm_parameters(parameters).version, 2);
  EXPECT_EQ(read_scm_slab(positions, 1, 16.0).atoms.size(), 2U);
}

// ---------------------------------------------------------------------------------------------
// Errors in params.dat
// ---------------------------------------------------------------------------------------------

// SCM D8: line 1 is "VERSION n", the tag case-sensitive in columns 1-7, n an integer.
TEST(ScmParameterFileErrors, ScmD8VersionTag) {
  expect_parameter_file_error(example_layer + "\n" + example_layer + "\n", {"line 1", "VERSION 2"});
  expect_parameter_file_error(parameter_file("version 2"), {"line 1 is 'version 2'", "upper case"});
  expect_parameter_file_error(parameter_file(" VERSION 2"), {"line 1", "column 1"});
  expect_parameter_file_error(parameter_file("VERSION"), {"line 1", "without a number"});
  expect_parameter_file_error(parameter_file("VERSION two"), {":1:", "'two' is not an integer"});
  expect_parameter_file_error(parameter_file("VERSION 2.0"), {":1:", "'2.0' is not an integer"});
  expect_parameter_file_error("", {"empty"});
}

// SCM D8: as in the Fortran, only version 2 (a file without the tag would be version 1).
TEST(ScmParameterFileErrors, ScmD8OnlyVersion2) {
  expect_parameter_file_error(parameter_file("VERSION 1"), {"VERSION 1", "only version 2"});
  expect_parameter_file_error(parameter_file("VERSION 3"), {"VERSION 3", "only version 2"});
}

TEST(ScmParameterFileErrors, ShortFileOrBadValue) {
  expect_parameter_file_error("VERSION 2\n" + example_layer + "\n", {"end of the file", "P1–P18 of layer 2"});
  expect_parameter_file_error("VERSION 2\n1 2 3\n", {"end of the file", "P1–P18 of layer 1"});
  expect_parameter_file_error("VERSION 2\n1 2 x3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18\n" + example_layer + "\n",
                              {":2:", "'x3' is not a real number", "P1–P18 of layer 1"});
}

// SCM D6: P17 < P18 (the Fortran does not check it; Setup_3Body_Rydberg_Par clamps with IF/ELSE IF).
TEST(ScmParameterFileErrors, ScmD6ClampRangeMustBeIncreasing) {
  const std::string equal = example_layer.substr(0, example_layer.rfind(" 0.725")) + " 1.0 1.0";
  const std::string inverted = example_layer.substr(0, example_layer.rfind(" 0.725")) + " 1.847 0.725";
  expect_parameter_file_error("VERSION 2\n" + equal + "\n" + example_layer + "\n", {"layer 1", "P17", "P18"});
  expect_parameter_file_error("VERSION 2\n" + example_layer + "\n" + inverted + "\n", {"layer 2", "P17", "P18"});
}

TEST(ScmParameterFileErrors, MissingFile) {
  const TemporaryDirectory directory("scm_io");
  const std::filesystem::path file = directory.path() / "missing.dat";
  expect_input_error([&] { (void)read_scm_parameters(file); }, "SCMParFilepath", {file.string(), "no such file"});
}

// ---------------------------------------------------------------------------------------------
// Errors in the position file and its selection
// ---------------------------------------------------------------------------------------------

TEST(ScmPositionFileErrors, Structure) {
  expect_position_file_error("s BCC110\n 1 1 2 1 3 2\n 4.8\n", 1, 16.0, "PosFilepath", {"SlabSurfStruct", "BCC110"});
  expect_position_file_error("s fcc111\n 1 1 2 1 3 2\n 4.8\n", 1, 16.0, "PosFilepath", {"SlabSurfStruct", "fcc111"});
}

// SCM D5: versions 3 and 2 only.
TEST(ScmPositionFileErrors, ScmD5OnlyVersions3And2) {
  expect_position_file_error("s FCC111\n 1 1 2 1 1 2\n 4.8\n", 1, 16.0, "PosFilepath",
                             {"PosFileVersion 1", "not supported"});
  expect_position_file_error("s FCC111\n 1 1 2 1 4 2\n 4.8\n", 1, 16.0, "PosFilepath", {"PosFileVersion 4"});
}

// As the Fortran: a configuration beyond NumSlabConfigs.
TEST(ScmPositionFileErrors, ConfigurationBeyondNumSlabConfigs) {
  expect_position_file_error(small_v3_file, 2, 16.0, "PosFileSlabNum", {"configuration 2", "NumSlabConfigs = 1"});
}

// SCM D6: a configuration number ≤ 0 (the Fortran reads configuration 1).
TEST(ScmPositionFileErrors, ScmD6ConfigurationAtLeast1) {
  for (int configuration : {0, -1}) {
    expect_input_error([&] { (void)read_slab_text(small_v3_file, configuration, 16.0); }, "PosFileSlabNum",
                       {"at least 1", std::to_string(configuration)});
  }
}

TEST(ScmPositionFileErrors, CutoffMustBePositive) {
  for (double cutoff : {0.0, -1.0}) {
    expect_input_error([&] { (void)read_slab_text(small_v3_file, 1, cutoff); }, "SCMCutoff", {"positive"});
  }
}

// SCM D14: nz ≥ 2 (as the Fortran) and nx, ny ≥ 1 in version 3 (two negative sizes would give
// the Fortran a positive atom count, here 2).
TEST(ScmPositionFileErrors, ScmD14NzAtLeast2NxNyAtLeast1) {
  expect_position_file_error("s FCC111\n 1 1 1 1 3 1\n 4.8\n 0.0 0.0 0.1 0.0 0.0 0.0 1\n", 1, 16.0, "PosFilepath",
                             {"nz = 1", "at least 2"});
  expect_position_file_error("s FCC111\n 0 1 2 1 3 2\n 4.8\n", 1, 16.0, "PosFilepath", {"nx = 0"});
  expect_position_file_error("s FCC111\n -1 -1 2 1 3 2\n 4.8\n", 1, 16.0, "PosFilepath", {"nx = -1, ny = -1"});
}

// SCM D6: MaxNumSurfAtoms ≥ 1 in version 2 (the Fortran silently uses no atom).
TEST(ScmPositionFileErrors, ScmD6MaxNumSurfAtomsAtLeast1) {
  expect_position_file_error("s FCC111\n 1 1 2 1 2 0\n 4.8\n", 1, 16.0, "PosFilepath", {"MaxNumSurfAtoms = 0"});
  expect_position_file_error("s FCC111\n 1 1 2 1 2 -3\n 4.8\n", 1, 16.0, "PosFilepath", {"MaxNumSurfAtoms = -3"});
}

// SCM D6: a layer index ≤ 0, in both versions, also for an atom outside the cut-off.
TEST(ScmPositionFileErrors, ScmD6LayerAtLeast1) {
  for (const std::string version : {"3", "2"}) {
    for (const std::string layer : {"0", "-1"}) {
      const std::string file = "s FCC111\n 1 1 2 1 " + version + " 2\n 4.8\n 0.0 0.0 0.1 0.0 0.0 0.0 1\n" +
                               " 1.0 1.0 -3.9 1.0 1.0 -4.0 " + layer + "\n";
      expect_position_file_error(file, 1, 16.0, "PosFilepath", {"atom 2 of configuration 1", "layer " + layer});
    }
  }
  const std::string outside =
      "s FCC111\n 1 1 2 1 3 2\n 4.8\n 0.0 0.0 0.1 0.0 0.0 0.0 1\n 50.0 0.0 0.0 50.0 0.0 0.0 0\n";
  expect_position_file_error(outside, 1, 16.0, "PosFilepath", {"atom 2 of configuration 1", "layer 0"});
}

TEST(ScmPositionFileErrors, ShortFileOrBadValue) {
  expect_position_file_error("s FCC111\n 1 1 2 1 3 2\n 4.8\n 0.0 0.0 0.1 0.0 0.0 0.0 1\n", 1, 16.0, "PosFilepath",
                             {"end of the file", "atom 2 of configuration 1"});
  // Configuration 3 of a file that ends inside configuration 2: lines 1-3 the header, then two
  // lines per configuration (4-5, 6-7, 8-9); the file has 6 lines.
  const std::string atom = " 0.0 0.0 0.1 0.0 0.0 0.0 1\n";
  expect_position_file_error("s FCC111\n 1 1 2 3 3 2\n 4.8\n" + atom + atom + atom, 3, 16.0, "PosFilepath",
                             {"ends after line 6", "configuration 3 would begin on line 8"});
  expect_position_file_error("s FCC111\n 1 1 seven 1 3 2\n 4.8\n", 1, 16.0, "PosFilepath",
                             {":2:", "'seven' is not an integer", "header"});
  expect_position_file_error("s FCC111\n 1 1 2 1 3 2\n 4.8\n 0.0 0.0 0.1 0.0 0.0 0.0 1.0\n", 1, 16.0, "PosFilepath",
                             {":4:", "'1.0' is not an integer", "atom 1"});
}

// As the Fortran ("No atoms in SCM range"): no atom of layers 1-2 inside the cut-off (version 3).
TEST(ScmPositionFileErrors, NoAtomInsideTheCutoff) {
  expect_position_file_error(small_v3_file, 1, 0.05, "PosFilepath", {"no atom", "cut-off"});
}

// SCM D6: a version 2 configuration without an atom of layers 1-2 (the Fortran silently gives ΔV = 0).
TEST(ScmPositionFileErrors, ScmD6Version2ConfigurationWithoutLayer1Or2Atom) {
  expect_position_file_error("s FCC111\n 1 1 2 1 2 1\n 4.8\n 0.0 0.0 0.1 0.0 0.0 0.0 3\n", 1, 16.0, "PosFilepath",
                             {"no atom", "layers 1–2"});
}

TEST(ScmPositionFileErrors, MissingFile) {
  const TemporaryDirectory directory("scm_io");
  const std::filesystem::path file = directory.path() / "missing.dat";
  expect_input_error([&] { (void)read_scm_slab(file, 1, 16.0); }, "PosFilepath", {file.string(), "no such file"});
}

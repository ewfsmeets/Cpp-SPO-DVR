// Tests: reading, writing and validating CRP data sets (src/crplib/dataset.*), on synthetic data
// sets (tests/crplib/crp_synthetic.hpp).
#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "crp_synthetic.hpp"
#include "crplib/dataset.hpp"
#include "crplib/errors.hpp"
#include "temporary_directory.hpp"

using namespace crplib;
using test_support::TemporaryDirectory;

namespace {

void expect_equal(const Table1D& a, const Table1D& b) {
  EXPECT_EQ(a.x, b.x);
  EXPECT_EQ(a.f, b.f);
}

// Field by field, exactly (the file image must survive a write/read round trip unchanged).
void expect_equal(const DataSet& a, const DataSet& b) {
  EXPECT_EQ(a.symmetry, b.symmetry);
  EXPECT_EQ(a.units, b.units);
  EXPECT_EQ(a.lattice_constant, b.lattice_constant);
  EXPECT_EQ(a.interlayer_spacings, b.interlayer_spacings);
  EXPECT_EQ(a.displacement_1, b.displacement_1);
  EXPECT_EQ(a.displacement_2, b.displacement_2);
  ASSERT_EQ(a.info.has_value(), b.info.has_value());
  if (a.info) {
    EXPECT_EQ(a.info->short_info, b.info->short_info);
    EXPECT_EQ(a.info->reference, b.info->reference);
  }
  expect_equal(a.corrugation, b.corrugation);
  EXPECT_EQ(a.corrugation_slope_first, b.corrugation_slope_first);
  EXPECT_EQ(a.corrugation_slope_last, b.corrugation_slope_last);
  ASSERT_EQ(a.atom_cuts.size(), b.atom_cuts.size());
  for (std::size_t k = 0; k < a.atom_cuts.size(); ++k) {
    EXPECT_EQ(a.atom_cuts[k].x, b.atom_cuts[k].x);
    EXPECT_EQ(a.atom_cuts[k].y, b.atom_cuts[k].y);
    EXPECT_EQ(a.atom_cuts[k].file_name, b.atom_cuts[k].file_name);
    EXPECT_EQ(a.atom_cuts[k].slope_first, b.atom_cuts[k].slope_first);
    EXPECT_EQ(a.atom_cuts[k].slope_last, b.atom_cuts[k].slope_last);
    expect_equal(a.atom_cuts[k].potential, b.atom_cuts[k].potential);
  }
  EXPECT_EQ(a.cuts_interpolation, b.cuts_interpolation);
  EXPECT_EQ(a.zso_min, b.zso_min);
  EXPECT_EQ(a.zso_max, b.zso_max);
  expect_equal(a.z_dependence, b.z_dependence);
  expect_equal(a.h2_potential, b.h2_potential);
  ASSERT_EQ(a.molecule_cuts.size(), b.molecule_cuts.size());
  for (std::size_t k = 0; k < a.molecule_cuts.size(); ++k) {
    const MoleculeCut &p = a.molecule_cuts[k], &q = b.molecule_cuts[k];
    EXPECT_EQ(p.present, q.present) << k;
    EXPECT_EQ(p.x, q.x);
    EXPECT_EQ(p.y, q.y);
    EXPECT_EQ(p.theta, q.theta);
    EXPECT_EQ(p.phi, q.phi);
    EXPECT_EQ(p.file_name, q.file_name);
    EXPECT_EQ(p.r, q.r);
    EXPECT_EQ(p.z, q.z);
    EXPECT_EQ(p.v, q.v);
  }
}

// read_data_set must fail with a message containing `fragment`.
void expect_data_error(const std::filesystem::path& directory, const std::string& fragment) {
  try {
    (void)read_data_set(directory);
    ADD_FAILURE() << "accepted; expected an error containing '" << fragment << "'";
  } catch (const DataError& e) {
    EXPECT_NE(std::string(e.what()).find(fragment), std::string::npos) << e.what();
  }
}

// Writes the synthetic C6v data set, modified by `edit`, and expects `fragment` in the error.
template <typename Edit>
void expect_rejected(Symmetry symmetry, Edit edit, const std::string& fragment) {
  TemporaryDirectory dir("dataset");
  DataSet d = crplib_test::pairwise_data_set(symmetry);
  edit(d);
  write_data_set(d, dir.path() / "set");
  expect_data_error(dir.path() / "set", fragment);
}

}  // namespace

class DataSetRoundTrip : public testing::TestWithParam<Symmetry> {};

TEST_P(DataSetRoundTrip, WriteAndReadBackExactly) {
  TemporaryDirectory dir("dataset");
  for (const DataSet& d : {crplib_test::pairwise_data_set(GetParam()),
                           crplib_test::in_angstrom_and_ev(crplib_test::pairwise_data_set(GetParam())),
                           crplib_test::impulse_data_set(GetParam(), 3)}) {
    write_data_set(d, dir.path() / "set");
    expect_equal(read_data_set(dir.path() / "set"), d);
    std::filesystem::remove_all(dir.path() / "set");
  }
}

INSTANTIATE_TEST_SUITE_P(Symmetries, DataSetRoundTrip, testing::Values(Symmetry::c3v, Symmetry::c4v, Symmetry::c6v));

TEST(DataSet, OptionalLinesAndFiles) {
  TemporaryDirectory dir("dataset");
  DataSet d = crplib_test::pairwise_data_set(Symmetry::c6v);
  d.info.reset();
  write_data_set(d, dir.path() / "set");
  // CRPData without the Disp line and with a single layer (an empty line for the spacings).
  dir.write("set/CRPData", "C6v\n0 2\n4.8\n1\n\n");
  const DataSet read = read_data_set(dir.path() / "set");
  EXPECT_FALSE(read.info.has_value());
  EXPECT_EQ(read.num_layers(), 1);
  EXPECT_EQ(read.displacement_1, 0.0);
  // An unreadable Disp line also gives 0 (Fortran: ERR=).
  dir.write("set/CRPData", "C6v\n0 2\n4.8\n2\n3.9\nabc def\n");
  EXPECT_EQ(read_data_set(dir.path() / "set").displacement_2, 0.0);
  dir.write("set/CRPData", "C6v\n0 2\n4.8\n2\n3.9\n0.1 -0.2\n");
  EXPECT_EQ(read_data_set(dir.path() / "set").displacement_2, -0.2);
  // Partial info file: the missing lines stay empty.
  dir.write("set/CRPData.info", "only the short description\n");
  EXPECT_EQ(read_data_set(dir.path() / "set").info->short_info, "only the short description");
}

TEST(DataSet, MoleculeCutRowsLastRowWins) {
  TemporaryDirectory dir("dataset");
  write_data_set(crplib_test::pairwise_data_set(Symmetry::c6v), dir.path() / "set");
  std::string text = "4\n4\n";
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      // r values of the first rows differ from the last ones; the Fortran keeps the last.
      const double r = 1.0 + i + (j < 3 ? 0.5 : 0.0), z = 2.0 * j + (i < 3 ? 0.25 : 0.0);
      text += std::format("{} {} {}\n", z, r, 10 * j + i);
    }
  }
  dir.write("set/Molecule/cut01.dat", text);
  const DataSet d = read_data_set(dir.path() / "set");
  EXPECT_EQ(d.molecule_cuts[0].r, (std::vector<double>{1.0, 2.0, 3.0, 4.0}));
  EXPECT_EQ(d.molecule_cuts[0].z, (std::vector<double>{0.0, 2.0, 4.0, 6.0}));
  EXPECT_EQ(d.molecule_cuts[0].v[1 + 4 * 2], 21.0);
}

TEST(DataSet, MissingDataSetAndUnsupportedHeader) {
  TemporaryDirectory dir("dataset");
  expect_data_error(dir.path() / "nothing", "does not exist");
  write_data_set(crplib_test::pairwise_data_set(Symmetry::c6v), dir.path() / "set");
  dir.write("set.dat", "1\n");
  expect_data_error(dir.path() / "set", "multi-PES");
}

TEST(DataSet, CrpDataErrors) {
  TemporaryDirectory dir("dataset");
  write_data_set(crplib_test::pairwise_data_set(Symmetry::c6v), dir.path() / "set");
  dir.write("set/CRPData", "c6v\n0 2\n4.8\n2\n3.9\n");
  expect_data_error(dir.path() / "set", "unknown symmetry");  // case-sensitive, as the Fortran
  dir.write("set/CRPData", "C6v\n2 2\n4.8\n2\n3.9\n");
  expect_data_error(dir.path() / "set", "unknown length unit");
  dir.write("set/CRPData", "C6v\n0 7\n4.8\n2\n3.9\n");
  expect_data_error(dir.path() / "set", "unknown energy unit");
  dir.write("set/CRPData", "C6v\n0 2\n0.0\n2\n3.9\n");
  expect_data_error(dir.path() / "set", "lattice constant");
  dir.write("set/CRPData", "C6v\n0 2\n4.8\n0\n");
  expect_data_error(dir.path() / "set", "number of layers");
}

// Unreadable text (the file reader's fortran_text::ReadError) reaches the caller as a DataError
// (expect_data_error catches nothing else), with the reader's message naming the file and line.
TEST(DataSet, ReaderErrorsAreDataErrors) {
  TemporaryDirectory dir("dataset");
  const auto fresh_set = [&](const std::string& name) {
    write_data_set(crplib_test::pairwise_data_set(Symmetry::c6v), dir.path() / name);
    return dir.path() / name;
  };
  // A unit code that is not an integer: the reader's message, with the file name put in front.
  const std::filesystem::path units = fresh_set("units");
  dir.write("units/CRPData", "C6v\nx 2\n4.8\n2\n3.9\n");
  const std::string crpdata = (units / "CRPData").string();
  expect_data_error(units, crpdata + ": " + crpdata + ":2: 'x' is not an integer");
  // A malformed number in an atom-surface file (the second value of the second point).
  const std::filesystem::path corrugation = fresh_set("corrugation");
  dir.write("corrugation/Atom/intrep.dat", "-0.036\n0\n4\n0 0.03\n0.25 abc\n");
  expect_data_error(corrugation, (corrugation / "Atom" / "intrep.dat").string() + ":5: 'abc' is not a real number");
  // A missing file of the molecule-surface PES.
  const std::filesystem::path h2 = fresh_set("h2");
  std::filesystem::remove(h2 / "Molecule" / "VH2.dat");
  expect_data_error(h2, (h2 / "Molecule" / "VH2.dat").string() + ": cannot open the file");
}

TEST(DataSet, LayerRequirements) {
  expect_rejected(Symmetry::c3v, [](DataSet& d) { d.interlayer_spacings.pop_back(); }, "C3v needs at least 4 layers");
  expect_rejected(Symmetry::c3v, [](DataSet& d) { d.displacement_1 = 0.1; }, "displacements");
  expect_rejected(Symmetry::c4v, [](DataSet& d) { d.interlayer_spacings.clear(); }, "C4v needs at least 2 layers");
}

TEST(DataSet, CutCounts) {
  expect_rejected(Symmetry::c6v, [](DataSet& d) { d.atom_cuts.pop_back(); }, "exactly 7 atom-surface cuts");
  expect_rejected(Symmetry::c6v, [](DataSet& d) { d.atom_cuts.push_back(d.atom_cuts.back()); }, "exactly 7");
  expect_rejected(Symmetry::c3v, [](DataSet& d) { d.atom_cuts.pop_back(); }, "exactly 10");
  expect_rejected(Symmetry::c4v, [](DataSet& d) { d.atom_cuts.pop_back(); }, "at least 6");
  expect_rejected(Symmetry::c6v, [](DataSet& d) { d.molecule_cuts.pop_back(); }, "at least 18 molecule-surface cuts");
  expect_rejected(Symmetry::c3v, [](DataSet& d) { d.molecule_cuts.pop_back(); }, "at least 29");
  expect_rejected(Symmetry::c4v, [](DataSet& d) { d.molecule_cuts.pop_back(); }, "at least 16");
}

TEST(DataSet, CutFiles) {
  TemporaryDirectory dir("dataset");
  write_data_set(crplib_test::pairwise_data_set(Symmetry::c6v), dir.path() / "set");
  std::filesystem::remove(dir.path() / "set/Molecule/cut05.dat");
  expect_data_error(dir.path() / "set", "molecule-surface cut file not found");
  write_data_set(crplib_test::pairwise_data_set(Symmetry::c6v), dir.path() / "set");
  std::filesystem::remove(dir.path() / "set/Atom/site3.dat");
  expect_data_error(dir.path() / "set", "atom-surface cut file not found");
  expect_rejected(
      Symmetry::c6v, [](DataSet& d) { d.molecule_cuts[2].file_name = std::string(51, 'x'); },
      "longer than 50 characters");
}

TEST(DataSet, SplineGrids) {
  expect_rejected(
      Symmetry::c6v, [](DataSet& d) { d.h2_potential.x.resize(3), d.h2_potential.f.resize(3); }, "at least 4 points");
  expect_rejected(
      Symmetry::c6v, [](DataSet& d) { std::swap(d.corrugation.x[2], d.corrugation.x[3]); }, "not strictly increasing");
  expect_rejected(
      Symmetry::c6v, [](DataSet& d) { d.molecule_cuts[0].z[4] = d.molecule_cuts[0].z[3]; }, "not strictly increasing");
  expect_rejected(Symmetry::c4v, [](DataSet& d) { d.atom_cuts[3].potential.x[5] += 1e-3; }, "same z values");
  expect_rejected(Symmetry::c6v, [](DataSet& d) { std::swap(d.zso_min, d.zso_max); }, "ZsoMin is larger");
}

TEST(DataSet, EqualSwitchingLimitsAreValid) {
  TemporaryDirectory dir("dataset");
  DataSet d = crplib_test::pairwise_data_set(Symmetry::c6v);
  d.zso_max = d.zso_min;  // the switch function returns before it divides
  write_data_set(d, dir.path() / "set");
  EXPECT_NO_THROW((void)read_data_set(dir.path() / "set"));
}

// An unused extra cut may lack its file: it is skipped, as in the Fortran. A used cut may not.
TEST(DataSet, MissingFileOfAnUnusedExtraCut) {
  TemporaryDirectory dir("dataset");
  const DataSet d = crplib_test::case_data_set("c3v_missing_extra_cut");
  write_data_set(d, dir.path() / "set");
  EXPECT_TRUE(std::filesystem::exists(dir.path() / "set/Molecule/cut30.dat"));
  EXPECT_FALSE(std::filesystem::exists(dir.path() / "set/Molecule/cut31.dat"));
  const DataSet back = read_data_set(dir.path() / "set");
  ASSERT_EQ(back.molecule_cuts.size(), 31u);
  EXPECT_TRUE(back.molecule_cuts[29].present);
  EXPECT_FALSE(back.molecule_cuts[30].present);
  expect_equal(back, d);
  const auto warnings = data_set_warnings(back);
  ASSERT_EQ(warnings.size(), 2u);
  EXPECT_NE(warnings[0].find("only the first 29"), std::string::npos);
  EXPECT_NE(warnings[1].find("no file for the unused cuts cut31.dat"), std::string::npos);

  DataSet used = d;  // a used cut without a file, also in a data set built in memory
  used.molecule_cuts[3].present = false;
  EXPECT_THROW(validate_data_set(used), DataError);
  std::filesystem::remove(dir.path() / "set/Molecule/cut04.dat");
  expect_data_error(dir.path() / "set", "molecule-surface cut file not found");
}

TEST(DataSet, Warnings) {
  DataSet d = crplib_test::pairwise_data_set(Symmetry::c6v);
  EXPECT_TRUE(data_set_warnings(d).empty());
  d.molecule_cuts.push_back(d.molecule_cuts.back());
  d.z_dependence.x.back() = d.zso_max - 0.1;
  const auto warnings = data_set_warnings(d);
  ASSERT_EQ(warnings.size(), 2u);
  EXPECT_NE(warnings[0].find("only the first 18"), std::string::npos);
  EXPECT_NE(warnings[1].find("below ZsoMax"), std::string::npos);
}

// --- the committed reference data sets -----------------------------------------------------------

namespace {

std::string file_content(const std::filesystem::path& file) {
  std::ifstream in(file, std::ios::binary);
  return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

// The files of a directory tree, relative to it, sorted.
std::vector<std::filesystem::path> files_below(const std::filesystem::path& dir) {
  std::vector<std::filesystem::path> files;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(dir))
    if (entry.is_regular_file()) files.push_back(std::filesystem::relative(entry.path(), dir));
  std::sort(files.begin(), files.end());
  return files;
}

}  // namespace

// Every data set in tests/crplib/reference/<case>/dataset (read by the Fortran reference runs and
// by the gate cases of tests/reference/crp_*) is byte for byte what the generator writes today.
TEST(ReferenceDataSets, MatchTheGenerator) {
  int cases = 0;
  for (const auto& entry : std::filesystem::directory_iterator(CRPLIB_REFERENCE_DIR)) {
    const std::filesystem::path committed = entry.path() / "dataset";
    if (!std::filesystem::is_directory(committed)) continue;
    const std::string name = entry.path().filename().string();
    TemporaryDirectory dir("reference_" + name);
    write_data_set(crplib_test::case_data_set(name), dir.path());
    const std::vector<std::filesystem::path> files = files_below(committed);
    ASSERT_EQ(files, files_below(dir.path())) << name;
    for (const std::filesystem::path& file : files)
      EXPECT_EQ(file_content(committed / file), file_content(dir.path() / file)) << name << ": " << file.string();
    ++cases;
  }
  EXPECT_EQ(cases, 20);
}

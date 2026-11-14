// Tests: the CRP PES in the program (src/potential/crp_potential.*, SystemID CRPlibrary): factory
// and system parameters, V equal to the library's QD wrapper, library errors as InputError naming
// CRPlibDataSet, the data-set summary and range report, the setup report, the surface-symmetry
// self-check on the synthetic data sets, and the potential grid (thread count, deterministic error).
// Data sets: the committed synthetic ones (tests/crplib/reference/<case>/dataset) and variants
// written by the library's test helpers (tests/crplib/crp_synthetic.hpp).
#include <gtest/gtest.h>
#include <omp.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include <nlohmann/json.hpp>

#include "comparison.hpp"
#include "core/constants.hpp"
#include "core/error.hpp"
#include "crp_synthetic.hpp"
#include "crplib/dataset.hpp"
#include "crplib/qd_wrapper.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "potential/crp_potential.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "report/crp_report.hpp"
#include "symmetry/grid_symmetry.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using crplib::Symmetry;
using spodvr::testing::nan_aware_max;
using test_support::TemporaryDirectory;
using Json = nlohmann::json;
namespace fs = std::filesystem;

namespace {

fs::path repository_root() { return fs::path(SPODVR_REFERENCE_DIR).parent_path().parent_path(); }

std::string lower_name(Symmetry symmetry) {
  return symmetry == Symmetry::c3v ? "c3v" : symmetry == Symmetry::c4v ? "c4v" : "c6v";
}

// The committed synthetic data set of the pairwise model.
fs::path pairwise_data_set(Symmetry symmetry) {
  return repository_root() / "tests/crplib/reference" / (lower_name(symmetry) + "_pairwise") / "dataset";
}

InputParameters crp_input(const fs::path& data_set) {
  InputParameters p;
  p.system_id = "CRPlibrary";
  p.crp_dataset = data_set.string();
  return p;
}

// The gate case tests/reference/crp_<sym>_synthetic_small with the data set given absolutely.
Json gate_case_json(Symmetry symmetry) {
  std::ifstream file(fs::path(SPODVR_REFERENCE_DIR) / ("crp_" + lower_name(symmetry) + "_synthetic_small") /
                     "input.json");
  Json json = Json::parse(file);
  json["CRPlibDataSet"] = (repository_root() / json["CRPlibDataSet"].get<std::string>()).string();
  return json;
}

// The C6v gate case on a smaller grid (fast enough for single-threaded runs) with another data set.
Json small_json(const fs::path& data_set) {
  Json json = gate_case_json(Symmetry::c6v);
  json["CRPlibDataSet"] = data_set.string();
  json["dZ"] = 0.25;
  json["nZ"] = 56;
  json["nZSpecular"] = 90;
  json["ZOptStart"] = 9.0;
  json["ZOptEnd"] = 11.75;
  json["ZSpOptStart"] = 16.75;
  json["ZSpOptEnd"] = 20.25;
  json["dR"] = 0.5;
  json["nR"] = 16;
  json["ROptStart"] = 4.3;
  json["ROptEnd"] = 8.3;
  json["nX"] = 4;
  json["nY"] = 4;
  json["JMax"] = 2;
  json["MJMax"] = 2;
  return json;
}

struct RunSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  explicit RunSetup(const Json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)) {}
};

// An impulse data set (atom-surface PES zero, molecular cut 1 equal to 1) whose atom cuts start at
// z = 0: a molecule near Zmin with a large tilt puts an atom below them, where the library stops.
fs::path write_shallow_atom_data_set(const fs::path& directory) {
  crplib_test::SyntheticGrids grids = crplib_test::impulse_grids();
  grids.atom_z = crplib_test::grid(0.0, 20.0, 5);
  crplib::write_data_set(crplib_test::impulse_data_set(Symmetry::c6v, 0, grids), directory);
  return directory;
}

std::string error_message(const std::function<void()>& f) {
  try {
    f();
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), "CRPlibDataSet") << e.what();
    return e.what();
  }
  ADD_FAILURE() << "expected an InputError";
  return "";
}

}  // namespace

class CrpPotentialTest : public ::testing::TestWithParam<Symmetry> {};

TEST_P(CrpPotentialTest, FactoryAndSystemParameters) {
  const Symmetry symmetry = GetParam();
  const std::unique_ptr<Potential> pes = make_potential(crp_input(pairwise_data_set(symmetry)));
  EXPECT_EQ(pes->name(), "CRPlibrary");
  const SystemParameters& p = pes->system_parameters();
  EXPECT_EQ(p.lattice_constant_x, crplib_test::lattice_constant(symmetry));  // bohr data
  EXPECT_EQ(p.lattice_constant_y, p.lattice_constant_x);
  EXPECT_EQ(p.z_asymptotic, 10.0);  // last Z of Molecule/extrapol.dat
  // CRP_Constants_PI / 3.0 or / 2.0: the same double as the program's pi
  EXPECT_EQ(p.skewing_angle, symmetry == Symmetry::c4v ? Constants::pi / 2.0 : Constants::pi / 3.0);
  EXPECT_EQ(p.point_group, crplib::symmetry_name(symmetry));
}

// V is the library's QD wrapper value, bitwise (the adapter adds nothing).
TEST_P(CrpPotentialTest, ValueIsTheLibraryValue) {
  const fs::path data_set = pairwise_data_set(GetParam());
  const CrpPotential pes(data_set);
  const crplib::QdWrapper wrapper(crplib::read_data_set(data_set));
  for (double z : {-1.0, 0.7, 2.3, 5.6, 7.1, 12.0})
    for (double r : {0.3, 1.4, 2.9, 6.0})
      for (double theta : {0.0, 0.8, Constants::pi / 2.0}) {
        const double x = 0.37 * z, y = 1.9 - 0.2 * r, phi = 2.1 * theta - 0.4;
        EXPECT_EQ(pes.value(x, y, z, r, theta, phi, 0, 0), wrapper.potential_ev(x, y, z, r, theta, phi))
            << z << ' ' << r << ' ' << theta;
      }
}

// The synthetic data sets interpolate an exactly symmetric model with a symmetric scheme: the PES
// self-check (UseSymmetry) passes with the full point group (nPhisOptimal) at rounding level.
TEST_P(CrpPotentialTest, SymmetryCheckPasses) {
  Json json = gate_case_json(GetParam());
  json["UseSymmetry"] = true;
  json["nPhisOptimal"] = true;
  const RunSetup s(json);
  const std::unique_ptr<Potential> pes = make_potential(s.config.input);
  const GridSymmetry symmetry = make_grid_symmetry(s.config, s.angular);
  EXPECT_EQ(symmetry.group_order(), GetParam() == Symmetry::c6v ? 12 : GetParam() == Symmetry::c4v ? 8 : 6);
  EXPECT_TRUE(symmetry.dropped_operations().empty());
  const SymmetryCheck check = check_potential_symmetry(*pes, symmetry, s.config, s.grids, s.angular);
  EXPECT_EQ(check.n_comparisons, check.n_samples * (symmetry.group_order() - 1));
  EXPECT_LT(check.max_deviation, 1e-11) << "eV";
  std::cout << "  max |V(q) - V(g q)| = " << check.max_deviation << " eV\n";
  EXPECT_TRUE(check.passed());
}

// Speed-up S1 (docs/crp_library.md): with UseSymmetry the PES is evaluated at the (x, y)
// representatives only. Gathered through the symmetry maps, the reduced grid equals the full grid
// bitwise at the representatives and to the PES's asymmetry (rounding) at their images.
TEST_P(CrpPotentialTest, SymmetricGridGathersTheFullGrid) {
  Json json = gate_case_json(GetParam());
  json["UseSymmetry"] = true;
  json["nPhisOptimal"] = true;
  const RunSetup s(json);
  const std::unique_ptr<Potential> pes = make_potential(s.config.input);
  const PotentialGrid reduced = compute_potential_grid(*pes, s.config, s.grids, s.angular);
  const PotentialGrid full = compute_potential_grid(
      *pes, s.config, s.grids, s.angular, GridSymmetry(s.config.input.nx, s.config.input.ny, s.angular.n_phi()));
  const GridSymmetry& g = reduced.symmetry;
  std::cout << "  (x, y) points computed: " << g.n_representatives() << " of " << g.n_xy() << '\n';
  EXPECT_LE(g.n_representatives() * 3, g.n_xy());
  EXPECT_EQ(reduced.asymptotic, full.asymptotic);
  const double one_ev = 1.0 / 27.2114;
  double worst = 0.0;
  for (std::size_t izr = 0; izr < full.index.n2; ++izr)
    for (int ixy = 0; ixy < g.n_xy(); ++ixy) {
      const bool representative = g.representative_xy(g.representative_of(ixy)) == ixy;
      for (std::size_t itp = 0; itp < full.index.n0; ++itp) {
        const double v = reduced(itp, ixy, izr), v_full = full(itp, ixy, izr);
        if (representative) {
          ASSERT_EQ(v, v_full) << ixy;
        }
        worst = nan_aware_max(worst, std::abs(v - v_full) / std::max(std::abs(v_full), one_ev));
      }
    }
  EXPECT_LT(worst, 5e-13);  // rounding: the images run through other fold and Fourier arithmetic
  std::cout << "  max scaled difference reduced vs full: " << worst << '\n';
}

INSTANTIATE_TEST_SUITE_P(Symmetries, CrpPotentialTest, ::testing::Values(Symmetry::c3v, Symmetry::c4v, Symmetry::c6v),
                         [](const ::testing::TestParamInfo<Symmetry>& param) { return lower_name(param.param); });

TEST(CrpPotential, ImplementedCaseInsensitive) {
  EXPECT_TRUE(potential_is_implemented("CRPlibrary"));
  EXPECT_TRUE(potential_is_implemented("crplibrary"));
}

TEST(CrpPotential, Summary) {
  const CrpDataSetSummary s = CrpPotential(pairwise_data_set(Symmetry::c3v)).summary();
  EXPECT_EQ(s.description, "synthetic test data set");
  EXPECT_EQ(s.symmetry, "C3v");
  EXPECT_EQ(s.file_units, "bohr, hartree");
  EXPECT_EQ(s.surface_layers, 4);
  EXPECT_EQ(s.atom_cuts, 10);
  EXPECT_EQ(s.molecule_cuts, 29);
  EXPECT_EQ(s.molecule_cuts_used, 29);
  EXPECT_FALSE(s.cuts_hold_interpolation);
  EXPECT_EQ(s.r_min, 0.6);
  EXPECT_EQ(s.r_max, 3.2);
  EXPECT_EQ(s.z_min, 0.5);
  EXPECT_EQ(s.zso_min, 5.0);
  EXPECT_EQ(s.zso_max, 6.5);
  EXPECT_EQ(s.z_asymptotic, 10.0);
  EXPECT_EQ(s.energy_shift, -27.2113961 * crplib_test::PairwiseModel::h2(8.0));  // -VH2(r_last) in eV
  EXPECT_TRUE(s.warnings.empty());
  EXPECT_EQ(s.cluster_atoms, 256);  // C3v MURII: 8 x 8 x 4
  EXPECT_EQ(s.cluster_atoms_within_reach, 54);

  const fs::path angstrom = repository_root() / "tests/crplib/reference/c6v_pairwise_angstrom/dataset";
  EXPECT_EQ(CrpPotential(angstrom).summary().file_units, "angstrom, eV");
}

// The C6v gate case: r = 0.8 + 0.25 k has 10 points in [Rmin, Rmax] = [0.6, 3.2] (22 clamped);
// Z = -2 + 0.125 k below Zmin = 0.5 (+1e-12) for k <= 20 (21 clamped of the 80 6D points).
TEST(CrpPotential, RangeReportAndSetupReport) {
  const RunSetup s(gate_case_json(Symmetry::c6v));
  const CrpPotential pes(s.config.input.crp_dataset);
  const auto z_6d = s.grids.z_scattering.points().first(s.grids.nz_potential);
  const CrpRangeReport range = pes.range_report(z_6d, s.grids.r.points());
  EXPECT_EQ(range.r_points, 32);
  EXPECT_EQ(range.r_clamped, 22);
  EXPECT_EQ(range.z_points, 80);
  EXPECT_EQ(range.z_clamped, 21);
  EXPECT_EQ(range.zr_points, 2560);
  EXPECT_EQ(range.zr_uncovered, 0);

  std::ostringstream out;
  print_crp_data_set(out, Style{}, pes.summary(), range);
  const std::string text = out.str();
  for (const char* expected : {"CRP DATA SET", "C6v", "18 (18)", "22 of 32", "21 of 80", "0 of 2560"})
    EXPECT_NE(text.find(expected), std::string::npos) << expected << " in\n" << text;

  // Analysis-only runs: no 6D points, only the r line.
  std::ostringstream analysis_only;
  print_crp_data_set(analysis_only, Style{}, pes.summary(), pes.range_report({}, s.grids.r.points()));
  EXPECT_EQ(analysis_only.str().find("(6D part)"), std::string::npos);
}

TEST(CrpPotential, UnusableDataSetNamesTheKey) {
  EXPECT_NE(error_message([] { (void)make_potential(crp_input("no/such/data/set")); }).find("does not exist"),
            std::string::npos);

  TemporaryDirectory dir("crp_potential");
  crplib::write_data_set(crplib_test::pairwise_data_set(Symmetry::c6v), dir.path());
  std::ofstream(dir.path() / "CRPData") << "C5v\n";
  (void)error_message([&] { (void)CrpPotential(dir.path()); });

  // Unreadable text (a reader error in the library) too: a lattice constant that is not a number.
  std::ofstream(dir.path() / "CRPData") << "C6v\n0 2\nabc\n2\n3.9\n";
  const std::string message = error_message([&] { (void)CrpPotential(dir.path()); });
  EXPECT_NE(message.find((dir.path() / "CRPData").string() + ":3: 'abc' is not a real number"), std::string::npos)
      << message;
}

// Where the Fortran library stops (an atom below the atom-surface cuts) the program reports an
// input error naming CRPlibDataSet and the grid point.
TEST(CrpPotential, EvaluationErrorNamesTheKeyAndThePoint) {
  TemporaryDirectory dir("crp_potential");
  const CrpPotential pes(write_shallow_atom_data_set(dir.path()));
  EXPECT_NO_THROW((void)pes.value(0.0, 0.0, 2.0, 1.4, 0.0, 0.0, 0, 0));
  const std::string message = error_message([&] { (void)pes.value(0.0, 0.0, 0.5, 3.0, 0.0, 0.0, 0, 0); });
  EXPECT_NE(message.find("cannot be evaluated at x = 0, y = 0, Z = 0.5, r = 3"), std::string::npos) << message;
}

// The library is const and allocation-free: the potential grid does not depend on the number of
// threads.
TEST(CrpPotentialGrid, IndependentOfThreadCount) {
  const RunSetup s(small_json(pairwise_data_set(Symmetry::c6v)));
  const CrpPotential pes(s.config.input.crp_dataset);
  omp_set_num_threads(1);
  const PotentialGrid serial = compute_potential_grid(pes, s.config, s.grids, s.angular);
  omp_set_num_threads(4);
  const PotentialGrid parallel = compute_potential_grid(pes, s.config, s.grids, s.angular);
  ASSERT_EQ(serial.values.size(), parallel.values.size());
  EXPECT_TRUE(std::equal(serial.values.begin(), serial.values.end(), parallel.values.begin()));
}

// An exception inside the parallel loop reaches the caller; the reported point is that of the
// lowest failing (Z, r) group, whatever the number of threads.
TEST(CrpPotentialGrid, ErrorIsDeterministic) {
  TemporaryDirectory dir("crp_potential");
  const RunSetup s(small_json(write_shallow_atom_data_set(dir.path())));
  const CrpPotential pes(s.config.input.crp_dataset);
  omp_set_num_threads(1);
  const std::string serial = error_message([&] { (void)compute_potential_grid(pes, s.config, s.grids, s.angular); });
  omp_set_num_threads(4);
  const std::string parallel = error_message([&] { (void)compute_potential_grid(pes, s.config, s.grids, s.angular); });
  // The first failing point (group 0): Z = -2 (evaluated at Zmin = 0.5), r = 1.3, the first theta.
  EXPECT_NE(serial.find("x = 0, y = 0, Z = -2, r = 1.3 bohr"), std::string::npos) << serial;
  EXPECT_EQ(serial, parallel);
}

// Regression gates 3 and 4 (docs/validation.md): the propagation and the on-the-fly analysis of
// tests/reference/hardwall_short (1000 steps) against the Fortran output of the same case
// (fortran_stdout.txt, fortran_averages.out, fortran_fort25.txt and, if present, fortran_fort21.bin, written by
// tools/fortran_reference/run_propagation_reference.sh). Skipped when the stdout or averages file is missing; a
// missing fortran_fort25.txt fails.
// The full 7000-step cases (hardwall_small, cu111_dailight_small) run only with
// SPODVR_LONG_TESTS=1. Each case skips while its Fortran reference files are missing, except the
// SCM case scm_cu111_srp48_small and cu111_dailight_oblique, which fail.
// LowZGateFailsWithoutTheScatteringAbsorberLong and ScmGateFailsWithoutTheCorrectionLong (also
// SPODVR_LONG_TESTS=1) are the negative tests of the low-Z and the SCM case. RegressionHarness: the data-file paths of
// a case's input (reference_case_input.hpp).
#include <cstdlib>
#include <fstream>
#include <string>

#include "potential/potential.hpp"
#include "potential/scm/scm_crp_potential.hpp"
#include "regression_gate.hpp"
#include "scm_example_input.hpp"

using namespace spodvr;
using namespace spodvr::testing;
namespace fs = std::filesystem;

TEST(RegressionPropagation, Gate3ShortRun) { run_gate("hardwall_short"); }

TEST(RegressionPropagation, Gate3FullRunLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  run_gate("hardwall_small");
}

// The corrugated Dai-Light H2/Cu(111) cases (60 deg cell, C6v): the first gates with theta/phi/x/y
// coupling. _sym: C6v surface symmetry in the potential (UseSymmetry, nPhisOptimal) in both codes.
TEST(RegressionPropagation, Gate3DaiLightShort) { run_gate("cu111_dailight_short"); }
TEST(RegressionPropagation, Gate3DaiLightShortSymmetric) { run_gate("cu111_dailight_short_sym"); }
TEST(RegressionPropagation, Gate3DaiLightFullRunLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  run_gate("cu111_dailight_small");
}

// The low-Z absorbers of upstream 5.24: cu111_dailight_short with ZOptLow… (energy form) and
// ZSpOptLow… (factor form), both from ZStart to Zana.
TEST(RegressionPropagation, Gate3DaiLightShortLowZ) { run_gate("cu111_dailight_short_lowz"); }

// Oblique incidence with J > 0: cu111_dailight_short with JIn 1, MJIn -1 (the odd-parity angular grid, mJ != 0),
// JMax 3, EParallel 0.1 eV and PhiParallel 17 deg (the momentum shift, the oblique channels); 1050 steps.
TEST(RegressionPropagation, Gate3DaiLightOblique) {
  run_gate("cu111_dailight_oblique", make_cpu_backend, References::required);
}

// Negative test of that gate: the C++ run without the scattering low-Z absorber must fail gates 3
// and 4 by at least 100 times their tolerance, so the gate sees that absorber. It is the weaker of
// the two (removing the specular one changes the averages about 100 times more: docs/validation.md).
// Long suite: it cannot catch a code regression the gate misses (a broken absorber fails the gate),
// its outcome changes only with the input or the reference, and a shorter run does not pay: the
// deviation passes 100 tolerances in both gates only from t = 1500 of 2000 (measured 2026-11-12).
TEST(RegressionPropagation, LowZGateFailsWithoutTheScatteringAbsorberLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  const fs::path dir = fs::path(SPODVR_REFERENCE_DIR) / "cu111_dailight_short_lowz";
  const auto ref = read_fortran_reference(dir / "fortran_stdout.txt", dir / "fortran_averages.out");
  if (!ref || !fs::exists(dir / "fortran_fort21.bin"))
    GTEST_SKIP() << "no Fortran reference in " << dir << "; run "
                 << "tools/fortran_reference/run_propagation_reference.sh cu111_dailight_short_lowz";
  std::ifstream file(dir / "input.json");
  nlohmann::json input = nlohmann::json::parse(file);
  for (const std::string key : {"ZOptLowStart", "ZOptLowEnd", "ZOptLowEnergy"}) ASSERT_EQ(input.erase(key), 1u) << key;
  const PropagationRun cpp = run_cpp_propagation(input);
  const auto bk_ref = read_fortran_bk_records(dir / "fortran_fort21.bin");
  ASSERT_EQ(cpp.averages.size(), ref->averages.size());
  ASSERT_EQ(cpp.bk.size(), bk_ref.size());

  // In the units of the gates: printed half units (gate 3 allows 1) and the BK tolerance (gate 4).
  const double averages = averages_deviation(cpp, *ref);
  const BkDeviation bk = bk_deviation(cpp, bk_ref);
  const double bk_in_tolerances = bk.worst / (bk_gate_tolerance * bk.scale);
  EXPECT_GT(averages, 100.0);
  EXPECT_GT(bk_in_tolerances, 100.0);
  std::cout << "  without the scattering low-Z absorber: averages " << averages << " x (half printed unit), BK "
            << bk_in_tolerances << " x the gate-4 tolerance\n";
}

// The SCM distorted-surface PES (SystemID CRPDistrtSurface, upstream 5.24): the SRP48 CRP data set
// plus the correction of the example slab, rescaled to the CRP lattice (so the x, y grid, the
// diffraction channels and the initial wave packet are those of the slab's lattice constant). The
// first CRP-based case of gates 3 and 4; 3000 steps, so that the wave reflected by the surface
// reaches the analysis line (why 3000: docs/validation.md, gate 3, SCM). A missing reference file
// fails it.
TEST(RegressionPropagation, Gate3ScmShort) {
  run_gate("scm_cu111_srp48_small", make_cpu_backend, References::required);
}

// The SCM potential of a case without its correction: the CRP at (u/s, v/s) on the slab's lattice,
// as ScmCrpPotential::value computes its CRP term (test-only, for the negative test below).
class ScmCrpTermOnly final : public Potential {
 public:
  explicit ScmCrpTermOnly(const Config& config)
      : scm_(config.input.crp_dataset, config.input.scm), crp_scale_(scm_.rescaling_factor().value_or(1.0)) {}
  [[nodiscard]] double value(double u, double v, double z, double r, double theta, double phi, int iz,
                             int ir) const override {
    return scm_.crp().value(u / crp_scale_, v / crp_scale_, z, r, theta, phi, iz, ir);
  }
  [[nodiscard]] const SystemParameters& system_parameters() const override { return scm_.system_parameters(); }
  [[nodiscard]] std::string name() const override { return scm_.name(); }

 private:
  ScmCrpPotential scm_;
  double crp_scale_;  // s = a_T/a0 with rescaling, else 1
};

// Negative test of the SCM gate: the C++ run with the CRP term alone (the correction removed, the
// lattice kept) must fail gates 3 and 4 by at least 100 times their tolerance, so the gates see the
// correction. Long suite, as the low-Z negative test: a broken correction fails the gate itself.
TEST(RegressionPropagation, ScmGateFailsWithoutTheCorrectionLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  const fs::path dir = fs::path(SPODVR_REFERENCE_DIR) / "scm_cu111_srp48_small";
  const auto ref = read_fortran_reference(dir / "fortran_stdout.txt", dir / "fortran_averages.out");
  ASSERT_TRUE(ref.has_value()) << "no Fortran reference in " << dir;
  ASSERT_TRUE(fs::exists(dir / "fortran_fort21.bin")) << "no fortran_fort21.bin in " << dir;
  std::ifstream file(dir / "input.json");
  const Config config = case_config(nlohmann::json::parse(file));
  const PropagationRun cpp = run_cpp_propagation(config, ScmCrpTermOnly(config));
  const auto bk_ref = read_fortran_bk_records(dir / "fortran_fort21.bin");
  ASSERT_EQ(cpp.averages.size(), ref->averages.size());
  ASSERT_EQ(cpp.bk.size(), bk_ref.size());

  // In the units of the gates: printed half units (gate 3 allows 1) and the BK tolerance (gate 4).
  const double averages = averages_deviation(cpp, *ref);
  const BkDeviation bk = bk_deviation(cpp, bk_ref);
  const double bk_in_tolerances = bk.worst / (bk_gate_tolerance * bk.scale);
  EXPECT_GT(averages, 100.0);
  EXPECT_GT(bk_in_tolerances, 100.0);
  std::cout << "  without the SCM correction: averages " << averages << " x (half printed unit), BK "
            << bk_in_tolerances << " x the gate-4 tolerance\n";
}

TEST(RegressionPropagation, ParsesTheFortranOutputFormats) {
  // The relevant lines in the formats of WaveFunction.f90 CoherentStateFunc (stdout) and
  // Averages.f90 (averages.out). Five-column lines in the stdout must not count as averages.
  const fs::path stdout_file = fs::temp_directory_path() / "spodvr_test_fortran_stdout.txt";
  const fs::path averages_file = fs::temp_directory_path() / "spodvr_test_fortran_averages.out";
  {
    std::ofstream out(stdout_file);
    out << " * Norm initial wave packet (Z,r) =  0.10000000E+01\n"
        << "   - at k = 0                                              0.48657986E-06\n"
        << "   - at the end/beginning of the grid                      0.79631124E-16\n"
        << "   - maximum value =                                       0.54436230E+00\n"
        << "   - Low boundary (asymptotic value of Z), iLow =                                               126\n"
        << "   - High boundary (start of opt. pot on spec. grid OR a specifically set ZAnaSp ), iHigh =     235\n"
        << " * Norm of the initial wave packet between iLow and iHigh  =  0.10000000E+01\n"
        << "  0    7     0    0  0.4317480E-32\n";
    std::ofstream avg(averages_file);
    avg << "********************************************************************************\n"
        << "               INFO ON PROPAGATING THE WAVE PACKET\n"
        << "********************************************************************************\n\n"
        << "      time             norm           z              r        rest norm\n\n"
        << "  0.00000000E+00  0.10000000E+01  0.11440000E+02  0.14708348E+01  0.10000000E+01\n\n"
        << "  0.10000000E+03  0.99999946E+00  0.11166768E+02  0.14708477E+01  0.99999945E+00\n";
  }
  const auto ref = read_fortran_reference(stdout_file, averages_file);
  fs::remove(stdout_file);
  fs::remove(averages_file);
  ASSERT_TRUE(ref.has_value());
  EXPECT_EQ(ref->norm, 1.0);
  EXPECT_EQ(ref->k_at_zero, 0.48657986E-06);
  EXPECT_EQ(ref->k_maximum, 0.54436230E+00);
  EXPECT_EQ(ref->i_low, 126);
  EXPECT_EQ(ref->i_high, 235);
  ASSERT_EQ(ref->averages.size(), 2u);
  EXPECT_EQ(ref->averages[1].time, 100.0);
  EXPECT_EQ(ref->averages[1].rest_norm, 0.99999945);
}

// --- the data-file paths of a case's input.json ---------------------------------------------------

// They are relative to the repository root; the key absent (empty) and an absolute path stay.
TEST(RegressionHarness, DataPathsResolveAgainstTheRepositoryRoot) {
  InputParameters input;
  resolve_repository_paths(input);
  EXPECT_TRUE(input.crp_dataset.empty());
  EXPECT_TRUE(input.scm.parameter_file.empty());
  EXPECT_TRUE(input.scm.position_file.empty());
  input.crp_dataset = "potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v";
  input.scm.parameter_file = "potentials/H2Cu111SRP_C3v.dir/params.dat";
  input.scm.position_file = "/an/absolute/slab.dat";
  resolve_repository_paths(input);
  EXPECT_EQ(input.crp_dataset, (repository_root() / "potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v").string());
  EXPECT_EQ(input.scm.parameter_file, (repository_root() / "potentials/H2Cu111SRP_C3v.dir/params.dat").string());
  EXPECT_EQ(input.scm.position_file, "/an/absolute/slab.dat");
}

// An SCM case whose input.json names the CRP data set and both SCM files relative to the repository
// root (the SCM example) builds through the harness of the gates (case_setup, used by
// run_cpp_propagation and so by run_gate): its potential reads all three from the resolved paths.
TEST(RegressionHarness, ScmCaseInputBuildsWithRepositoryRelativePaths) {
  nlohmann::json json = test_support::scm_example_json(test_support::scm_example_directory / "examplePosfile.dat");
  for (const std::string key : {"CRPlibDataSet", "SCMParFilepath", "PosFilepath"})
    json[key] = fs::relative(json[key].get<std::string>(), repository_root()).string();
  ASSERT_EQ(json["SCMParFilepath"], "potentials/H2Cu111SRP_C3v.dir/params.dat");
  const CaseSetup setup = case_setup(json);
  const Config& config = setup.config;
  EXPECT_TRUE(setup.pes->scm_potential());
  EXPECT_EQ(config.input.crp_dataset, (repository_root() / "potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v").string());
  EXPECT_EQ(config.input.scm.parameter_file, (repository_root() / "potentials/H2Cu111SRP_C3v.dir/params.dat").string());
  EXPECT_EQ(config.input.scm.position_file,
            (repository_root() / "potentials/H2Cu111SRP_C3v.dir/examplePosfile.dat").string());
  EXPECT_TRUE(config.derived.scm.has_value());
}

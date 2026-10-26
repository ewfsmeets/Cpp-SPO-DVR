// Regression gates 3 and 4 (CLAUDE.md): the propagation and the on-the-fly analysis of
// tests/reference/hardwall_short (1000 steps) against the Fortran output of the same case
// (fortran_stdout.txt, fortran_averages.out and, if present, fortran_fort21.bin, written by
// tools/fortran_reference/run_propagation_reference.sh). Skipped when missing.
// The full 7000-step case (hardwall_small) runs only with SPODVR_LONG_TESTS=1.
#include <cstdlib>

#include "analysis_reference.hpp"
#include "propagation_reference.hpp"

using namespace spodvr;
using namespace spodvr::testing;
namespace fs = std::filesystem;

namespace {

void run_gate(const std::string& case_name) {
  const fs::path dir = fs::path(SPODVR_REFERENCE_DIR) / case_name;
  const auto ref = read_fortran_reference(dir / "fortran_stdout.txt", dir / "fortran_averages.out");
  if (!ref) GTEST_SKIP() << "no Fortran reference in " << dir << "; run "
                         << "tools/fortran_reference/run_propagation_reference.sh " << case_name;
  ASSERT_FALSE(ref->averages.empty()) << "no averages found in fortran_averages.out";
  const PropagationRun cpp = run_cpp_propagation(dir / "input.json");
  compare_with_fortran(cpp, *ref);  // gate 3

  // Gate 4: the printed probabilities, and the BK coefficients of every analysis step.
  compare_probabilities_with_fortran(cpp, read_fortran_probabilities(dir / "fortran_stdout.txt"));
  if (fs::exists(dir / "fortran_fort21.bin"))
    compare_bk_with_fortran(cpp, read_fortran_bk_records(dir / "fortran_fort21.bin"), 1e-10);
  else
    std::cout << "  no fortran_fort21.bin: BK coefficients not compared\n";

  // Restart gate: a Fortran run split by a restart gives the probabilities of the unsplit run
  // (tools/fortran_reference/run_restart_reference.sh). The C++ restart is bitwise exact (test_restart).
  if (fs::exists(dir / "fortran_restart_stdout.txt"))
    compare_probabilities_with_fortran(cpp, read_fortran_probabilities(dir / "fortran_restart_stdout.txt"));
}

}  // namespace

TEST(RegressionPropagation, Gate3ShortRun) { run_gate("hardwall_short"); }

TEST(RegressionPropagation, Gate3FullRunLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  run_gate("hardwall_small");
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

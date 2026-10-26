// Regression gate 3 (CLAUDE.md): the propagation of tests/reference/hardwall_short (1000 steps)
// against the Fortran stdout of the same case (fortran_stdout.txt, written by
// tools/fortran_reference/run_propagation_reference.sh). Skipped when the file is missing.
// The full 7000-step case (hardwall_small) runs only with SPODVR_LONG_TESTS=1.
#include <cstdlib>

#include "propagation_reference.hpp"

using namespace spodvr;
using namespace spodvr::testing;
namespace fs = std::filesystem;

namespace {

void run_gate(const std::string& case_name) {
  const fs::path dir = fs::path(SPODVR_REFERENCE_DIR) / case_name;
  const auto ref = read_fortran_reference(dir / "fortran_stdout.txt");
  if (!ref) GTEST_SKIP() << "no Fortran reference " << (dir / "fortran_stdout.txt") << "; run "
                         << "tools/fortran_reference/run_propagation_reference.sh " << case_name;
  ASSERT_FALSE(ref->averages.empty()) << "no averages table found in the Fortran stdout";
  compare_with_fortran(run_cpp_propagation(dir / "input.json"), *ref);
}

}  // namespace

TEST(RegressionPropagation, Gate3ShortRun) { run_gate("hardwall_short"); }

TEST(RegressionPropagation, Gate3FullRunLong) {
  if (std::getenv("SPODVR_LONG_TESTS") == nullptr) GTEST_SKIP() << "set SPODVR_LONG_TESTS=1 to run";
  run_gate("hardwall_small");
}

TEST(RegressionPropagation, ParsesTheFortranStdoutFormat) {
  // The relevant lines in the formats of WaveFunction.f90 CoherentStateFunc and Averages.f90.
  const fs::path file = fs::temp_directory_path() / "spodvr_test_fortran_stdout.txt";
  {
    std::ofstream out(file);
    out << " * Norm initial wave packet (Z,r) =  0.10000000E+01\n"
        << " * Absolute value of the wave function in momentum space:\n"
        << "   - at k = 0                                              0.48657986E-06\n"
        << "   - at the end/beginning of the grid                      0.79631124E-16\n"
        << "   - maximum value =                                       0.54436230E+00\n"
        << "   - Low boundary (asymptotic value of Z), iLow =                                               126\n"
        << "   - High boundary (start of opt. pot on spec. grid OR a specifically set ZAnaSp ), iHigh =     235\n"
        << " * Norm of the initial wave packet between iLow and iHigh  =  0.10000000E+01\n"
        << "      time             norm           z              r        rest norm\n\n"
        << "  0.00000000E+00  0.10000000E+01  0.11440000E+02  0.14708348E+01  0.10000000E+01\n\n"
        << "  ==> ZAnalysis on specular channel fixed at Z =  0.16720000E+02\n\n"
        << "  0.10000000E+03  0.99999946E+00  0.11166768E+02  0.14708477E+01  0.99999945E+00\n\n"
        << " CPU TIMES \n  0.1 0.2 0.3 0.4 0.5\n";
  }
  const auto ref = read_fortran_reference(file);
  fs::remove(file);
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

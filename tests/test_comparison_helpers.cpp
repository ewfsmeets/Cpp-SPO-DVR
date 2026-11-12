// Tests: the NaN-aware maximum of the test comparisons (tests/comparison.hpp) and the gate
// comparisons built on it: a NaN on either side fails gate 3 (averages) and gate 4 (printed
// probabilities, BK coefficients); for finite values the maximum is std::max, bit for bit.
#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <vector>

#include "analysis_reference.hpp"
#include "comparison.hpp"
#include "propagation_reference.hpp"
#include "setup_reference.hpp"

using namespace spodvr;
using namespace spodvr::testing;
using cplx = std::complex<double>;

namespace {

constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();
constexpr double infinity = std::numeric_limits<double>::infinity();

// The largest value, reduced as the gates reduce their deviations.
double largest_deviation(const std::vector<double>& deviations) {
  double worst = 0.0;
  for (const double d : deviations) worst = nan_aware_max(worst, d);
  return worst;
}

// A C++ run whose initial wave packet block and one averages row equal `fortran_reference()`.
PropagationRun run_with_averages() {
  PropagationRun run;
  run.initial = {1.0, 0.5, 0.0, 0.8, 10, 20, 0.99};  // norm, k(0), k(edge), k_max, iLow, iHigh, inside
  run.averages = {{10.0, 1.0, 5.0, 1.4, 0.5}};       // time, norm, <Z>, <r>, rest norm
  return run;
}

FortranPropagationReference fortran_reference() {
  FortranPropagationReference ref;
  ref.norm = 1.0;
  ref.k_at_zero = 0.5;
  ref.k_at_edge = 0.0;
  ref.k_maximum = 0.8;
  ref.i_low = 11;  // 1-based
  ref.i_high = 21;
  ref.norm_inside = 0.99;
  ref.averages = {{10.0, 1.0, 5.0, 1.4, 0.5}};
  return ref;
}

// One energy with the printed total energy, sum of probabilities and reaction probability.
PropagationRun run_with_probabilities() {
  EnergyResult e{};
  e.total_energy = 0.45837825;
  e.p_total = 0.95;
  e.p_reaction = 0.05;
  PropagationRun run;
  run.scattering.energies = {e};
  return run;
}

std::vector<FortranEnergyBlock> fortran_probabilities() {
  FortranEnergyBlock block;
  block.total_energy = "0.45837825D+00";
  block.p_total = "0.9500000D+00";
  block.p_reaction = "0.5000000D-01";
  return {block};
}

// Two analysis steps of two BK coefficients.
std::vector<std::vector<cplx>> bk_records() { return {{{0.5, 0.1}, {0.2, -0.3}}, {{0.4, 0.0}, {-0.1, 0.25}}}; }

PropagationRun run_with_bk() {
  PropagationRun run;
  run.bk = bk_records();
  return run;
}

}  // namespace

// --- the maximum ------------------------------------------------------------------------------

// The defect guarded against: std::max(x, NaN) returns x, so a NaN deviation disappears.
TEST(NanAwareMax, StdMaxDropsNaN) { EXPECT_EQ(std::max(1.0, not_a_number), 1.0); }

TEST(NanAwareMax, IsStdMaxBitwiseForValuesThatAreNotNaN) {
  for (const double a : {-infinity, -1.5, -0.0, 0.0, 2e-300, 2.0, 1e308, infinity})
    for (const double b : {-infinity, -1.5, -0.0, 0.0, 2e-300, 2.0, 1e308, infinity})
      EXPECT_EQ(std::bit_cast<std::uint64_t>(nan_aware_max(a, b)), std::bit_cast<std::uint64_t>(std::max(a, b)))
          << a << ", " << b;
}

TEST(NanAwareMax, IsNaNWhenEitherArgumentIsNaN) {
  EXPECT_TRUE(std::isnan(nan_aware_max(not_a_number, 1.0)));
  EXPECT_TRUE(std::isnan(nan_aware_max(1.0, not_a_number)));
  EXPECT_TRUE(std::isnan(nan_aware_max(not_a_number, not_a_number)));
}

// Wherever the NaN is in the sequence, the reduced maximum stays NaN.
TEST(NanAwareMax, ReductionKeepsNaNAtAnyPosition) {
  EXPECT_TRUE(std::isnan(largest_deviation({not_a_number, 1e-15, 2e-15})));
  EXPECT_TRUE(std::isnan(largest_deviation({1e-15, not_a_number, 2e-15})));
  EXPECT_TRUE(std::isnan(largest_deviation({1e-15, 2e-15, not_a_number})));
  EXPECT_EQ(largest_deviation({1e-15, 3e-15, 2e-15}), 3e-15);
}

// The check of the gates on a reduced maximum (EXPECT_LE(worst, tolerance)) fails on NaN.
TEST(NanAwareMax, ToleranceCheckOfAReductionWithANaNFails) {
  const std::vector<double> deviations = {1e-15, not_a_number, 2e-15};
  const double worst = largest_deviation(deviations);
  EXPECT_NONFATAL_FAILURE(EXPECT_LE(worst, 1e-13), "worst");
}

// --- gates 1 and 2: the relative differences reduced by the setup regression tests -------------

TEST(GateComparisons, SetupRelativeDifferenceIsNaNWithANaNValue) {
  EXPECT_TRUE(std::isnan(spodvr::test::relative_difference(not_a_number, 1.0)));
  EXPECT_TRUE(std::isnan(spodvr::test::relative_difference(1.0, not_a_number)));
  EXPECT_TRUE(std::isnan(spodvr::test::relative_difference(not_a_number, not_a_number)));
}

// --- gate 3: initial wave packet block and averages ---------------------------------------------

TEST(GateComparisons, AveragesPassWhenEqual) { compare_with_fortran(run_with_averages(), fortran_reference()); }

TEST(GateComparisons, AveragesFailOnACppNaN) {
  PropagationRun run = run_with_averages();
  run.averages[0].z_mean = not_a_number;
  EXPECT_NONFATAL_FAILURE(compare_with_fortran(run, fortran_reference()), "z_mean");
}

TEST(GateComparisons, AveragesFailOnAFortranNaN) {
  FortranPropagationReference ref = fortran_reference();
  ref.averages[0].rest_norm = not_a_number;
  EXPECT_NONFATAL_FAILURE(compare_with_fortran(run_with_averages(), ref), "rest_norm");
}

// --- gate 4: printed probabilities --------------------------------------------------------------

TEST(GateComparisons, ProbabilitiesPassWhenEqual) {
  compare_probabilities_with_fortran(run_with_probabilities(), fortran_probabilities());
}

TEST(GateComparisons, ProbabilitiesFailOnACppNaN) {
  PropagationRun run = run_with_probabilities();
  run.scattering.energies[0].p_total = not_a_number;
  EXPECT_NONFATAL_FAILURE(compare_probabilities_with_fortran(run, fortran_probabilities()), "sum of the probabilities");
}

// gfortran prints a NaN as "NaN".
TEST(GateComparisons, ProbabilitiesFailOnAFortranNaN) {
  std::vector<FortranEnergyBlock> ref = fortran_probabilities();
  ref[0].p_reaction = "NaN";
  EXPECT_NONFATAL_FAILURE(compare_probabilities_with_fortran(run_with_probabilities(), ref), "reaction probability");
}

// --- gate 4: BK coefficients --------------------------------------------------------------------

TEST(GateComparisons, BkCoefficientsPassWhenEqual) { compare_bk_with_fortran(run_with_bk(), bk_records(), 1e-10); }

// A known-wrong finite coefficient fails: the comparison is live.
TEST(GateComparisons, BkCoefficientsFailOnAWrongValue) {
  PropagationRun run = run_with_bk();
  run.bk[0][1] += 1e-6;
  EXPECT_NONFATAL_FAILURE(compare_bk_with_fortran(run, bk_records(), 1e-10), "worst");
}

TEST(GateComparisons, BkCoefficientsFailOnACppNaN) {
  PropagationRun run = run_with_bk();
  run.bk[1][1] = {not_a_number, 0.25};
  EXPECT_NONFATAL_FAILURE(compare_bk_with_fortran(run, bk_records(), 1e-10), "worst");
}

// A NaN in the reference also makes the reference scale NaN instead of leaving it finite.
TEST(GateComparisons, BkCoefficientsFailOnAFortranNaN) {
  std::vector<std::vector<cplx>> ref = bk_records();
  ref[0][0] = {not_a_number, 0.1};
  EXPECT_NONFATAL_FAILURE(compare_bk_with_fortran(run_with_bk(), ref, 1e-10), "worst");
}

// --- the deviations a negative gate test asserts on ---------------------------------------------

// In printed half units, the largest over the four averages; NaN with a NaN value.
TEST(GateComparisons, AveragesDeviationInPrintedHalfUnits) {
  PropagationRun run = run_with_averages();
  EXPECT_EQ(averages_deviation(run, fortran_reference()), 0.0);
  run.averages[0].z_mean = 5.0 + 3.0 * printed_half_unit(5.0);
  // 5 + 3 half units is rounded to ulp(5) = 8.9e-16, i.e. 2e-8 half units: hence 1e-6.
  EXPECT_NEAR(averages_deviation(run, fortran_reference()), 3.0, 1e-6);
  run.averages[0].norm = not_a_number;
  EXPECT_TRUE(std::isnan(averages_deviation(run, fortran_reference())));
}

// The reference's largest |c| (0.5 + 0.1i), its scale and the largest |C++ - Fortran|; NaN with a NaN.
TEST(GateComparisons, BkDeviationGivesScaleAndLargestDifference) {
  PropagationRun run = run_with_bk();
  run.bk[1][0] += 1e-6;
  const BkDeviation d = bk_deviation(run, bk_records());
  EXPECT_EQ(d.largest, std::abs(cplx{0.5, 0.1}));
  EXPECT_EQ(d.scale, d.largest);
  EXPECT_NEAR(d.worst, 1e-6, 1e-15);
  run.bk[0][1] = {not_a_number, 0.0};
  EXPECT_TRUE(std::isnan(bk_deviation(run, bk_records()).worst));
}

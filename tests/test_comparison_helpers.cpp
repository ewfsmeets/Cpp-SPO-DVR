// Tests: the NaN-aware maximum of tests/comparison.hpp (std::max bit for bit on finite values) and
// the gates built on it: a NaN on either side fails gate 1b (the 6D potential, pes_gate.hpp), gate 3
// (averages) and gate 4 (the per-energy blocks, BK coefficients); gate 4's comparisons of the per-energy blocks
// and the analysis limits (analysis_reference.hpp) and of the S-matrix elements with the Fortran's fort.25
// (fort25_reference.hpp), and gate 5's counts (fine_analysis_reference.hpp), on synthetic data with exact counts.
#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#include "analysis_reference.hpp"
#include "comparison.hpp"
#include "core/unit_conversion.hpp"
#include "fine_analysis_reference.hpp"
#include "fort25_reference.hpp"
#include "pes_gate.hpp"
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

// The Fortran reference that run_with_averages() matches (i_low and i_high 1-based).
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

// A mini analysis basis: rotational states (j, mJ) = (0, 0), (2, 0), diffraction channels (n, m) = (0, 0) of order
// 0 and (1, 0) of order 1, levels v = 0, 1, so 4 states per level; one energy at which v = 0 is open and v = 1
// closed. Per state |S|^2 0.5, 0.1, 0.2, 0.05 in v = 0 (sum 0.85) with E_dif = 0.17 eV of order 1 and E_rot = 0.0442
// eV of j = 2: av edif = 0.17 * 0.15 / 0.85 = 0.03 eV and av erot = 0.0442 * 0.25 / 0.85 = 0.013 eV.
PropagationRun run_with_energy_block() {
  PropagationRun run;
  run.rotation_labels = {{0, 0}, {2, 0}};
  run.states = {{0, 0, 0, 0, 0}, {0, 0, 1, 1, 0}, {2, 0, 0, 0, 0}, {2, 0, 1, 1, 0}};
  run.n_states = 4;
  run.n_v = 2;
  run.n_orders = 2;
  run.order_block_printed = true;
  run.level_energy_ev = {{-4.467141, -4.46, -4.45}, {-3.9, -3.89, -3.88}};  // E(v, j), j <= 2
  run.limits = {2, 2, 2, 4, 1};                                               // jana, nviba, nana, nstata, inchna
  EnergyResult e{};
  e.total_energy = -4.267141;
  e.collision_energy = 0.2;
  e.z_energy = 0.2;
  e.total_momentum = 7.3491474;
  e.p_vibration = {0.85, 0.0};
  e.e_diffraction_v = {0.03, 0.0};
  e.e_rotation_v = {0.013, 0.0};
  e.p_rotation = {0.6, 0.25, 0.0, 0.0};                         // irotana + nana * v
  e.p_order = {0.5, 0.1, 0.2, 0.05, 0.0, 0.0, 0.0, 0.0};        // order + n_orders * (irotana + nana * v)
  e.p_total = 0.85;
  e.p_reaction = 0.15;
  e.e_vibration = 0.0;
  e.e_rotation = 0.013;
  e.e_diffraction = 0.03;
  run.scattering.energies = {e};
  return run;
}

// The block the Fortran prints for run_with_energy_block() (irot 1-based; the order block for v = 0 only).
std::vector<FortranEnergyBlock> fortran_energy_block() {
  FortranEnergyBlock block;
  block.total_energy = "-0.42671410D+01";
  block.collision_energy = "0.20000000D+00";
  block.z_energy = "0.20000000D+00";
  block.total_momentum = "0.73491474D+01";
  block.vibration = {{0, "0.8500000E+00", "0.3000000E-01", "0.1300000E-01"},
                     {1, "0.0000000E+00", "0.0000000E+00", "0.0000000E+00"}};
  block.rotation = {{0, 1, 0, 0, "0.6000000E+00"},
                    {0, 2, 2, 0, "0.2500000E+00"},
                    {1, 1, 0, 0, "0.0000000E+00"},
                    {1, 2, 2, 0, "0.0000000E+00"}};
  block.order_headers = 2;
  block.order = {{0, 0, 0, 0, "0.5000000E+00"},
                 {0, 0, 0, 1, "0.1000000E+00"},
                 {0, 2, 0, 0, "0.2000000E+00"},
                 {0, 2, 0, 1, "0.5000000E-01"}};
  block.p_total = "0.8500000";
  block.e_vibration = "0.000000";
  block.e_rotation = "0.013000";
  block.e_diffraction = "0.030000";
  block.p_reaction = "0.1500000D+00";
  return {block};
}

// Two analysis steps of two BK coefficients.
std::vector<std::vector<cplx>> bk_records() { return {{{0.5, 0.1}, {0.2, -0.3}}, {{0.4, 0.0}, {-0.1, 0.25}}}; }

// A run with the BK coefficients bk_records().
PropagationRun run_with_bk() {
  PropagationRun run;
  run.bk = bk_records();
  return run;
}

// Four points of a pes_reference.txt (hartree): the second lies above VMax and is clipped, the
// fourth (2.7 meV) below the 1 eV floor of the scaled difference.
test::PesReference pes_reference() {
  test::PesReference ref;
  ref.ntheta = 2;
  ref.nphi = 2;
  ref.points = {{0, 0, 0, 0, 0, 0, -0.1}, {0, 1, 0, 0, 0, 0, 0.5}, {1, 0, 0, 0, 0, 0, 0.2}, {1, 1, 0, 0, 0, 0, 1e-4}};
  return ref;
}

// The C++ potential at those points: the reference clipped to VMax = 0.3 hartree.
const PotentialLimits pes_limits{std::nullopt, 0.3};
std::vector<double> pes_values() { return {-0.1, 0.3, 0.2, 1e-4}; }
constexpr double pes_floor = units::legacy.ev_to_hartree(1.0);
constexpr double pes_tolerance = 1e-13;  // the 6D gate's (potential_tolerance, test_regression_pes.cpp)

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

// --- gate 1b: the 6D potential ------------------------------------------------------------------

TEST(GateComparisons, PesGatePassesWhenEqual) {
  expect_pes_gate("synthetic", pes_reference(), pes_values(), pes_limits, pes_floor, pes_tolerance);
  const PesGateResult result = compare_pes_values(pes_reference(), pes_values(), pes_limits, pes_floor, pes_tolerance);
  EXPECT_EQ(result.clipped, 1);
  EXPECT_EQ(result.max_difference, 0.0);
  EXPECT_TRUE(result.failing.empty());
}

// Below 1 eV the difference is scaled by the floor: 1e-15 hartree at V = 1e-4 hartree is 1e-11
// relative but 2.7e-14 of 1 eV, so it passes, and fails without the floor (only that point).
TEST(GateComparisons, PesGateScalesSmallValuesByTheFloor) {
  std::vector<double> values = pes_values();
  values[3] += 1e-15;
  const PesGateResult result = compare_pes_values(pes_reference(), values, pes_limits, pes_floor, pes_tolerance);
  EXPECT_TRUE(result.failing.empty());
  EXPECT_LT(result.max_difference, 0.5 * pes_tolerance);
  const PesGateResult without_floor = compare_pes_values(pes_reference(), values, pes_limits, 0.0, pes_tolerance);
  EXPECT_EQ(without_floor.failing, std::vector<std::size_t>{3});
  EXPECT_GT(without_floor.max_difference, 50.0 * pes_tolerance);
}

// With VMin the reference is clipped from below too: the clipped C++ value passes, the unclipped
// one fails at that point only.
TEST(GateComparisons, PesGateClipsTheReferenceToVMin) {
  const PotentialLimits limits{-0.05, 0.3};
  std::vector<double> values = pes_values();
  values[0] = -0.05;
  const PesGateResult result = compare_pes_values(pes_reference(), values, limits, pes_floor, pes_tolerance);
  EXPECT_EQ(result.clipped, 2);
  EXPECT_TRUE(result.failing.empty());
  EXPECT_NONFATAL_FAILURE(expect_pes_gate("synthetic", pes_reference(), pes_values(), limits, pes_floor, pes_tolerance),
                          "itheta 0 iphi 0");
  EXPECT_EQ(compare_pes_values(pes_reference(), pes_values(), limits, pes_floor, pes_tolerance).failing,
            std::vector<std::size_t>{0});
}

// A known-wrong finite value fails, and only its point (the comparison is live); the largest of
// two deviations is the worst point.
TEST(GateComparisons, PesGateFailsOnAWrongValue) {
  std::vector<double> values = pes_values();
  values[0] *= 1.0 + 5e-13;
  values[2] *= 1.0 + 5e-14;
  EXPECT_NONFATAL_FAILURE(expect_pes_gate("synthetic", pes_reference(), values, pes_limits, pes_floor, pes_tolerance),
                          "itheta 0 iphi 0");
  const PesGateResult result = compare_pes_values(pes_reference(), values, pes_limits, pes_floor, pes_tolerance);
  EXPECT_EQ(result.failing, std::vector<std::size_t>{0});
  EXPECT_EQ(result.worst, 0U);
  EXPECT_GT(result.max_difference, pes_tolerance);
}

TEST(GateComparisons, PesGateFailsOnACppNaN) {
  std::vector<double> values = pes_values();
  values[2] = not_a_number;
  EXPECT_NONFATAL_FAILURE(expect_pes_gate("synthetic", pes_reference(), values, pes_limits, pes_floor, pes_tolerance),
                          "itheta 1 iphi 0");
  const PesGateResult result = compare_pes_values(pes_reference(), values, pes_limits, pes_floor, pes_tolerance);
  EXPECT_TRUE(std::isnan(result.max_difference));
  EXPECT_EQ(result.failing, std::vector<std::size_t>{2});
  EXPECT_EQ(result.worst, 2U);
  EXPECT_EQ(result.clipped, 1);
}

// A NaN in the reference fails too, also where VMax would have clipped it, and does not count as
// clipped; with two NaN points both fail and the worst point is the first.
TEST(GateComparisons, PesGateFailsOnAFortranNaN) {
  test::PesReference ref = pes_reference();
  ref.points[1].v = not_a_number;
  EXPECT_NONFATAL_FAILURE(expect_pes_gate("synthetic", ref, pes_values(), pes_limits, pes_floor, pes_tolerance),
                          "itheta 0 iphi 1");
  ref.points[2].v = not_a_number;
  const PesGateResult result = compare_pes_values(ref, pes_values(), pes_limits, pes_floor, pes_tolerance);
  EXPECT_TRUE(std::isnan(result.max_difference));
  EXPECT_EQ(result.failing, (std::vector<std::size_t>{1, 2}));
  EXPECT_EQ(result.worst, 1U);
  EXPECT_EQ(result.clipped, 0);
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

// The momentum-space values of the initial packet: max(half a printed unit, 1e-14) absolute.
TEST(GateComparisons, MomentumValueOffByMoreThanItsHalfPrintedUnitFails) {
  PropagationRun run = run_with_averages();
  run.initial.k_space_maximum = 0.8 + 1e-8;  // half a printed unit of 0.8: 5.05e-9
  EXPECT_NONFATAL_FAILURE(compare_with_fortran(run, fortran_reference()), "k_space_maximum");
}

TEST(GateComparisons, MomentumValueAtNoiseLevelOffBy2e14Fails) {
  PropagationRun run = run_with_averages();
  run.initial.k_space_at_edge = 2e-14;  // the reference value is 0
  EXPECT_NONFATAL_FAILURE(compare_with_fortran(run, fortran_reference()), "k_space_at_edge");
}

TEST(GateComparisons, MomentumValueAtNoiseLevelOffBy5e15Passes) {
  // The printed noise of scm_cu111_srp48_small at k = 0; half a printed unit of it is 5.05e-23.
  FortranPropagationReference ref = fortran_reference();
  ref.k_at_zero = 0.27311442e-14;
  PropagationRun run = run_with_averages();
  run.initial.k_space_at_zero = ref.k_at_zero + 5e-15;
  compare_with_fortran(run, ref);
}

// The floor itself: at the noise-level references 0 (the edge) and 2.73e-15 (k = 0), whose printed
// half units are 0 and 5e-23, an error of 1.01e-14 fails and one of 0.99e-14 passes.
TEST(GateComparisons, MomentumValueAtNoiseLevelOffByJustAboveTheFloorFails) {
  FortranPropagationReference ref = fortran_reference();
  ref.k_at_zero = 0.27311442e-14;
  PropagationRun run = run_with_averages();
  run.initial.k_space_at_zero = ref.k_at_zero + 1.01e-14;
  EXPECT_NONFATAL_FAILURE(compare_with_fortran(run, ref), "k_space_at_zero");
  run = run_with_averages();
  run.initial.k_space_at_edge = 1.01e-14;  // the reference value is 0
  EXPECT_NONFATAL_FAILURE(compare_with_fortran(run, fortran_reference()), "k_space_at_edge");
}

TEST(GateComparisons, MomentumValueAtNoiseLevelOffByJustBelowTheFloorPasses) {
  FortranPropagationReference ref = fortran_reference();
  ref.k_at_zero = 0.27311442e-14;
  PropagationRun run = run_with_averages();
  run.initial.k_space_at_zero = ref.k_at_zero + 0.99e-14;
  run.initial.k_space_at_edge = 0.99e-14;  // the reference value is 0
  compare_with_fortran(run, ref);
}

TEST(GateComparisons, MomentumValueFailsOnACppNaN) {
  PropagationRun run = run_with_averages();
  run.initial.k_space_at_zero = not_a_number;
  EXPECT_NONFATAL_FAILURE(compare_with_fortran(run, fortran_reference()), "k_space_at_zero");
}

TEST(GateComparisons, MomentumValueFailsOnAFortranNaN) {
  FortranPropagationReference ref = fortran_reference();
  ref.k_at_edge = not_a_number;
  EXPECT_NONFATAL_FAILURE(compare_with_fortran(run_with_averages(), ref), "k_space_at_edge");
}

// --- gate 4: the per-energy blocks and the analysis limits -------------------------------------------------------

// Every section compared in full, the averages of the closed level (P(v) = 0) skipped.
TEST(GateComparisons, ProbabilitiesPassWhenEqual) {
  const EnergyBlockComparison c = compare_energy_blocks(run_with_energy_block(), fortran_energy_block());
  EXPECT_EQ(c.failed, 0);
  EXPECT_EQ(c.blocks, 1);
  EXPECT_EQ(c.scalars.compared, 6);
  EXPECT_EQ(c.averages.compared, 3);
  EXPECT_EQ(c.vibration.compared, 2);
  EXPECT_EQ(c.vibration_means.compared, 1);
  EXPECT_EQ(c.vibration_means.skipped, 1);
  EXPECT_EQ(c.rotation.compared, 4);
  EXPECT_EQ(c.order_headers.expected, 2);
  EXPECT_EQ(c.order.compared, 4);
  EXPECT_EQ(c.order.expected, 4);
  compare_energy_blocks_with_fortran(run_with_energy_block(), fortran_energy_block());
}

TEST(GateComparisons, ProbabilitiesFailOnACppNaN) {
  PropagationRun run = run_with_energy_block();
  run.scattering.energies[0].p_total = not_a_number;
  EXPECT_NONFATAL_FAILURE(compare_energy_blocks_with_fortran(run, fortran_energy_block()), "sum of the probabilities");
  EXPECT_EQ(compare_energy_blocks(run, fortran_energy_block()).failed, 1);
}

// gfortran prints a NaN as "NaN".
TEST(GateComparisons, ProbabilitiesFailOnAFortranNaN) {
  std::vector<FortranEnergyBlock> ref = fortran_energy_block();
  ref[0].p_reaction = "NaN";
  EXPECT_NONFATAL_FAILURE(compare_energy_blocks_with_fortran(run_with_energy_block(), ref), "reaction probability");
  EXPECT_TRUE(std::isnan(compare_energy_blocks(run_with_energy_block(), ref).worst));
}

// A NaN in a new section fails it once, also where it would decide a skip: a NaN P(v) does not skip its averages.
TEST(GateComparisons, EnergyBlockNaNFailsEverySection) {
  const auto failed_with = [](const std::function<void(PropagationRun&, FortranEnergyBlock&)>& plant) {
    PropagationRun run = run_with_energy_block();
    std::vector<FortranEnergyBlock> ref = fortran_energy_block();
    plant(run, ref[0]);
    return compare_energy_blocks(run, ref).failed;
  };
  EXPECT_EQ(
      failed_with([](PropagationRun& r, FortranEnergyBlock&) { r.scattering.energies[0].p_order[1] = not_a_number; }),
      1);
  EXPECT_EQ(failed_with([](PropagationRun&, FortranEnergyBlock& f) { f.order[1].p = "NaN"; }), 1);
  EXPECT_EQ(failed_with([](PropagationRun& r, FortranEnergyBlock&) {
              r.scattering.energies[0].e_rotation_v[0] = not_a_number;
            }),
            1);
  EXPECT_EQ(failed_with([](PropagationRun&, FortranEnergyBlock& f) { f.vibration[0].e_diffraction = "NaN"; }), 1);
  EXPECT_EQ(
      failed_with([](PropagationRun& r, FortranEnergyBlock&) { r.scattering.energies[0].e_vibration = not_a_number; }),
      1);
  EXPECT_EQ(failed_with([](PropagationRun&, FortranEnergyBlock& f) { f.e_rotation = "NaN"; }), 1);
  // P(v) of the closed level printed NaN: P(v) fails, and its averages (0 and 0) are compared, not skipped.
  PropagationRun run = run_with_energy_block();
  std::vector<FortranEnergyBlock> ref = fortran_energy_block();
  ref[0].vibration[1].p = "NaN";
  const EnergyBlockComparison c = compare_energy_blocks(run, ref);
  EXPECT_EQ(c.failed, 1);
  EXPECT_EQ(c.vibration_means.compared, 2);
  EXPECT_EQ(c.vibration_means.skipped, 0);
}

// P(v = 0, j = 0, mJ = 0, order 1) = 0.1000000E+00: half a printed unit 5.05e-8; 6e-8 off fails once, 4e-8 passes.
TEST(GateComparisons, DiffractionOrderJustAboveToleranceFailsOnce) {
  PropagationRun run = run_with_energy_block();
  run.scattering.energies[0].p_order[1] = 0.1 + 6e-8;
  const EnergyBlockComparison c = compare_energy_blocks(run, fortran_energy_block());
  EXPECT_EQ(c.failed, 1);
  EXPECT_EQ(c.order.compared, 4);
  ASSERT_EQ(c.failures.size(), 1u);
  EXPECT_NE(c.failures[0].find("(v, j, mJ, order) = (0, 0, 0, 1)"), std::string::npos) << c.failures[0];
  run.scattering.energies[0].p_order[1] = 0.1 + 4e-8;
  EXPECT_EQ(compare_energy_blocks(run, fortran_energy_block()).failed, 0);
  // An order row of a closed level (v = 1 lies above the total energy) fails, and so do its block's row count, its
  // header count (a row without a header) and its own count of rows (1 of 2 orders).
  std::vector<FortranEnergyBlock> ref = fortran_energy_block();
  ref[0].order.push_back({1, 0, 0, 0, "0.0000000E+00"});
  EXPECT_EQ(compare_energy_blocks(run_with_energy_block(), ref).failed, 4);
}

// av edif of v = 0 = 0.3000000E-01: half a printed unit 5.05e-9; 6e-9 off fails once, 4e-9 passes. The averages of
// v = 1 (P(v) = 0, noise) are skipped: any value there passes.
TEST(GateComparisons, MeanEnergyJustAboveToleranceFailsOnce) {
  PropagationRun run = run_with_energy_block();
  run.scattering.energies[0].e_diffraction_v[0] = 0.03 + 6e-9;
  const EnergyBlockComparison c = compare_energy_blocks(run, fortran_energy_block());
  EXPECT_EQ(c.failed, 1);
  ASSERT_EQ(c.failures.size(), 1u);
  EXPECT_NE(c.failures[0].find("v = 0 av edif"), std::string::npos) << c.failures[0];
  run.scattering.energies[0].e_diffraction_v[0] = 0.03 + 4e-9;
  EXPECT_EQ(compare_energy_blocks(run, fortran_energy_block()).failed, 0);
  run = run_with_energy_block();
  run.scattering.energies[0].e_rotation_v[0] = 0.013 + 6e-9;  // half a unit of 0.1300000E-01: 5.05e-9
  EXPECT_EQ(compare_energy_blocks(run, fortran_energy_block()).failed, 1);
  run = run_with_energy_block();
  run.scattering.energies[0].e_rotation_v[1] = 1.0;
  const EnergyBlockComparison skipped = compare_energy_blocks(run, fortran_energy_block());
  EXPECT_EQ(skipped.failed, 0);
  EXPECT_EQ(skipped.vibration_means.skipped, 1);
  // Above the noise threshold the same row is compared: P(v = 1) printed as 2 thresholds fails its average (and
  // its P(v)).
  std::vector<FortranEnergyBlock> ref = fortran_energy_block();
  ref[0].vibration[1].p = std::format("{:.7E}", 2.0 * average_energy_noise_threshold);
  const EnergyBlockComparison compared = compare_energy_blocks(run, ref);
  EXPECT_EQ(compared.vibration_means.skipped, 0);
  EXPECT_EQ(compared.failed, 2);
}

// The three averages (f12.6: half a unit 5.05e-7): 6e-7 off fails once each, 4e-7 passes; where the printed P(v)
// add up to noise they are skipped (and so are the averages of each v).
TEST(GateComparisons, AverageEnergiesJustAboveToleranceFailOnce) {
  for (double EnergyResult::* average :
       {&EnergyResult::e_vibration, &EnergyResult::e_rotation, &EnergyResult::e_diffraction}) {
    PropagationRun run = run_with_energy_block();
    EnergyResult& e = run.scattering.energies[0];
    e.*average = e.*average + 6e-7;
    EXPECT_EQ(compare_energy_blocks(run, fortran_energy_block()).failed, 1);
    e.*average = e.*average - 2e-7;
    EXPECT_EQ(compare_energy_blocks(run, fortran_energy_block()).failed, 0);
  }
  // The probabilities of hardwall_short (~1e-29), whose sum prints as 0.0000000: the averages (0.5 eV against the
  // printed 0.013) are skipped, the three totals and those of v = 0.
  PropagationRun run = run_with_energy_block();
  EnergyResult& e = run.scattering.energies[0];
  e.p_vibration[0] = 1e-29;
  e.p_total = 1e-29;
  e.p_reaction = 1.0;
  e.e_rotation = 0.5;
  e.e_rotation_v[0] = 0.5;
  std::vector<FortranEnergyBlock> ref = fortran_energy_block();
  ref[0].vibration[0].p = "0.1000000E-28";
  ref[0].p_total = "0.0000000";
  ref[0].p_reaction = "0.1000000D+01";
  EnergyBlockComparison c = compare_energy_blocks(run, ref);
  EXPECT_EQ(c.failed, 0);
  EXPECT_EQ(c.averages.skipped, 3);
  EXPECT_EQ(c.averages.compared, 0);
  EXPECT_EQ(c.vibration_means.skipped, 2);
  // P(v = 0) twice the threshold: both averages of 0.5 eV are compared and fail, the total and that of v = 0.
  ref[0].vibration[0].p = std::format("{:.7E}", 2.0 * average_energy_noise_threshold);
  e.p_vibration[0] = 2.0 * average_energy_noise_threshold;
  c = compare_energy_blocks(run, ref);
  EXPECT_EQ(c.averages.compared, 3);
  EXPECT_EQ(c.vibration_means.skipped, 1);
  EXPECT_EQ(c.failed, 2);
}

// A rotation row whose (irot, j, mJ) is not a C++ rotational state fails, once; an order header with an unknown
// (j, mJ) fails each of its rows.
TEST(GateComparisons, UnmatchedRotationLabelFails) {
  std::vector<FortranEnergyBlock> ref = fortran_energy_block();
  ref[0].rotation[1].j = 4;  // irot 2 is (2, 0)
  EXPECT_EQ(compare_energy_blocks(run_with_energy_block(), ref).failed, 1);
  ref = fortran_energy_block();
  ref[0].rotation[1].irot = 1;  // (2, 0) is irot 2, not 1
  EXPECT_EQ(compare_energy_blocks(run_with_energy_block(), ref).failed, 1);
  ref = fortran_energy_block();
  ref[0].order[2].mj = 1;
  ref[0].order[3].mj = 1;  // the header (0, 2, 1): no such state
  const EnergyBlockComparison c = compare_energy_blocks(run_with_energy_block(), ref);
  EXPECT_EQ(c.failed, 2);
  EXPECT_EQ(c.order.compared, 2);
  ref = fortran_energy_block();
  ref[0].vibration[1].v = 2;  // no level v = 2
  EXPECT_EQ(compare_energy_blocks(run_with_energy_block(), ref).failed, 1);
}

// A section left empty (or a row missing) fails the block's count once, and the gate's count of the section.
TEST(GateComparisons, EmptySectionFailsTheCount) {
  const std::vector<std::string> sections = {"P(v)", "P(v, j, mJ)", "P(v, j, mJ, order)"};
  for (const std::string& section : sections) {
    std::vector<FortranEnergyBlock> ref = fortran_energy_block();
    if (section == "P(v)") ref[0].vibration.clear();
    if (section == "P(v, j, mJ)") ref[0].rotation.clear();
    if (section == "P(v, j, mJ, order)") ref[0].order.clear();
    const EnergyBlockComparison c = compare_energy_blocks(run_with_energy_block(), ref);
    // An empty order section also leaves its 2 headers without rows (the per-header count).
    const std::size_t failures = section == "P(v, j, mJ, order)" ? 2 : 1;
    EXPECT_EQ(c.failed, static_cast<int>(failures)) << section;
    ASSERT_EQ(c.failures.size(), failures) << section;
    EXPECT_NE(c.failures[0].find("0 " + section + " rows read, the print rule gives"), std::string::npos)
        << c.failures[0];
  }
  // One header fewer than the print rule: the header count, and 1 header for 2 (v, j, mJ) with rows.
  std::vector<FortranEnergyBlock> ref = fortran_energy_block();
  ref[0].order_headers = 1;
  EXPECT_EQ(compare_energy_blocks(run_with_energy_block(), ref).failed, 2);
  ref = fortran_energy_block();
  ref[0].e_diffraction.clear();  // an average not printed: the count and the value
  EXPECT_EQ(compare_energy_blocks(run_with_energy_block(), ref).failed, 2);
  ref = fortran_energy_block();
  ref[0].rotation.push_back(ref[0].rotation[0]);  // a repeated row
  EXPECT_EQ(compare_energy_blocks(run_with_energy_block(), ref).failed, 1);
  // Without writeRDProbs the Fortran prints no order block: the rows of one fail the count.
  PropagationRun run = run_with_energy_block();
  run.order_block_printed = false;
  EXPECT_EQ(compare_energy_blocks(run, fortran_energy_block()).failed, 2);  // the headers and the rows
  ref = fortran_energy_block();
  ref[0].order.clear();
  ref[0].order_headers = 0;
  EXPECT_EQ(compare_energy_blocks(run, ref).failed, 0);
  // An energy block missing.
  EXPECT_EQ(compare_energy_blocks(run_with_energy_block(), {}).failed, 1);
}

// The gate's assertions (expect_energy_blocks_gate), each failing alone.
TEST(GateComparisons, EnergyBlocksGateAssertsItsCounts) {
  const EnergyBlockComparison c = compare_energy_blocks(run_with_energy_block(), fortran_energy_block());
  expect_energy_blocks_gate(c);
  EnergyBlockComparison failed = c;
  failed.failed = 1;
  EXPECT_NONFATAL_FAILURE(expect_energy_blocks_gate(failed), "failed values, rows or counts");
  EnergyBlockComparison short_read = c;
  short_read.rotation.read = 3;
  short_read.rotation.compared = 3;
  EXPECT_NONFATAL_FAILURE(expect_energy_blocks_gate(short_read), "P(v, j, mJ) rows read vs the print rule");
  EnergyBlockComparison not_compared = c;
  not_compared.order.compared = 3;
  EXPECT_NONFATAL_FAILURE(expect_energy_blocks_gate(not_compared), "P(v, j, mJ, order) rows compared or skipped");
  EnergyBlockComparison not_skipped = c;
  not_skipped.vibration_means.skipped = 0;
  EXPECT_NONFATAL_FAILURE(expect_energy_blocks_gate(not_skipped), "av edif / av erot rows compared or skipped");
  EnergyBlockComparison headers = c;
  headers.order_headers.read = 1;
  EXPECT_NONFATAL_FAILURE(expect_energy_blocks_gate(headers), "order headers read vs the print rule");
  EnergyBlockComparison empty = c;
  empty.blocks = 0;
  EXPECT_NONFATAL_FAILURE(expect_energy_blocks_gate(empty), "no energy block read");
}

// The limits each against the C++ basis: one wrong value fails once; a limit not printed, printed twice or printed
// where the stdout must not print it fails once too.
TEST(GateComparisons, AnalysisLimitsMismatchFails) {
  const AnalysisLimits cpp = run_with_energy_block().limits;
  const FortranAnalysisLimits printed = {
      {"jana", {"2"}}, {"nviba", {"2"}}, {"nana", {"2"}}, {"nstata", {"4"}}, {"inchna", {"1"}}};
  const LimitsComparison equal = compare_analysis_limits(cpp, printed, limits_of_a_run);
  EXPECT_EQ(equal.failed, 0);
  EXPECT_EQ(equal.compared, 5);
  expect_analysis_limits(equal, limits_of_a_run);
  for (const auto& [name, values] : printed) {
    FortranAnalysisLimits wrong = printed;
    wrong[name] = {"3"};
    EXPECT_EQ(compare_analysis_limits(cpp, wrong, limits_of_a_run).failed, 1) << name;
    FortranAnalysisLimits missing = printed;
    missing.erase(name);
    EXPECT_EQ(compare_analysis_limits(cpp, missing, limits_of_a_run).failed, 1) << name;
    FortranAnalysisLimits twice = printed;
    twice[name].push_back(values[0]);
    EXPECT_EQ(compare_analysis_limits(cpp, twice, limits_of_a_run).failed, 1) << name;
  }
  // A set-up stdout prints nana, nstata and inchna only.
  EXPECT_EQ(compare_analysis_limits(cpp, printed, limits_of_a_setup).failed, 2);
  FortranAnalysisLimits unreadable = printed;
  unreadable["nana"] = {"***"};
  const LimitsComparison c = compare_analysis_limits(cpp, unreadable, limits_of_a_run);
  EXPECT_EQ(c.failed, 1);
  EXPECT_NONFATAL_FAILURE(expect_analysis_limits(c, limits_of_a_run), "nana: C++ 2, Fortran ***");
  LimitsComparison not_compared = equal;
  not_compared.compared = 4;
  EXPECT_NONFATAL_FAILURE(expect_analysis_limits(not_compared, limits_of_a_run), "analysis limits compared");
}

// A printed E-format zero is exact (a nonzero value prints a nonzero mantissa), an F-format zero is not (f10.6
// prints 3e-7 as 0.000000); in the energy blocks a C++ 1e-15 under a printed zero fails once (half a unit, 5.05e-8,
// let it pass before).
TEST(GateComparisons, PrintedZeroIsExactInEFormatOnly) {
  EXPECT_EQ(printed_half_unit("0.0000000E+00", 1e-20), 1e-20);
  EXPECT_EQ(printed_half_unit("-0.000000D+00", 0.0), 0.0);
  EXPECT_EQ(printed_half_unit("0.000000", 0.0), 1.01 * 0.5e-6);
  EXPECT_EQ(printed_half_unit("0.1000000E+00", 0.0), printed_text_tolerance("0.1000000E+00", 0.0));
  PropagationRun run = run_with_energy_block();
  run.scattering.energies[0].p_rotation[2] = 1e-15;  // printed 0.0000000E+00: (v, irot, j, mJ) = (1, 1, 0, 0)
  const EnergyBlockComparison c = compare_energy_blocks(run, fortran_energy_block());
  EXPECT_EQ(c.failed, 1);
  ASSERT_EQ(c.failures.size(), 1u);
  EXPECT_NE(c.failures[0].find("(v, irot, j, mJ) = (1, 1, 0, 0)"), std::string::npos) << c.failures[0];
}

// The order rows are counted per header: one header with a row too many and one with a row too few (the block's
// total unchanged) fail once each, and so does a header printed twice.
TEST(GateComparisons, OrderRowsAreCountedPerHeader) {
  const auto order_count_failures = [](const EnergyBlockComparison& c) {
    return std::ranges::count_if(c.failures, [](const std::string& s) {
      return s.find("order rows under") != std::string::npos || s.find("order headers for") != std::string::npos;
    });
  };
  EXPECT_EQ(order_count_failures(compare_energy_blocks(run_with_energy_block(), fortran_energy_block())), 0);
  std::vector<FortranEnergyBlock> moved = fortran_energy_block();
  moved[0].order[3] = {0, 0, 0, 1, "0.1000000E+00"};  // the second header's order-1 row under the first header
  EXPECT_EQ(order_count_failures(compare_energy_blocks(run_with_energy_block(), moved)), 2);
  std::vector<FortranEnergyBlock> twice = fortran_energy_block();
  twice[0].order_headers = 3;  // a header printed twice (its rows read under one key)
  EXPECT_EQ(order_count_failures(compare_energy_blocks(run_with_energy_block(), twice)), 1);
}

// The averages skipped as noise are not compared, but a non-finite C++ value there fails once: av erot of the closed
// level (P(v) = 0), and the three totals when every printed P(v) is 0.
TEST(GateComparisons, SkippedAverageMustBeFinite) {
  PropagationRun run = run_with_energy_block();
  run.scattering.energies[0].e_rotation_v[1] = std::numeric_limits<double>::quiet_NaN();
  const EnergyBlockComparison per_v = compare_energy_blocks(run, fortran_energy_block());
  EXPECT_EQ(per_v.failed, 1);
  EXPECT_EQ(per_v.vibration_means.skipped, 1);

  std::vector<FortranEnergyBlock> closed = fortran_energy_block();
  closed[0].vibration[0].p = "0.0000000E+00";  // every printed P(v) 0: the totals are skipped
  const int base = compare_energy_blocks(run_with_energy_block(), closed).failed;
  PropagationRun total = run_with_energy_block();
  total.scattering.energies[0].e_rotation = std::numeric_limits<double>::infinity();
  const EnergyBlockComparison totals = compare_energy_blocks(total, closed);
  EXPECT_EQ(totals.averages.skipped, 3);
  EXPECT_EQ(totals.failed, base + 1);
}

// The fixed specular line (gate 3): equal to half a printed unit (e15.8 of 10.88: 5e-7, with the 1 % margin 5.05e-7)
// passes; just beyond it, a NaN, a line fixed in only one code, a line printed twice or an unreadable value fails
// once; a line fixed in neither passes.
TEST(GateComparisons, FixedLineMismatchFails) {
  const std::vector<std::string> printed = {"0.10880000E+02"};
  expect_fixed_line(10.88, printed);
  expect_fixed_line(10.88 + 4.9e-7, printed);
  expect_fixed_line(std::nullopt, {});
  EXPECT_NONFATAL_FAILURE(expect_fixed_line(10.88 + 5.2e-7, printed), "fixed specular line: C++ Z");
  EXPECT_NONFATAL_FAILURE(expect_fixed_line(std::numeric_limits<double>::quiet_NaN(), printed),
                          "fixed specular line: C++ Z");
  EXPECT_NONFATAL_FAILURE(expect_fixed_line(std::nullopt, printed), "specular line fixed: C++ no, Fortran yes");
  EXPECT_NONFATAL_FAILURE(expect_fixed_line(10.88, {}), "specular line fixed: C++ yes, Fortran no");
  EXPECT_NONFATAL_FAILURE(expect_fixed_line(10.88, {printed[0], printed[0]}),
                          "printed its fixed specular line 2 times");
  EXPECT_NONFATAL_FAILURE(expect_fixed_line(10.88, {"***************"}), "is not a finite number");
}

// The reader of the limits: the two print formats, a name inside another word not taken.
TEST(GateComparisons, AnalysisLimitsReader) {
  const std::filesystem::path path = std::filesystem::path(::testing::TempDir()) / "limits_synthetic.txt";
  std::ofstream(path) << " * number of rotational states analysed, nana =              15\n"
                         " * Number of rotation-diffraction analysed, nstata =              240\n"
                         " * Index of initial rot-dif channel in basis, inchna =             97\n"
                         " * Maximum j in analysis:                        jana =   4\n"
                         " * Maximum v in analysis:                        nviba =   2\n"
                         " xnana = 3\n";
  const FortranAnalysisLimits printed = read_fortran_analysis_limits(path);
  const FortranAnalysisLimits expected = {
      {"jana", {"4"}}, {"nviba", {"2"}}, {"nana", {"15"}}, {"nstata", {"240"}}, {"inchna", {"97"}}};
  EXPECT_EQ(printed, expected);
}

// The reader of the energy blocks on the Fortran's formats (Analyzor.f90 smatrix_combine) gives the fixture.
TEST(GateComparisons, EnergyBlockReader) {
  const std::filesystem::path path = std::filesystem::path(::testing::TempDir()) / "energy_block_synthetic.txt";
  std::ofstream(path)
      << "   0  0.1  0.2  0.3\n"  // before the first block: not read
      << "\n total energy = -0.42671410D+01\n trans energy =  0.20000000D+00\n"
      << " trans energy in z =  0.20000000D+00\n total momentum =  0.73491474D+01\n\n"
      << " ************************************\n rotation-diffraction summed probabilities\n"
      << " ************************************\n\n v    probability         av edif         av erot\n\n"
      << "   0  0.8500000E+00  0.3000000E-01  0.1300000E-01\n\n"
      << "   1  0.0000000E+00  0.0000000E+00  0.0000000E+00\n\n"
      << " ************************************\n diffraction summed probabilities\n"
      << " ************************************\n\n  v   irot   j   m  probability\n\n"
      << "  0    1     0    0  0.6000000E+00\n\n  0    2     2    0  0.2500000E+00\n\n"
      << "  1    1     0    0  0.0000000E+00\n\n  1    2     2    0  0.0000000E+00\n\n"
      << " ************************************\n probabilities as function of diffraction order\n"
      << " ************************************\n\n v, j, mj =   0   0   0\n\n"
      << "   0  0.5000000E+00\n\n   1  0.1000000E+00\n\n v, j, mj =   0   2   0\n\n"
      << "   0  0.2000000E+00\n\n   1  0.5000000E-01\n\n"
      << " sum of the probabilities = 0.8500000\n\n"
      << " average vibrational energy (ev) =     0.000000\n\n"
      << " average rotational  energy (ev) =     0.013000\n\n"
      << " average diffraction energy (ev) =     0.030000\n\n"
      << " the reaction probability =  0.1500000D+00\n\n"
      << "   1  0.2\n";  // after the tables: not an order row
  const std::vector<FortranEnergyBlock> blocks = read_fortran_energy_blocks(path);
  ASSERT_EQ(blocks.size(), 1u);
  EXPECT_TRUE(blocks[0].problems.empty());
  const EnergyBlockComparison c = compare_energy_blocks(run_with_energy_block(), blocks);
  EXPECT_EQ(c.failed, 0);
  EXPECT_EQ(c.order.compared, 4);
  EXPECT_EQ(c.rotation.compared, 4);
  EXPECT_EQ(c.vibration.compared, 2);
  EXPECT_EQ(c.averages.compared, 3);
  // Lines the reader cannot use are problems, each one failure: an order row before its header, an unreadable
  // header (its row has no header either) and an unreadable rotation label.
  std::ofstream(path) << " total energy = -0.42671410D+01\n probabilities as function of diffraction order\n"
                      << "   0  0.5000000E+00\n v, j, mj =   0   *   0\n   1  0.1000000E+00\n"
                      << " diffraction summed probabilities\n  0    1     x    0  0.6000000E+00\n";
  const std::vector<FortranEnergyBlock> faulty = read_fortran_energy_blocks(path);
  ASSERT_EQ(faulty.size(), 1u);
  EXPECT_EQ(faulty[0].problems.size(), 4u);
  EXPECT_EQ(faulty[0].order_headers, 1);
  EXPECT_TRUE(faulty[0].order.empty());
}

// --- gate 5: the counts -----------------------------------------------------------------------------------------

// The counts of the Fortran's write rules on the mini basis with nen95 = 3: fil10 3 * (7 + 3 * 2 + 2 * 2 * (1 + 2))
// = 75 values, unit 3 3 * 4 * 2 = 24, fort.25 the 4 open states of v = 0, three values of each not skipped (1 here);
// each count one short fails once.
TEST(GateComparisons, FineCountsFailWhenShort) {
  const FineCounts expected = expected_fine_counts(run_with_energy_block(), 3, 1);
  EXPECT_EQ(expected.f10, 75);
  EXPECT_EQ(expected.u3, 24);
  EXPECT_EQ(expected.p25, 4);
  EXPECT_EQ(expected.s25, 9);
  expect_fine_counts(expected, expected);
  FineCounts f10 = expected;
  f10.f10 = f10.f10 - 1;
  EXPECT_NONFATAL_FAILURE(expect_fine_counts(f10, expected), "fil10 values compared");
  FineCounts u3 = expected;
  u3.u3 = u3.u3 - 1;
  EXPECT_NONFATAL_FAILURE(expect_fine_counts(u3, expected), "unit 3 values compared");
  FineCounts p25 = expected;
  p25.p25 = p25.p25 - 1;
  EXPECT_NONFATAL_FAILURE(expect_fine_counts(p25, expected), "fort.25 probabilities compared");
  FineCounts s25 = expected;
  s25.s25 = s25.s25 - 1;
  EXPECT_NONFATAL_FAILURE(expect_fine_counts(s25, expected), "fort.25 S and phases compared");
  // Without writeRDProbs fil10 has no P per diffraction order: 3 * (7 + 6 + 4) = 51.
  PropagationRun run = run_with_energy_block();
  run.order_block_printed = false;
  EXPECT_EQ(expected_fine_counts(run, 3, 1).f10, 51);
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

// ---------------------------------------------------------------------------------------------------------------
// Gate 4, the S-matrix elements (fort25_reference.hpp): a run of two states at one energy and the fort.25 the
// Fortran prints for it, each fault planted on top.

namespace {

// x as the Fortran's D13.6 prints it (0.dddddd and the exponent).
std::string fortran_d(double x) {
  if (x == 0.0) return "0.000000D+00";
  int exponent = static_cast<int>(std::floor(std::log10(std::abs(x)))) + 1;
  double mantissa = std::round(x / std::pow(10.0, exponent) * 1e6) / 1e6;
  if (std::abs(mantissa) >= 1.0) {
    mantissa = mantissa / 10.0;
    exponent = exponent + 1;
  }
  return std::format("{:.6f}D{}{:02d}", mantissa, exponent < 0 ? '-' : '+', std::abs(exponent));
}

// The run: states (j, mJ, order, n, m) = (0, 0, 0, 0, 0) and (0, 0, 1, 1, 0), one level E(0, 0) below the total
// energy, so the print rule expects both rows.
struct Fort25Fixture {
  PropagationRun run;
  std::vector<cplx> s = {cplx(-0.470660e-03, -0.850680e-03), cplx(0.100000e-01, 0.200000e-01)};

  Fort25Fixture() {
    run.n_states = 2;
    run.n_v = 1;
    run.states = {{0, 0, 0, 0, 0}, {0, 0, 1, 1, 0}};
    run.level_energy_ev = {{-4.467141}};
    EnergyResult e{};
    e.total_energy = -4.267141;
    run.scattering.energies = {e};
    set_cpp(s);
  }
  // The C++ S, P and phase of the two states.
  void set_cpp(const std::vector<cplx>& values) {
    EnergyResult& e = run.scattering.energies[0];
    e.s_matrix = values;
    e.probability = {std::norm(values[0]), std::norm(values[1])};
    e.phase = {std::arg(values[0]), std::arg(values[1])};
  }
  // The row the Fortran prints for state k (Re S, Im S, P and the phase from s[k]).
  std::string row(int k, int n, int m) const {
    return std::format("{:4d}{:4d}  {:>13} {:>13} {:>13} {:10.6f}", n, m, fortran_d(s[k].real()),
                       fortran_d(s[k].imag()), fortran_d(std::norm(s[k])), std::arg(s[k]));
  }
  std::string text() const {
    return " total energy = -0.42671410D+01\n"
           " v j mj dif.order are:   0   0   0   0\n"
           "  m1  m2        s-matrix element          prob          phase\n" +
           row(0, 0, 0) +
           "\n"
           " v j mj dif.order are:   0   0   0   1\n"
           "  m1  m2        s-matrix element          prob          phase\n" +
           row(1, 1, 0) + "\n";
  }
};

// The comparison of run with the file text, tolerance's code term s_scale_fraction.
Fort25Comparison compare_text(const PropagationRun& run, const std::string& text, double s_scale_fraction = 1e-10) {
  const std::filesystem::path path = std::filesystem::path(::testing::TempDir()) / "fort25_synthetic.txt";
  std::ofstream(path) << text;
  Fort25Tolerance tolerance;
  tolerance.s_scale_fraction = s_scale_fraction;
  return compare_fort25(run, read_fortran_fort25_rows(path), tolerance);
}

std::string replaced(std::string text, const std::string& old, const std::string& replacement) {
  const auto at = text.find(old);
  EXPECT_NE(at, std::string::npos) << old;
  return at == std::string::npos ? text : text.replace(at, old.size(), replacement);
}

}  // namespace

TEST(GateComparisons, Fort25PassesWhenEqual) {
  const Fort25Fixture f;
  const Fort25Comparison c = compare_text(f.run, f.text());
  EXPECT_EQ(c.failed, 0);
  EXPECT_EQ(c.rows, 2);
  EXPECT_EQ(c.expected_rows, 2);
  EXPECT_EQ(c.row_like_lines, 2);
  EXPECT_EQ(c.compared, 2);
  EXPECT_EQ(c.phase_compared, 2);
}

// Re S of the C++ off by 1e-9: above the printed half units (7.1e-10) plus the code term (1e-10 * 0.022).
TEST(GateComparisons, Fort25ValueJustAboveToleranceFailsOnce) {
  Fort25Fixture f;
  const std::string text = f.text();
  f.set_cpp({f.s[0] + cplx(1e-9, 0.0), f.s[1]});
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
}

TEST(GateComparisons, Fort25FailsOnACppNaN) {
  Fort25Fixture f;
  const std::string text = f.text();
  f.set_cpp({cplx(not_a_number, 0.0), f.s[1]});
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
}

TEST(GateComparisons, Fort25FailsOnAFortranNaNToken) {
  const Fort25Fixture f;
  const Fort25Comparison c = compare_text(f.run, replaced(f.text(), "0.100000D-01", "NaN"));
  EXPECT_EQ(c.failed, 1);
  EXPECT_EQ(c.compared, 1);
}

TEST(GateComparisons, Fort25FailsOnAnInfinityOrAsteriskToken) {
  const Fort25Fixture f;
  EXPECT_EQ(compare_text(f.run, replaced(f.text(), "0.100000D-01", "Infinity")).failed, 1);
  EXPECT_EQ(compare_text(f.run, replaced(f.text(), "0.100000D-01", "*************")).failed, 1);
}

TEST(GateComparisons, Fort25ReadsAThreeDigitExponent) {
  EXPECT_EQ(fortran_number("0.100000-119"), 0.1e-119);
  EXPECT_EQ(fortran_number("-0.250000+123"), -0.25e123);
  EXPECT_EQ(fortran_number("0.500000D+00"), 0.5);
  EXPECT_FALSE(fortran_number("NaN"));
  EXPECT_FALSE(fortran_number("Infinity"));
  EXPECT_FALSE(fortran_number("*********"));
  EXPECT_EQ(parse_fortran_real("0.100000-119"), 0.1e-119);
  EXPECT_EQ(parse_fortran_real("1.5e-03"), 1.5e-3);
  // The printed digits of such a value: half a unit of its 6th digit, not of an exponent 0.
  EXPECT_EQ(fortran_text_with_e("0.100000-119"), "0.100000E-119");
  EXPECT_EQ(printed_text_tolerance(fortran_text_with_e("0.100000-119"), 0.0), 1.01 * 0.5 * std::pow(10.0, -125));
}

// printed_text_tolerance reads every exponent form gfortran prints: a three-digit exponent without its letter (it
// read that as exponent 0 before: half a unit of 1e-7 instead of 1e-107), E, D, none; a text without a decimal point
// (NaN, a non-number) gets a NaN tolerance, so every comparison with it fails.
TEST(GateComparisons, PrintedToleranceReadsEveryExponentForm) {
  EXPECT_DOUBLE_EQ(printed_text_tolerance("0.1234567-100", 0.0), 5.05e-108);
  EXPECT_DOUBLE_EQ(printed_text_tolerance("0.1234567+100", 0.0), 5.05e92);
  EXPECT_DOUBLE_EQ(printed_text_tolerance("-0.1234567E-05", 0.0), 5.05e-13);
  EXPECT_DOUBLE_EQ(printed_text_tolerance("0.1234567D+03", 0.0), 5.05e-5);
  EXPECT_DOUBLE_EQ(printed_text_tolerance("0.1234567", 0.0), 5.05e-8);
  EXPECT_TRUE(std::isnan(printed_text_tolerance("NaN", 0.0)));
  EXPECT_TRUE(std::isnan(printed_text_tolerance("*********", 1e-20)));
}

// The second block's header with diffraction order 2: its row matches no C++ state of that order.
TEST(GateComparisons, Fort25WrongOrderHeaderFails) {
  const Fort25Fixture f;
  EXPECT_EQ(compare_text(f.run, replaced(f.text(), "   0   0   0   1", "   0   0   0   2")).failed, 1);
}

TEST(GateComparisons, Fort25DuplicateKeyFails) {
  const Fort25Fixture f;
  const std::string text = f.text() + f.row(1, 1, 0) + "\n";
  const Fort25Comparison c = compare_text(f.run, text);
  EXPECT_EQ(c.failed, 1);
  EXPECT_EQ(c.rows, 3);
  EXPECT_NE(c.rows, c.expected_rows);
}

TEST(GateComparisons, Fort25UnmatchedRowFails) {
  const Fort25Fixture f;
  EXPECT_EQ(compare_text(f.run, replaced(f.text(), f.row(1, 1, 0), f.row(1, 2, 0))).failed, 1);
}

TEST(GateComparisons, Fort25MissingRowsFailTheCount) {
  const Fort25Fixture f;
  const Fort25Comparison c = compare_text(f.run, replaced(f.text(), f.row(1, 1, 0) + "\n", ""));
  EXPECT_EQ(c.failed, 0);
  EXPECT_EQ(c.rows, 1);
  EXPECT_EQ(c.row_like_lines, 1);
  EXPECT_EQ(c.expected_rows, 2);
}

// Without the second header its row belongs to the first block (order 0), where the C++ has no such state.
TEST(GateComparisons, Fort25CollapsedHeaderFails) {
  const Fort25Fixture f;
  const std::string header =
      " v j mj dif.order are:   0   0   0   1\n  m1  m2        s-matrix element          prob          phase\n";
  EXPECT_EQ(compare_text(f.run, replaced(f.text(), header, "")).failed, 1);
}

// arg S of the C++ just below pi, the Fortran's printed phase just above -pi: the same angle.
TEST(GateComparisons, Fort25PhaseWrapsAround) {
  Fort25Fixture f;
  f.s = {cplx(-0.01, 1e-12), f.s[1]};
  f.set_cpp(f.s);
  std::string text = f.text();
  text = replaced(text, "0.100000D-11", "-0.100000D-11");
  text = replaced(text, std::format("{:10.6f}", std::arg(f.s[0])), std::format("{:10.6f}", -std::numbers::pi));
  EXPECT_EQ(compare_text(f.run, text, 1e-9).failed, 0);
}

// The Fortran prints open channels only: a row of a state whose level lies above the total energy fails.
TEST(GateComparisons, Fort25ClosedChannelRowFails) {
  Fort25Fixture f;
  const std::string text = f.text();
  f.run.level_energy_ev = {{-4.0}};  // E(0, 0) above the total energy -4.267141 eV
  const Fort25Comparison c = compare_text(f.run, text);
  EXPECT_EQ(c.expected_rows, 0);
  EXPECT_EQ(c.failed, 2);
}

// A block's printed total energy one unit of its 8th digit away from the C++ energy fails once.
TEST(GateComparisons, Fort25WrongTotalEnergyFails) {
  const Fort25Fixture f;
  EXPECT_EQ(compare_text(f.run, replaced(f.text(), "-0.42671410D+01", "-0.42671411D+01")).failed, 1);
}

// Each of the three checks of a row on its own: S, P or the phase of the C++ off, the other two as printed.
TEST(GateComparisons, Fort25SOnlyDeviationFails) {
  Fort25Fixture f;
  const std::string text = f.text();
  f.run.scattering.energies[0].s_matrix[0] += cplx(1e-9, 0.0);  // above 7.1e-10 of digits + 2.2e-12 of code
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
}

TEST(GateComparisons, Fort25ProbabilityOnlyDeviationFails) {
  Fort25Fixture f;
  const std::string text = f.text();
  f.run.scattering.energies[0].probability[0] += 1.2e-12;  // P 0.945177D-06: half a unit 5.05e-13
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
}

TEST(GateComparisons, Fort25PhaseOnlyDeviationFails) {
  Fort25Fixture f;
  const std::string text = f.text();
  f.run.scattering.energies[0].phase[0] += 1e-5;  // the tolerance: 5.05e-7 + code / |S| = 5.07e-7
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
}

// A NaN in the C++ P or phase alone fails, also on a row whose phase is not compared (Re S printed 0).
TEST(GateComparisons, Fort25NaNProbabilityOrPhaseFails) {
  Fort25Fixture f;
  std::string text = f.text();
  f.run.scattering.energies[0].probability[0] = not_a_number;
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
  f.set_cpp(f.s);
  f.run.scattering.energies[0].phase[0] = not_a_number;
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
  f.s[0] = cplx(0.0, -0.850680e-03);
  f.set_cpp(f.s);
  text = f.text();
  const Fort25Comparison clean = compare_text(f.run, text);
  EXPECT_EQ(clean.failed, 0);
  EXPECT_EQ(clean.phase_compared, 1);
  f.run.scattering.energies[0].phase[0] = not_a_number;
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
}

// The gate's assertions (expect_fort25_gate), each failing alone.
TEST(GateComparisons, Fort25GateAssertsItsCounts) {
  Fort25Comparison c;
  c.rows = c.expected_rows = c.row_like_lines = c.compared = c.phase_compared = 2;
  expect_fort25_gate(c, 0, {2, 0});
  EXPECT_NONFATAL_FAILURE(expect_fort25_gate(c, 0, {1, 0}), "phases compared vs the case's pin");
  EXPECT_NONFATAL_FAILURE(expect_fort25_gate(c, 1, {2, 0}), "referred to Zana");
  Fort25Comparison failed = c;
  failed.failed = 1;
  EXPECT_NONFATAL_FAILURE(expect_fort25_gate(failed, 0, {2, 0}), "failed rows");
  Fort25Comparison not_compared = c;
  not_compared.compared = 1;
  EXPECT_NONFATAL_FAILURE(expect_fort25_gate(not_compared, 0, {2, 0}), "rows compared vs rows read");
  Fort25Comparison empty;
  EXPECT_NONFATAL_FAILURE(expect_fort25_gate(empty, 0, {0, 0}), "no rows read");
  // From files: a missing row, and a row that the independent pattern does not count (an E exponent).
  const Fort25Fixture f;
  const Fort25Comparison missing = compare_text(f.run, replaced(f.text(), f.row(1, 1, 0) + "\n", ""));
  EXPECT_NONFATAL_FAILURE(expect_fort25_gate(missing, 0, {1, 0}), "rows read vs the Fortran's print rule");
  const Fort25Comparison e_letter = compare_text(f.run, replaced(f.text(), "0.100000D-01", "0.100000E-01"));
  EXPECT_NONFATAL_FAILURE(expect_fort25_gate(e_letter, 0, {2, 0}), "rows read vs an independent count of the file");
}

TEST(GateComparisons, Fort25RowOutsideABlockFails) {
  const Fort25Fixture f;
  const Fort25Comparison c = compare_text(f.run, f.row(0, 0, 0) + "\n" + f.text());
  EXPECT_EQ(c.failed, 1);
  ASSERT_EQ(c.failures.size(), 1u);
  EXPECT_NE(c.failures[0].find("a row outside a block"), std::string::npos) << c.failures[0];
}

TEST(GateComparisons, Fort25ExtraEnergyBlockFails) {
  const Fort25Fixture f;
  EXPECT_EQ(compare_text(f.run, f.text() + " total energy = -0.42671410D+01\n").failed, 1);
}

TEST(GateComparisons, Fort25EmptyFileFails) {
  const Fort25Fixture f;
  const Fort25Comparison c = compare_text(f.run, "");
  EXPECT_EQ(c.failed, 1);
  EXPECT_EQ(c.rows, 0);
  EXPECT_EQ(c.expected_rows, 2);
}

// A block header with a level the run does not analyse (v = 1 of n_v = 1, or v = -1).
TEST(GateComparisons, Fort25LevelOutsideTheRunFails) {
  const Fort25Fixture f;
  EXPECT_EQ(compare_text(f.run, replaced(f.text(), "   0   0   0   1", "   1   0   0   1")).failed, 1);
  EXPECT_EQ(compare_text(f.run, replaced(f.text(), "   0   0   0   1", "  -1   0   0   1")).failed, 1);
}

// Elements that are all rounding noise (|S| ~ 1e-15): the scale is the floor 1e-3, the code term 1e-13.
TEST(GateComparisons, Fort25CodeTermAtTheFloor) {
  Fort25Fixture f;
  f.s = {cplx(1e-15, 2e-15), cplx(-1e-15, 1e-15)};
  f.set_cpp(f.s);
  const std::string text = f.text();
  f.set_cpp({f.s[0] + cplx(0.9e-13, 0.0), f.s[1]});
  EXPECT_EQ(compare_text(f.run, text).failed, 0);
  f.set_cpp({f.s[0] + cplx(1.1e-13, 0.0), f.s[1]});
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
}

// A printed zero is exact: a C++ element of 1e-9 there fails (half a unit of 0.000000D+00 would allow 5e-7).
TEST(GateComparisons, Fort25PrintedZeroIsExact) {
  Fort25Fixture f;
  f.s[1] = cplx(0.0, 0.0);
  f.set_cpp(f.s);
  const std::string text = f.text();
  EXPECT_EQ(compare_text(f.run, text).failed, 0);
  f.set_cpp({f.s[0], cplx(1e-9, 0.0)});
  EXPECT_EQ(compare_text(f.run, text).failed, 1);
}

// With one diffraction order the Fortran prints one header per v for all rotational states: refused loudly.
TEST(GateComparisons, Fort25SingleDiffractionOrderIsNotComparable) {
  Fort25Fixture f;
  f.run.states[1].order = 0;
  const Fort25Comparison c = compare_text(f.run, f.text());
  EXPECT_EQ(c.failed, 2);  // the refusal, and the row under the order-1 header
  ASSERT_FALSE(c.failures.empty());
  EXPECT_NE(c.failures[0].find("one diffraction order"), std::string::npos) << c.failures[0];
}

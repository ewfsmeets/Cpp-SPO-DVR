// Tests: memory estimate (array sizes from the input) and the byte formatter.
#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

#include "analysis/fine_analysis.hpp"
#include "analysis/flux_analysis.hpp"
#include "app/fine_analysis_run.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"
#include "report/format.hpp"

using namespace spodvr;

namespace {

nlohmann::json reference_json(const std::string& name) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + name + "/input.json");
  return nlohmann::json::parse(file);
}

struct RunGrids {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;

  explicit RunGrids(const nlohmann::json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)) {}

  [[nodiscard]] MemoryEstimate estimate(int n_threads) const {
    return estimate_memory(config, grids, angular, n_threads);
  }
};

constexpr std::size_t real_bytes = 8;
constexpr std::size_t complex_bytes = 16;

}  // namespace

// hardwall_small: 3 x 10 angular points, 4 x 4 (x, y), 113 of 176 Z points in the interaction
// region, 56 r points, 280 specular Z points, 15 rotational states.
TEST(MemoryEstimate, ReferenceCaseByHand) {
  const MemoryEstimate e = RunGrids(reference_json("hardwall_small")).estimate(4);
  ASSERT_EQ(e.items.size(), 5U);

  const std::size_t potential = 3 * 10 * 16 * 113 * 56;
  EXPECT_EQ(e.items[0].values, potential);
  EXPECT_EQ(e.items[0].bytes, potential * real_bytes);
  EXPECT_TRUE(e.items[0].during_setup);
  EXPECT_FALSE(e.items[0].during_propagation);  // freed: exp(-i V dt) is precomputed
  EXPECT_EQ(e.items[1].bytes, potential * complex_bytes);

  const std::size_t scattering = 15 * 16 * 176 * 56;
  EXPECT_EQ(e.items[2].values, scattering);
  EXPECT_EQ(e.items[3].bytes, 2 * 56 * 280 * complex_bytes);

  // T_Z + T_r unit (chunk 4 of n_rot * ny = 60, nr * nz) is larger than the (x, y) group
  // (all 56 r points: 15 * 16 * 56) plus the DVR blocks ((10 + 1) * 3 * 16).
  const std::size_t work = 4 * 56 * 176;
  EXPECT_EQ(e.items[4].values, work);
  EXPECT_EQ(e.items[4].bytes, 4 * work * complex_bytes);
  EXPECT_FALSE(e.items[4].during_setup);

  const std::size_t wave_functions = (scattering + 2 * 56 * 280) * complex_bytes;
  EXPECT_EQ(e.setup_bytes, potential * (real_bytes + complex_bytes) + wave_functions);
  EXPECT_EQ(e.propagation_bytes, potential * complex_bytes + wave_functions + 4 * work * complex_bytes);
  EXPECT_EQ(e.peak_bytes, std::max(e.setup_bytes, e.propagation_bytes));
}

// The estimate agrees with the arrays the program actually allocates.
TEST(MemoryEstimate, MatchesAllocations) {
  const RunGrids setup(reference_json("hardwall_short"));
  const MemoryEstimate e = setup.estimate(1);
  PotentialGrid potential =
      compute_potential_grid(*make_potential(setup.config.input), setup.config, setup.grids, setup.angular);
  EXPECT_EQ(e.items[0].values, potential.values.size());

  const FftwFft fft;
  const GasPhase gas =
      compute_gas_phase(setup.config, setup.grids, setup.angular, potential.asymptotic, fft, LapackEigensolver());
  WavePacketPropagation propagation(setup.config, setup.grids, setup.angular, std::move(potential), gas, fft,
                                    make_cpu_backend, e.flux_slot_count);
  EXPECT_EQ(e.items[2].values, propagation.wave_packet_on_host().scattering.values.size());
  EXPECT_EQ(e.items[3].bytes, 2 * propagation.wave_packet_on_host().specular.values.size() * complex_bytes);
}

// Above 4 GiB the exp(-i V dt) table is not stored and the potential is kept instead.
TEST(MemoryEstimate, LargePotentialComputedOnTheFly) {
  nlohmann::json json = reference_json("hardwall_small");
  json["nX"] = 64;
  json["nY"] = 64;
  const MemoryEstimate e = RunGrids(json).estimate(1);
  EXPECT_GT(e.items[0].values * complex_bytes, PotentialOperator::default_table_limit);
  EXPECT_TRUE(e.items[0].during_propagation);
  EXPECT_EQ(e.items[1].bytes, 0U);
}

TEST(MemoryEstimate, PotentialOnlyRun) {
  nlohmann::json json = reference_json("hardwall_small");
  json["Potential"] = true;
  const MemoryEstimate e = RunGrids(json).estimate(4);
  ASSERT_EQ(e.items.size(), 1U);
  EXPECT_EQ(e.peak_bytes, e.items[0].bytes);
  EXPECT_EQ(e.propagation_bytes, 0U);
}

namespace {

// hardwall_small with the flux analysis at r_c = 2.0 bohr, 0.2 .. 0.6 eV by 0.15 eV: 3 energies (flux Q3).
nlohmann::json hardwall_small_with_flux() {
  nlohmann::json json = reference_json("hardwall_small");
  json.update({{"doFluxAnalysis", true},
               {"rFluxAnalysis", 2.0},
               {"EMinFluxAna", 0.2},
               {"EMaxFluxAna", 0.6},
               {"DeltaEFluxAna", 0.15}});
  return json;
}

}  // namespace

// The flux analysis adds its coefficients and the slots, 2 complex values (psi, d psi / dr) per block
// and plane below Z_asymptotic each, during the time loop only: 15 x 16 blocks, 113 planes, 3
// energies, and 8 slots (the budget holds far more of this small grid).
TEST(MemoryEstimate, FluxAnalysisByHand) {
  const MemoryEstimate without = RunGrids(reference_json("hardwall_small")).estimate(4);
  const MemoryEstimate e = RunGrids(hardwall_small_with_flux()).estimate(4);
  EXPECT_EQ(without.flux_slot_count, 0);
  ASSERT_EQ(e.items.size(), 7U);
  const std::size_t surface = 2 * 15 * 16 * 113;
  EXPECT_EQ(e.items[5].values, surface * 3);
  EXPECT_EQ(e.items[5].bytes, surface * 3 * complex_bytes);
  EXPECT_EQ(e.flux_slot_count, max_flux_slots);
  EXPECT_EQ(e.items[6].values, surface * max_flux_slots);
  EXPECT_EQ(e.items[6].bytes, surface * max_flux_slots * complex_bytes);
  for (const std::size_t i : {5U, 6U}) {
    EXPECT_FALSE(e.items[i].during_setup);
    EXPECT_TRUE(e.items[i].during_propagation);
  }
  EXPECT_EQ(e.setup_bytes, without.setup_bytes);
  EXPECT_EQ(e.propagation_bytes, without.propagation_bytes + e.items[5].bytes + e.items[6].bytes);
}

// The flux arrays of the estimate are those the propagation allocates with the estimate's slot count.
TEST(MemoryEstimate, FluxMatchesAllocations) {
  const RunGrids setup(hardwall_small_with_flux());
  const MemoryEstimate e = setup.estimate(1);
  PotentialGrid potential =
      compute_potential_grid(*make_potential(setup.config.input), setup.config, setup.grids, setup.angular);
  const FftwFft fft;
  const GasPhase gas =
      compute_gas_phase(setup.config, setup.grids, setup.angular, potential.asymptotic, fft, LapackEigensolver());
  WavePacketPropagation propagation(setup.config, setup.grids, setup.angular, std::move(potential), gas, fft,
                                    make_cpu_backend, e.flux_slot_count);
  ASSERT_TRUE(propagation.flux().has_value());
  const FluxCoefficients& c = propagation.flux()->coefficients();
  EXPECT_EQ(e.items[5].values, c.phi.size() + c.phi_derivative.size());
  EXPECT_EQ(propagation.flux()->slot_count(), e.flux_slot_count);
  EXPECT_EQ(e.items[6].values, std::size_t(e.flux_slot_count) * 2 * c.n_blocks() * c.n_planes());
}

// The slot budget on the DZ grid (inputs/thesis_fig5_3_dz, which must exist, with r_c = 2.5 bohr and
// 0.4 .. 0.95 eV by 0.025 eV): 435 x 144 blocks, 125 planes, so one slot is 250.56 MB, the coefficients
// of 22 energies 22 slots, and the budget holds 8 slots (the comment at flux_slot_budget_bytes).
TEST(MemoryEstimate, FluxSlotsOfTheDzGrid) {
  const std::string path = std::string(SPODVR_REFERENCE_DIR) + "/../../inputs/thesis_fig5_3_dz/input.json";
  std::ifstream file(path);
  ASSERT_TRUE(file.is_open()) << path;
  nlohmann::json json = nlohmann::json::parse(file);
  json.update({{"doFluxAnalysis", true},
               {"rFluxAnalysis", 2.5},
               {"EMinFluxAna", 0.4},
               {"EMaxFluxAna", 0.95},
               {"DeltaEFluxAna", 0.025}});
  const MemoryEstimate e = RunGrids(json).estimate(1);
  const std::size_t slot_bytes = 2 * std::size_t(435 * 144) * 125 * complex_bytes;
  EXPECT_EQ(slot_bytes, 250'560'000U);
  ASSERT_GE(e.items.size(), 2U);
  const MemoryItem& coefficients = e.items[e.items.size() - 2];
  const MemoryItem& slots = e.items.back();
  EXPECT_EQ(coefficients.bytes, 22 * slot_bytes);
  EXPECT_EQ(e.flux_slot_count, 8);
  EXPECT_EQ(slots.bytes, 8 * slot_bytes);
  EXPECT_LE(slots.bytes, flux_slot_budget_bytes);
}

// The slot count: what the budget holds, 1 to max_flux_slots, without a known physical memory or when
// the run is too large without the flux analysis too; within the physical memory as many as fit next
// to the base and the coefficients, and "does not fit" only when not even one slot does.
TEST(PlanFluxMemory, SlotsFromTheBudgetAndThePhysicalMemory) {
  constexpr std::size_t gib = std::size_t(1) << 30;
  constexpr std::size_t mib = std::size_t(1) << 20;
  auto slots = [](std::size_t slot_bytes, std::optional<std::size_t> physical = std::nullopt) {
    return plan_flux_memory(10 * gib, 4 * gib, slot_bytes, physical);
  };
  EXPECT_EQ(slots(1).slot_count, max_flux_slots);
  EXPECT_EQ(slots(250 * mib).slot_count, 8);
  EXPECT_EQ(slots(gib).slot_count, 2);
  EXPECT_EQ(slots(3 * gib).slot_count, 1);  // above the budget: still one slot
  EXPECT_TRUE(slots(3 * gib).fits);

  EXPECT_EQ(slots(gib, 20 * gib).slot_count, 2);                      // the budget binds
  EXPECT_EQ(slots(250 * mib, 15 * gib).slot_count, 4);                // 1 GiB free for slots
  const FluxMemoryPlan one = slots(250 * mib, 14 * gib + 250 * mib);  // exactly one slot
  EXPECT_EQ(one.slot_count, 1);
  EXPECT_TRUE(one.fits);
  const FluxMemoryPlan none = slots(250 * mib, 14 * gib + 249 * mib);
  EXPECT_EQ(none.slot_count, 1);
  EXPECT_FALSE(none.fits);
  EXPECT_FALSE(slots(250 * mib, 12 * gib).fits);  // not even the coefficients
  const FluxMemoryPlan too_large_anyway = slots(250 * mib, 9 * gib);
  EXPECT_EQ(too_large_anyway.slot_count, 8);
  EXPECT_TRUE(too_large_anyway.fits);  // the flux analysis is not the cause

  EXPECT_THROW((void)plan_flux_memory(0, 0, 0, std::nullopt), std::invalid_argument);
}

TEST(FormatBytes, Units) {
  EXPECT_EQ(format_bytes(0), "0 B");
  EXPECT_EQ(format_bytes(1023), "1023 B");
  EXPECT_EQ(format_bytes(1024), "1.0 KiB");
  EXPECT_EQ(format_bytes(std::size_t(3) << 19), "1.5 MiB");
  EXPECT_EQ(format_bytes(std::size_t(4) << 30), "4.0 GiB");
  EXPECT_EQ(format_bytes(std::size_t(2048) << 30), "2.0 TiB");
}

// Fine analysis (hardwall_small, 1400 records: ndused 1386, 65 grid points; 240 states x 2 v).
TEST(MemoryEstimate, FineAnalysisByHandAndAllocations) {
  nlohmann::json json = reference_json("hardwall_small");
  json["onlyAnalysis"] = true;
  json["RestartDir"] = "unused";
  const RunGrids setup(json);
  const FftwFft fft;
  const GasPhase gas =
      compute_gas_phase(setup.config, setup.grids, setup.angular,
                        compute_asymptotic_potential(*make_potential(setup.config.input), setup.config, setup.grids),
                        fft, LapackEigensolver());
  const AnalysisSetup analysis(setup.config, setup.angular, gas);
  const FineAnalysisGrid grid = make_fine_analysis_grid(setup.config, 1400);
  const MemoryEstimate e = estimate_fine_analysis_memory(grid, analysis.basis);
  ASSERT_EQ(e.items.size(), 4U);
  const std::size_t n_series = 240 * 2;
  EXPECT_EQ(e.items[0].values, n_series * 1387);
  EXPECT_EQ(e.items[1].values, n_series * 65);
  EXPECT_EQ(e.items[3].values, std::size_t(2772) * 64);
  EXPECT_EQ(e.peak_bytes, e.items[0].bytes + e.items[1].bytes + e.items[2].bytes + e.items[3].bytes);
  EXPECT_EQ(e.propagation_bytes, 0U);

  // Against the arrays that are allocated.
  const BkSeries series(grid, analysis.basis, 0.0);
  EXPECT_EQ(e.items[0].values, std::size_t(series.n_series()) * series.length());
}

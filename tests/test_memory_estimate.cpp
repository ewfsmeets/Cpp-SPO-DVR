// Tests: memory estimate (array sizes from the input) and the byte formatter.
#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

#include "analysis/fine_analysis.hpp"
#include "analysis/flux_analysis.hpp"
#include "app/fine_analysis_run.hpp"
#include "core/error.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"
#include "report/flux_report.hpp"
#include "report/format.hpp"
#include "report/sections.hpp"

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
                                    make_cpu_backend, e.flux_slot_count());
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
  EXPECT_EQ(without.flux_slot_count(), 0);
  ASSERT_EQ(e.items.size(), 7U);
  const std::size_t surface = 2 * 15 * 16 * 113;
  EXPECT_EQ(e.items[5].values, surface * 3);
  EXPECT_EQ(e.items[5].bytes, surface * 3 * complex_bytes);
  EXPECT_EQ(e.flux_slot_count(), max_flux_slots);
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
                                    make_cpu_backend, e.flux_slot_count());
  ASSERT_TRUE(propagation.flux().has_value());
  const FluxCoefficients& c = propagation.flux()->coefficients();
  EXPECT_EQ(e.items[5].values, c.phi.size() + c.phi_derivative.size());
  EXPECT_EQ(propagation.flux()->slot_count(), e.flux_slot_count());
  EXPECT_EQ(e.items[6].values, std::size_t(e.flux_slot_count()) * 2 * c.n_blocks() * c.n_planes());
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
  EXPECT_EQ(e.flux_slot_count(), 8);
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

namespace {

constexpr std::size_t gib = std::size_t(1) << 30;
constexpr std::size_t mib = std::size_t(1) << 20;

// A run that needs 10 GiB of host memory during the set-up and 12 GiB in the time loop on the CPU, and
// 6 GiB and 9 GiB as a GPU run; a flux analysis of 4 energies of 1 GiB each and slots of 250 MiB.
MemoryEstimate made_up_estimate(HostMemoryBase gpu_host = {6 * gib, 9 * gib}) {
  MemoryEstimate e;
  e.setup_bytes = 10 * gib;
  e.propagation_bytes = 12 * gib;
  e.peak_bytes = 12 * gib;
  e.gpu_host = gpu_host;
  return e;
}
const FluxMemoryItems made_up_flux{gib, 4, 250 * mib};

FluxMemoryCheck cpu_check(std::optional<std::size_t> physical) {
  return check_flux_memory(made_up_estimate(), made_up_flux, physical, Device::Cpu);
}

// The message of the InputError of require_flux_memory, or "" if there is none; fails on another error.
std::string flux_memory_error(const FluxMemoryCheck& check) {
  try {
    require_flux_memory(check);
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), "DeltaEFluxAna");
    return e.what();
  }
  return "";
}

}  // namespace

// CPU: the slots shrink to what the physical memory leaves next to the time loop's 12 GiB and the
// coefficients' 4 GiB (the budget's 8, then 2, then 1), and the flux analysis "does not fit", an input
// error naming DeltaEFluxAna, only when not even one slot fits. The most energies that fit with one slot
// is reported in each case.
TEST(FluxMemoryRule, OnTheCpuTheSlotsShrinkAndOnlyTheFluxAnalysisCausesTheError) {
  const FluxMemoryCheck roomy = cpu_check(32 * gib);
  EXPECT_EQ(roomy.base.setup, 10 * gib);
  EXPECT_EQ(roomy.base.time_loop, 12 * gib);
  EXPECT_TRUE(roomy.base_fits);
  EXPECT_EQ(roomy.plan.slot_count, 8);
  EXPECT_TRUE(roomy.plan.fits);
  EXPECT_EQ(roomy.max_energies, 19U);  // (20 GiB - 250 MiB) / 1 GiB
  EXPECT_EQ(flux_memory_error(roomy), "");

  const FluxMemoryCheck two = cpu_check(16 * gib + 600 * mib);
  EXPECT_EQ(two.plan.slot_count, 2);
  EXPECT_TRUE(two.plan.fits);
  const FluxMemoryCheck one = cpu_check(16 * gib + 250 * mib);
  EXPECT_EQ(one.plan.slot_count, 1);
  EXPECT_TRUE(one.plan.fits);
  EXPECT_EQ(one.max_energies, 4U);
  EXPECT_EQ(flux_memory_error(one), "");

  const FluxMemoryCheck too_small = cpu_check(16 * gib + 249 * mib);
  EXPECT_TRUE(too_small.base_fits);
  EXPECT_FALSE(too_small.plan.fits);
  EXPECT_EQ(too_small.plan.slot_count, 1);
  EXPECT_EQ(too_small.max_energies, 3U);
  const std::string error = flux_memory_error(too_small);
  EXPECT_NE(error.find("DeltaEFluxAna: the flux analysis needs 4.2 GiB"), std::string::npos) << error;
  EXPECT_NE(error.find("at most 3 energies fit"), std::string::npos) << error;
}

// The flux analysis is not the cause when the run does not fit without it: in the time loop (11 GiB for
// 12 GiB) or only during the set-up (a 20 GiB set-up in 18 GiB, where the time loop with the flux
// analysis would fit). Then the slots are the budget's and there is no error; so with an unknown
// physical memory. Through estimate_memory and print_memory_estimate (hardwall_small with the flux
// analysis, the physical memory 1 byte below the set-up of a CPU run or of a GPU run's host memory,
// with room for the time loop and the flux arrays): the estimate's warning is printed, and no error.
TEST(FluxMemoryRule, NoErrorWhenTheRunIsTooLargeWithoutItOrTheMemoryIsUnknown) {
  for (const FluxMemoryCheck& check :
       {cpu_check(11 * gib), cpu_check(std::nullopt),
        check_flux_memory(made_up_estimate({20 * gib, 9 * gib}), made_up_flux, 18 * gib, Device::Cuda)}) {
    EXPECT_FALSE(check.base_fits);
    EXPECT_EQ(check.plan.slot_count, 8);
    EXPECT_TRUE(check.plan.fits);
    EXPECT_FALSE(check.max_energies.has_value());
    EXPECT_EQ(flux_memory_error(check), "");
  }

  const RunGrids setup(hardwall_small_with_flux());
  for (const Device device : {Device::Cpu, Device::Cuda}) {
    SCOPED_TRACE(to_string(device));
    auto estimate = [&](std::optional<std::size_t> physical) {
      return estimate_memory(setup.config, setup.grids, setup.angular, make_grid_symmetry(setup.config, setup.angular),
                             4, device, physical);
    };
    const FluxMemoryCheck unlimited = estimate(std::nullopt).flux_memory.value();
    const std::size_t flux_bytes =
        unlimited.items.coefficient_bytes() + unlimited.plan.slot_count * unlimited.items.slot_bytes;
    const std::size_t physical = unlimited.base.setup - 1;
    ASSERT_LT(unlimited.base.time_loop + flux_bytes, physical);
    const MemoryEstimate e = estimate(physical);
    EXPECT_FALSE(e.flux_memory->base_fits);
    EXPECT_EQ(e.flux_slot_count(), max_flux_slots);
    EXPECT_EQ(flux_memory_error(*e.flux_memory), "");
    std::ostringstream out;
    print_memory_estimate(out, Style{false}, e, physical);
    EXPECT_NE(out.str().find("Warning: the run needs"), std::string::npos) << out.str();
  }
}

// --device cuda: the base is the GPU run's host memory (MemoryEstimate::gpu_host), the larger of its set-up
// and time loop, and the flux analysis comes next to its time loop. In 14 GiB + 250 MiB a CPU run's flux
// analysis does not fit (12 GiB + 4 GiB), a GPU run's does with 5 slots (9 GiB + 4 GiB + 5 x 250 MiB).
TEST(FluxMemoryRule, OnTheGpuTheBaseIsTheHostMemoryOfAGpuRun) {
  const std::size_t physical = 14 * gib + 250 * mib;
  EXPECT_FALSE(cpu_check(physical).plan.fits);
  const FluxMemoryCheck gpu = check_flux_memory(made_up_estimate(), made_up_flux, physical, Device::Cuda);
  EXPECT_EQ(gpu.device, Device::Cuda);
  EXPECT_EQ(gpu.base.setup, 6 * gib);
  EXPECT_EQ(gpu.base.time_loop, 9 * gib);
  EXPECT_TRUE(gpu.base_fits);
  EXPECT_TRUE(gpu.plan.fits);
  EXPECT_EQ(gpu.plan.slot_count, 5);
  EXPECT_EQ(gpu.max_energies, 5U);  // (5 GiB + 250 MiB - 250 MiB) / 1 GiB
  const FluxMemoryCheck gpu_short =
      check_flux_memory(made_up_estimate({6 * gib, 11 * gib}), made_up_flux, physical, Device::Cuda);
  EXPECT_FALSE(gpu_short.plan.fits);
  EXPECT_NE(flux_memory_error(gpu_short).find("host memory of a GPU run"), std::string::npos);
}

// max_energies is the largest count that fits: that many energies fit with one slot and one more do
// not, for physical memories in steps of 37 MiB from what the base alone needs.
TEST(FluxMemoryRule, MaxEnergiesIsTheLargestCountThatFits) {
  for (std::size_t physical = 12 * gib; physical < 20 * gib; physical += 37 * mib) {
    const std::optional<std::size_t> n = cpu_check(physical).max_energies;
    ASSERT_TRUE(n.has_value());
    auto fits = [&](std::size_t energies) {
      const FluxMemoryItems items{made_up_flux.bytes_per_energy, static_cast<int>(energies), made_up_flux.slot_bytes};
      return check_flux_memory(made_up_estimate(), items, physical, Device::Cpu).plan.fits;
    };
    if (*n > 0) {
      EXPECT_TRUE(fits(*n)) << physical;
    }
    EXPECT_FALSE(fits(*n + 1)) << physical;
  }
}

// The printout below the memory estimate shows the numbers of the error: the time loop's base, the
// physical memory, the size of an energy and of a slot, and the most energies that fit.
TEST(FluxMemoryRule, PrintoutAndErrorShowTheSameNumbers) {
  const FluxMemoryCheck check = cpu_check(16 * gib + 249 * mib);
  std::ostringstream out;
  print_flux_memory(out, Style{false}, check);
  const std::string printed = out.str();
  const std::string error = flux_memory_error(check);
  const std::string base = format_bytes(12 * gib);
  const std::string physical = format_bytes(16 * gib + 249 * mib);
  const std::string energy = format_bytes(gib);
  const std::string slot = format_bytes(250 * mib);
  // The same number in the phrase of each text: (printout, error).
  for (const auto& [in_printout, in_error] : {std::pair{"/ " + base, "minus " + base + " for the rest of the run"},
                                              {"Physical memory ", "(" + physical + " minus"},
                                              {"4 energies x " + energy, "4 energies of " + energy},
                                              {"1 x " + slot, "one slot of " + slot},
                                              {"3 (with 1 slot)", "at most 3 energies fit"}}) {
    EXPECT_NE(printed.find(in_printout), std::string::npos) << in_printout << " in\n" << printed;
    EXPECT_NE(error.find(in_error), std::string::npos) << in_error << " in\n" << error;
  }
  std::istringstream lines(printed);
  bool physical_line = false;
  for (std::string line; std::getline(lines, line);)
    if (line.find("Physical memory") != std::string::npos) physical_line = line.ends_with(" " + physical);
  EXPECT_TRUE(physical_line) << printed;
  EXPECT_NE(printed.find("does not fit"), std::string::npos) << printed;
}

// The host memory of a GPU run of hardwall_small (the numbers of ReferenceCaseByHand): V and exp(i A)
// with the two specular waves during the set-up; psi and a third specular wave in the time loop; the
// read-back plane in both
// (15 rotational functions x 16 (x, y) x 56 r) in the time loop.
TEST(MemoryEstimate, GpuHostMemoryByHand) {
  const MemoryEstimate e = RunGrids(reference_json("hardwall_small")).estimate(4);
  const std::size_t potential = 3 * 10 * 16 * 113 * 56;
  const std::size_t specular = 56 * 280 * complex_bytes;
  const std::size_t plane = 15 * 16 * 56 * complex_bytes;
  EXPECT_EQ(e.gpu_host.setup, potential * (real_bytes + complex_bytes) + 2 * specular + plane);
  EXPECT_EQ(e.gpu_host.time_loop, 15 * 16 * 176 * 56 * complex_bytes + 3 * specular + plane);
}

// The memory warning of a GPU run compares its own peak host memory (gpu_host_peak_bytes: the larger of
// the set-up and the time loop with the flux arrays) with the physical memory, not a CPU run's peak:
// with the physical memory between the two, a GPU run gets no warning and a CPU run does.
TEST(MemoryEstimate, GpuRunWarnsOnItsOwnHostMemory) {
  const RunGrids setup(hardwall_small_with_flux());
  auto estimate = [&](Device device, std::optional<std::size_t> physical) {
    return estimate_memory(setup.config, setup.grids, setup.angular, make_grid_symmetry(setup.config, setup.angular), 4,
                           device, physical);
  };
  const MemoryEstimate gpu = estimate(Device::Cuda, std::nullopt);
  ASSERT_TRUE(gpu.gpu_host_peak_bytes.has_value());
  const FluxMemoryCheck& flux = gpu.flux_memory.value();
  const std::size_t flux_bytes = flux.items.coefficient_bytes() + flux.plan.slot_count * flux.items.slot_bytes;
  EXPECT_EQ(*gpu.gpu_host_peak_bytes, std::max(gpu.gpu_host.setup, gpu.gpu_host.time_loop + flux_bytes));
  EXPECT_FALSE(estimate(Device::Cpu, std::nullopt).gpu_host_peak_bytes.has_value());

  const std::size_t physical = (*gpu.gpu_host_peak_bytes + gpu.peak_bytes) / 2;
  ASSERT_LT(*gpu.gpu_host_peak_bytes, physical);
  ASSERT_GT(gpu.peak_bytes, physical);
  auto printout = [&](Device device) {
    std::ostringstream out;
    print_memory_estimate(out, Style{false}, estimate(device, physical), physical);
    return out.str();
  };
  const std::string gpu_text = printout(Device::Cuda);
  EXPECT_EQ(gpu_text.find("Warning: the run needs"), std::string::npos) << gpu_text;
  EXPECT_NE(gpu_text.find("Host memory of the GPU run (peak)"), std::string::npos) << gpu_text;
  EXPECT_NE(printout(Device::Cpu).find("Warning: the run needs " + format_bytes(gpu.peak_bytes)), std::string::npos);
}

// estimate_memory applies the rule with the physical memory: with 64 threads the time loop of
// hardwall_small with the flux analysis needs more than its set-up, and the slots shrink to the 3 that
// fit next to it and the coefficients; 1 byte less than the coefficients and one slot: does not fit.
TEST(MemoryEstimate, FluxSlotsFollowThePhysicalMemory) {
  const RunGrids setup(hardwall_small_with_flux());
  auto estimate = [&](std::optional<std::size_t> physical) {
    return estimate_memory(setup.config, setup.grids, setup.angular, make_grid_symmetry(setup.config, setup.angular),
                           64, Device::Cpu, physical);
  };
  const MemoryEstimate unlimited = estimate(std::nullopt);
  ASSERT_TRUE(unlimited.flux_memory.has_value());
  const FluxMemoryItems flux = unlimited.flux_memory->items;
  EXPECT_EQ(flux.n_energies, 3);
  EXPECT_EQ(flux.slot_bytes, 2 * std::size_t(15 * 16) * 113 * complex_bytes);
  EXPECT_EQ(flux.coefficient_bytes(), unlimited.items[5].bytes);
  const std::size_t time_loop = unlimited.flux_memory->base.time_loop;
  EXPECT_EQ(time_loop, unlimited.propagation_bytes - unlimited.items[5].bytes - unlimited.items[6].bytes);
  ASSERT_GT(time_loop, unlimited.setup_bytes);

  const MemoryEstimate three = estimate(time_loop + flux.coefficient_bytes() + 3 * flux.slot_bytes);
  EXPECT_EQ(three.flux_slot_count(), 3);
  EXPECT_EQ(three.items[6].bytes, 3 * flux.slot_bytes);
  EXPECT_TRUE(three.flux_memory->plan.fits);
  EXPECT_NO_THROW(require_flux_memory(*three.flux_memory));

  const MemoryEstimate none = estimate(time_loop + flux.coefficient_bytes() + flux.slot_bytes - 1);
  EXPECT_FALSE(none.flux_memory->plan.fits);
  EXPECT_EQ(none.flux_slot_count(), 1);
  EXPECT_THROW(require_flux_memory(*none.flux_memory), InputError);
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

// Tests: memory estimate (array sizes from the input) and the byte formatter.
#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

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
        angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear, config.input.j_in,
                config.input.mj_in) {}

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
      compute_potential_grid(*make_potential(setup.config.input.system_id), setup.config, setup.grids, setup.angular);
  EXPECT_EQ(e.items[0].values, potential.values.size());

  const FftwFft fft;
  const GasPhase gas =
      compute_gas_phase(setup.config, setup.grids, setup.angular, potential.asymptotic, fft, LapackEigensolver());
  const WavePacketPropagation propagation(setup.config, setup.grids, setup.angular, std::move(potential), gas, fft);
  EXPECT_EQ(e.items[2].values, propagation.wave_packet().scattering.values.size());
  EXPECT_EQ(e.items[3].bytes, 2 * propagation.wave_packet().specular.values.size() * complex_bytes);
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

TEST(FormatBytes, Units) {
  EXPECT_EQ(format_bytes(0), "0 B");
  EXPECT_EQ(format_bytes(1023), "1023 B");
  EXPECT_EQ(format_bytes(1024), "1.0 KiB");
  EXPECT_EQ(format_bytes(std::size_t(3) << 19), "1.5 MiB");
  EXPECT_EQ(format_bytes(std::size_t(4) << 30), "4.0 GiB");
  EXPECT_EQ(format_bytes(std::size_t(2048) << 30), "2.0 TiB");
}

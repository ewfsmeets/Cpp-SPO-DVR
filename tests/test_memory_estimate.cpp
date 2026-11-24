// Tests: memory estimate (array sizes from the input) and the byte formatter.
#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/fine_analysis.hpp"
#include "app/fine_analysis_run.hpp"
#include "designer_cases.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"
#include "propagation/propagation_tables.hpp"
#include "report/format.hpp"

using namespace spodvr;

namespace {

nlohmann::json reference_json(const std::string& name) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + name + "/input.json");
  return nlohmann::json::parse(file);
}

struct RunGrids {
  std::unique_ptr<Potential> pes;  // built once, as the program does
  Config config;
  SpatialGrids grids;
  AngularGrid angular;

  explicit RunGrids(const nlohmann::json& json)
      : pes(make_potential(parse_input(json))),
        config(build_config(parse_input(json), *pes)),
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
  PotentialGrid potential = compute_potential_grid(*setup.pes, setup.config, setup.grids, setup.angular);
  EXPECT_EQ(e.items[0].values, potential.values.size());

  const FftwFft fft;
  const GasPhase gas =
      compute_gas_phase(setup.config, setup.grids, setup.angular, potential.asymptotic, fft, LapackEigensolver());
  WavePacketPropagation propagation(setup.config, setup.grids, setup.angular, std::move(potential), gas, fft);
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
                        compute_asymptotic_potential(*setup.pes, setup.config, setup.grids), fft, LapackEigensolver());
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

// The wave-packet designer (tools/wave_packet_designer) copies the grid indices, the angular sizes, the
// symmetry's use in the estimate and the CPU memory estimate; its expected_cases.json holds its results. When
// this fails, update the designer's copy and regenerate the file (designer_cases.hpp).
// The estimate's rule for the T_xy work array is the operator's: hardwall_small on 64 x 64 (x, y) points (one r
// point per T_xy group: 15 x 4096 values of 16 bytes fit 1 MiB once), without and with parallel momentum. Without
// the shift T_xy runs in place and the work item is max((r, Z) unit, DVR blocks); with it the group array stays.
TEST(MemoryEstimate, InPlaceRuleIsTheOperators) {
  for (const double e_parallel : {0.0, 0.05}) {
    nlohmann::json json = reference_json("hardwall_small");
    json["nX"] = 64;
    json["nY"] = 64;
    json["EParallel"] = e_parallel;
    json["PhiParallel"] = 30.0;
    const RunGrids setup(json);
    const InputParameters& p = setup.config.input;
    const RotationalBasis basis(setup.angular, p.j_in, p.mj_in);
    const ScatteringShape shape{basis.size(), p.nx, p.ny, p.nr, p.nz};
    const PropagationTables tables =
        make_propagation_tables(setup.config, setup.grids, basis, std::vector<double>(p.nr, 0.0));
    const bool shift = xy_momentum_shift_needed(momentum_shift_forward(setup.config));
    EXPECT_EQ(shift, e_parallel != 0.0);
    EXPECT_EQ(shift, !std::all_of(tables.shift_forward.begin(), tables.shift_forward.end(),
                                  [](std::complex<double> c) { return c == std::complex<double>(1.0, 0.0); }));
    const FftwFft fft;
    const ScatteringKineticOperator kinetic(shape, tables, fft);
    ASSERT_EQ(kinetic.group_points(), 1);
    EXPECT_EQ(kinetic.transforms_in_place(), xy_in_place_possible(shape, xy_group_points(shape), shift));
    EXPECT_EQ(kinetic.transforms_in_place(), !shift);
    const std::size_t group = kinetic.transforms_in_place() ? 0 : std::size_t(shape.n_rot) * p.nx * p.ny;
    EXPECT_EQ(kinetic.group_work_size(), group);
    const std::size_t zr_unit = std::size_t(zr_chunk_size(shape.n_rot * shape.ny)) * p.nr * p.nz;
    const std::size_t dvr_blocks =
        std::size_t(AngularTransform::dvr_block_size(setup.angular.n_theta(), setup.angular.n_phi())) *
        AngularTransform::columns_per_pass(p.nx * p.ny, setup.angular.n_theta(), setup.angular.n_phi());
    const MemoryEstimate e = setup.estimate(1);
    const auto work =
        std::find_if(e.items.begin(), e.items.end(), [](const MemoryItem& i) { return i.name == "Work arrays"; });
    ASSERT_NE(work, e.items.end());
    EXPECT_EQ(work->values, std::max(zr_unit, group + dvr_blocks)) << "EParallel " << e_parallel;
  }
}

TEST(MemoryEstimate, MatchesTheWavePacketDesignerCases) {
  using namespace spodvr::testing;
  for (const nlohmann::json& designer_case : designer_cases()) {
    SCOPED_TRACE("designer case " + designer_case.at("name").get<std::string>());
    const RunGrids run(designer_case_input(designer_case));
    const nlohmann::json& grids = designer_case.at("grids");
    EXPECT_EQ(run.grids.iz_analysis, grids.at("iz_analysis").get<int>());
    EXPECT_EQ(run.grids.nz_potential, grids.at("nz_potential").get<int>());
    expect_designer_range(run.grids.z_optical, grids.at("z_optical"), "z_optical");
    expect_designer_range(run.grids.z_specular_optical, grids.at("z_specular_optical"), "z_specular_optical");
    expect_designer_range(run.grids.r_optical, grids.at("r_optical"), "r_optical");
    expect_designer_range(run.grids.z_optical_low, grids.at("z_optical_low"), "z_optical_low");
    expect_designer_range(run.grids.z_specular_optical_low, grids.at("z_specular_optical_low"),
                          "z_specular_optical_low");

    const nlohmann::json& angular = designer_case.at("angular");
    EXPECT_EQ(run.angular.n_theta(), angular.at("n_theta").get<int>());
    EXPECT_EQ(run.angular.n_phi(), angular.at("n_phi").get<int>());
    const InputParameters& p = run.config.input;
    EXPECT_EQ(RotationalBasis(run.angular, p.j_in, p.mj_in).size(), angular.at("rotational_basis_size").get<int>());

    // The designer takes these two counts from the file (the program's SURFACE SYMMETRY block prints them).
    const GridSymmetry symmetry = make_grid_symmetry(run.config, run.angular);
    EXPECT_EQ(symmetry.n_representatives(), designer_case.at("symmetry").at("representatives").get<int>());
    EXPECT_EQ(symmetry.operations().size(), designer_case.at("symmetry").at("operations").get<std::size_t>());

    expect_designer_estimate(run.estimate(designer_case.at("threads").get<int>()), designer_case.at("cpu"));
  }
}

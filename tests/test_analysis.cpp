// Tests: analysis basis (diffraction channels, analysed rotational and rotation-diffraction
// states, channel energies).
#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <set>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/analysis_basis.hpp"
#include "core/error.hpp"
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
#include "propagation/propagation_tables.hpp"

using namespace spodvr;

namespace {

nlohmann::json reference_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return nlohmann::json::parse(file);
}

// Setup up to the analysis basis.
struct AnalysisSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  PotentialGrid potential;
  FftwFft fft;
  GasPhase gas;
  RotationalBasis basis;
  PropagationTables tables;
  AnalysisBasis analysis;

  explicit AnalysisSetup(const nlohmann::json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear, config.input.j_in,
                config.input.mj_in),
        potential(compute_potential_grid(*make_potential(config.input.system_id), config, grids, angular)),
        gas(compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver())),
        basis(angular, config.input.j_in, config.input.mj_in),
        tables(make_propagation_tables(config, grids, basis, potential.asymptotic)),
        analysis(config, angular, basis, gas, tables.xy_kinetic_energies) {}
};

const AnalysisSetup& reference_setup() {
  static const AnalysisSetup setup(reference_json());
  return setup;
}

}  // namespace

TEST(DiffractionOrder, PointGroups) {
  EXPECT_EQ(diffraction_order("C2v", 2, -1), 3);
  EXPECT_EQ(diffraction_order(" c4v ", -2, -2), 4);
  // Hexagonal lattices: n m > 0 counts max(|n|, |m|), otherwise |n| + |m|.
  EXPECT_EQ(diffraction_order("C6v", 2, 1), 2);
  EXPECT_EQ(diffraction_order("C3v", -1, -3), 3);
  EXPECT_EQ(diffraction_order("C3", 2, -1), 3);
  EXPECT_EQ(diffraction_order("C6v", 0, -2), 2);
  EXPECT_THROW((void)diffraction_order("D4h", 0, 0), NotSupportedError);
}

// Fortran stdout of the reference case: "INFO FROM SETTING UP THE DIFFRACTION / ROTATION /
// ROTATION-DIFFRACTION BASIS" and "INFO FROM SETTING UP THE ANALYSIS".
TEST(AnalysisBasis, ReferenceCaseMatchesFortran) {
  const AnalysisBasis& a = reference_setup().analysis;
  EXPECT_EQ(a.diffraction_channels().size(), 16u);  // ndftow
  EXPECT_EQ(a.n_diffraction(), 16);                  // ndifa
  EXPECT_EQ(a.highest_order(), 4);                   // maxord
  EXPECT_EQ(a.max_channels_per_order(), 6);          // maxnum
  EXPECT_EQ(a.max_order(), 4);                       // mdifa
  EXPECT_TRUE(a.diffraction_basis_too_small());      // "Diffractional basis too small" warning
  EXPECT_EQ(a.n_rotation(), 15);                     // nana
  EXPECT_EQ(a.n_states(), 240);                      // nstata
  EXPECT_EQ(a.initial_state() + 1, 97);              // inchna
  EXPECT_EQ(a.j_max(), 4);                           // jana
  EXPECT_EQ(a.n_v(), 2);                             // nviba
}

TEST(AnalysisBasis, DiffractionChannelsAreSortedStablyByOrder) {
  const auto& channels = reference_setup().analysis.diffraction_channels();
  EXPECT_EQ(channels.front().n, 0);
  EXPECT_EQ(channels.front().m, 0);
  for (std::size_t i = 1; i < channels.size(); ++i) {
    EXPECT_LE(channels[i - 1].order, channels[i].order);
    // Within one order the Fortran generation order (n outer, m inner) is kept.
    if (channels[i - 1].order == channels[i].order) {
      EXPECT_LT(std::make_pair(channels[i - 1].n, channels[i - 1].m), std::make_pair(channels[i].n, channels[i].m));
    }
  }
  std::set<std::pair<int, int>> slots;
  for (const DiffractionChannel& c : channels) {
    EXPECT_EQ(c.order, std::abs(c.n) + std::abs(c.m));
    EXPECT_EQ(c.ix, c.n >= 0 ? c.n : 4 + c.n);
    EXPECT_EQ(c.iy, c.m >= 0 ? c.m : 4 + c.m);
    slots.insert({c.ix, c.iy});
  }
  EXPECT_EQ(slots.size(), 16u);  // every FFT slot exactly once
}

TEST(AnalysisBasis, StatesCombineRotationOuterAndDiffractionInner) {
  const AnalysisSetup& s = reference_setup();
  const AnalysisBasis& a = s.analysis;
  for (int is = 0; is < a.n_states(); ++is) {
    const int irot = a.rotation_index(a.state_rotation(is));
    EXPECT_EQ(a.j(is), s.basis.j(irot));
    EXPECT_EQ(a.mj(is), s.basis.mj(irot));
    EXPECT_EQ(&a.channel(is), &a.diffraction_channels()[is % a.n_diffraction()]);
  }
  const int is0 = a.initial_state();
  EXPECT_EQ(a.j(is0), s.config.input.j_in);
  EXPECT_EQ(a.mj(is0), s.config.input.mj_in);
  EXPECT_EQ(a.channel(is0).n, 0);
  EXPECT_EQ(a.channel(is0).m, 0);
}

TEST(AnalysisBasis, ChannelEnergiesAreDiffractionPlusRovibrationalEnergy) {
  const AnalysisSetup& s = reference_setup();
  const AnalysisBasis& a = s.analysis;
  for (int v = 0; v < a.n_v(); ++v)
    for (int is = 0; is < a.n_states(); ++is) {
      const DiffractionChannel& c = a.channel(is);
      const double e_dif = s.tables.xy_kinetic_energies[c.iy + s.config.input.ny * c.ix];
      EXPECT_EQ(a.channel_energy(is, v), e_dif + s.gas.states.energy(v, a.j(is)));
    }
  EXPECT_EQ(a.diffraction_energy(0), s.tables.parallel_energy);
}

TEST(AnalysisBasis, RegularOrderIsJThenMjAscending) {
  const AnalysisBasis& a = reference_setup().analysis;
  const std::vector<int>& rotations = a.regular_rotation_order();
  ASSERT_EQ(static_cast<int>(rotations.size()), a.n_rotation());
  std::vector<std::pair<int, int>> labels;
  for (const int irotana : rotations) labels.push_back({a.j(irotana * a.n_diffraction()), a.mj(irotana * a.n_diffraction())});
  EXPECT_TRUE(std::is_sorted(labels.begin(), labels.end()));
  EXPECT_EQ(labels.front(), std::make_pair(0, 0));
  EXPECT_EQ(labels.back(), std::make_pair(4, 4));

  const std::vector<int> states = a.regular_state_order();
  ASSERT_EQ(static_cast<int>(states.size()), a.n_states());
  EXPECT_EQ(std::set<int>(states.begin(), states.end()).size(), states.size());  // a permutation
  for (std::size_t k = 0; k < states.size(); ++k)
    EXPECT_EQ(a.state_diffraction(states[k]), static_cast<int>(k) % a.n_diffraction());
}

TEST(AnalysisBasis, DifOrdAnaLimitsTheAnalysedChannels) {
  nlohmann::json json = reference_json();
  json["DifOrdAna"] = 1;
  const AnalysisSetup s(json);
  EXPECT_EQ(s.analysis.max_order(), 1);
  EXPECT_EQ(s.analysis.n_diffraction(), 5);  // (0,0), (0,±1), (±1,0)
  EXPECT_EQ(s.analysis.n_states(), 5 * 15);
  EXPECT_EQ(s.analysis.initial_state(), 6 * 5);
}

// Tests: initial wave packet (CoherentStateFunc incl. the tolsp check) and the analysis surface.
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>

#include <nlohmann/json.hpp>

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
#include "propagation/analysis_surface.hpp"
#include "propagation/initial_wave_packet.hpp"

using namespace spodvr;

namespace {

nlohmann::json reference_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return nlohmann::json::parse(file);
}

// Setup up to and including the gas-phase states.
struct GasPhaseSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  PotentialGrid potential;
  FftwFft fft;
  GasPhase gas;

  explicit GasPhaseSetup(const nlohmann::json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear, config.input.j_in,
                config.input.mj_in),
        potential(compute_potential_grid(*make_potential(config.input.system_id), config, grids, angular)),
        gas(compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver())) {}
};

const GasPhaseSetup& reference_setup() {
  static const GasPhaseSetup setup(reference_json());
  return setup;
}

double specular_norm(const SpecularWave& psi, double volume) {
  double sum = 0.0;
  for (const auto& c : psi.values) sum += std::norm(c);
  return sum * volume;
}

}  // namespace

TEST(InitialWavePacket, NormalisedGaussianTimesRovibrationalState) {
  const GasPhaseSetup& s = reference_setup();
  const InitialWavePacket wp = make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);
  EXPECT_NEAR(wp.info.norm, 1.0, 1e-12);
  EXPECT_GT(wp.info.norm_inside, 1.0 - s.config.input.tol_specular);
  // The momentum distribution peaks away from k = 0 and vanishes at the grid edge.
  EXPECT_GT(wp.info.k_space_maximum, wp.info.k_space_at_zero);
  EXPECT_LT(wp.info.k_space_at_edge, 1e-10);

  // Boundaries: first point at Zana, last point before the specular absorber.
  EXPECT_EQ(s.grids.z_specular[wp.info.i_low], s.grids.z_analysis);
  EXPECT_EQ(wp.info.i_high, s.grids.z_specular_optical->first - 1);

  // psi(Z, r) = g(Z) phi(r) with phi normalised on the r grid.
  const auto phi = s.gas.states.wave_function(0, 0);
  const int iz = 170;
  const std::complex<double> ratio = wp.specular(10, iz) / phi[10];
  for (int ir = 0; ir < s.config.input.nr; ir += 5)
    EXPECT_NEAR(std::abs(wp.specular(ir, iz) - ratio * phi[ir]), 0.0, 1e-15);
  EXPECT_NEAR(specular_norm(wp.specular, s.config.input.dz * s.config.input.dr), 1.0, 1e-6);
}

TEST(InitialWavePacket, TooMuchNormOutsideTheAllowedRegionIsAnError) {
  // A narrow energy range gives a wide packet that does not fit between Zana and the absorber.
  nlohmann::json json = reference_json();
  json["EMinWP"] = 0.39;
  json["EMaxWP"] = 0.41;
  const GasPhaseSetup s(json);
  try {
    (void)make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);
    FAIL() << "expected InputError";
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), "TolSp");
  }
}

TEST(AnalysisSurface, DynamicSurfaceMovesInwardsAndFixes) {
  const GasPhaseSetup& s = reference_setup();
  const InitialWavePacket wp = make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);
  AnalysisSurface surface(s.config, s.grids, wp.specular);
  const int before_absorber = s.grids.z_specular_optical->first - 1;
  EXPECT_FALSE(surface.fixed());
  EXPECT_EQ(surface.index(), before_absorber);
  EXPECT_GE(surface.tolerance(), s.config.input.tol_specular * (1.0 - 1e-14));

  // Norm concentrated at one row: the surface moves just beyond it.
  SpecularWave psi(s.config.input.nr, s.grids.z_specular.size());
  psi(3, 150) = 1.0;
  EXPECT_FALSE(surface.update(psi).has_value());
  EXPECT_EQ(surface.index(), 151);
  EXPECT_FALSE(surface.fixed());

  // The norm moving back out fixes the surface where it was found.
  psi(3, 150) = 0.0;
  psi(3, 160) = 1.0;
  const auto fixed_at = surface.update(psi);
  ASSERT_TRUE(fixed_at.has_value());
  EXPECT_TRUE(surface.fixed());
  EXPECT_EQ(surface.index(), 161);
  EXPECT_FALSE(surface.update(psi).has_value());  // no-op once fixed
}

TEST(AnalysisSurface, FixesAtTheAnalysisPointWhenNoNormIsLeft) {
  const GasPhaseSetup& s = reference_setup();
  const InitialWavePacket wp = make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);
  AnalysisSurface surface(s.config, s.grids, wp.specular);
  const SpecularWave empty(s.config.input.nr, s.grids.z_specular.size());
  const auto fixed_at = surface.update(empty);
  ASSERT_TRUE(fixed_at.has_value());
  EXPECT_EQ(surface.index(), s.grids.iz_analysis);
  EXPECT_EQ(*fixed_at, s.grids.z_specular[s.grids.iz_analysis]);
}

TEST(AnalysisSurface, FixedByInput) {
  nlohmann::json json = reference_json();
  json["DynamZSpAna"] = false;
  const GasPhaseSetup s(json);
  const InitialWavePacket wp = make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);
  const AnalysisSurface surface(s.config, s.grids, wp.specular);
  EXPECT_TRUE(surface.fixed());
  EXPECT_EQ(surface.index(), s.grids.z_specular_optical->first);  // start of the specular absorber
}

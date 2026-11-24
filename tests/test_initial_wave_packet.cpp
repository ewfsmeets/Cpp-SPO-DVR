// Tests: initial wave packet (CoherentStateFunc incl. the tolsp check) and the analysis surface.
#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <utility>

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
#include "propagation/propagation.hpp"
#include "report/format.hpp"
#include "report/propagation_report.hpp"

using namespace spodvr;

namespace {

nlohmann::json reference_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return nlohmann::json::parse(file);
}

// Setup up to and including the gas-phase states.
struct GasPhaseSetup {
  std::unique_ptr<Potential> pes;  // built once, as the program does
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  PotentialGrid potential;
  FftwFft fft;
  GasPhase gas;

  explicit GasPhaseSetup(const nlohmann::json& json)
      : pes(make_potential(parse_input(json))),
        config(build_config(parse_input(json), *pes)),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)),
        potential(compute_potential_grid(*pes, config, grids, angular)),
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

namespace {

// The hard-wall input with the analysis surface fixed by the input (DynamZSpAna false), at ZAnaSp if
// given (the Fortran's ZAnalysisSp > 0).
nlohmann::json fixed_surface_json(std::optional<double> z_ana_sp) {
  nlohmann::json json = reference_json();
  json["DynamZSpAna"] = false;
  if (z_ana_sp) json["ZAnaSp"] = *z_ana_sp;
  return json;
}

// The set-up report of a propagation. noipa: with a Style known at compile time, GCC 15's LTO
// specialises print_propagation_setup and then warns falsely (-Wstringop-overflow) about the
// heading string of print_section.
[[gnu::noipa]] std::string propagation_setup_report(const SpatialGrids& grids, const WavePacketPropagation& propagation,
                                                    bool color) {
  std::ostringstream out;
  print_propagation_setup(out, Style{color}, grids, propagation);
  return out.str();
}

}  // namespace

// ZAnaSp fixes the surface one point before the first specular point at or above it (within the
// 1e-12 of .RGE.): 16 bohr is the point 225, 16 - 5e-13 counts as on it, 16.03 lies before the point
// 226. The set-up report prints that point ("fixed at Z").
TEST(AnalysisSurface, FixedAtZAnaSp) {
  const GasPhaseSetup s(fixed_surface_json(16.0));
  const UniformGrid& z = s.grids.z_specular;
  ASSERT_NEAR(z[225], 16.0, 1e-12);
  const InitialWavePacket wp = make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);
  const AnalysisSurface surface(s.config, s.grids, wp.specular);
  EXPECT_TRUE(surface.fixed());
  EXPECT_EQ(surface.index(), 224);
  EXPECT_EQ(surface.fixed_line(), z[224]);
  EXPECT_LT(surface.index(), s.grids.z_specular_optical->first - 1);  // not the absorber's point

  for (const auto& [z_ana_sp, index] : {std::pair{16.0 - 5e-13, 224}, std::pair{16.03, 225}}) {
    const GasPhaseSetup other(fixed_surface_json(z_ana_sp));
    const InitialWavePacket other_wp = make_initial_wave_packet(other.config, other.grids, other.gas, other.fft);
    EXPECT_EQ(AnalysisSurface(other.config, other.grids, other_wp.specular).index(), index) << z_ana_sp;
  }

  // The report, without and with colours.
  const WavePacketPropagation propagation(s.config, s.grids, s.angular, s.potential, s.gas, s.fft);
  for (const bool color : {false, true}) {
    const std::string text = propagation_setup_report(s.grids, propagation, color);
    EXPECT_NE(text.find("fixed at Z = " + fixed(z[224], 6)), std::string::npos) << text;
    EXPECT_EQ(text.find("dynamic, starts at Z"), std::string::npos) << text;
  }
}

// The last point allowed to carry the initial wave packet (CoherentStateFunc iHigh): before the
// specular absorber when the surface is dynamic or fixed at the absorber, before ZAnaSp when given.
TEST(AnalysisSurface, LastPointInside) {
  const int before_absorber = reference_setup().grids.z_specular_optical->first - 1;
  const GasPhaseSetup& dynamic = reference_setup();
  EXPECT_EQ(make_initial_wave_packet(dynamic.config, dynamic.grids, dynamic.gas, dynamic.fft).info.i_high,
            before_absorber);

  const GasPhaseSetup at_absorber(fixed_surface_json(std::nullopt));
  EXPECT_EQ(
      make_initial_wave_packet(at_absorber.config, at_absorber.grids, at_absorber.gas, at_absorber.fft).info.i_high,
      before_absorber);

  for (const auto& [z_ana_sp, last] : {std::pair{16.0, 224}, std::pair{16.0 - 5e-13, 224}, std::pair{16.03, 225}}) {
    const GasPhaseSetup at_z(fixed_surface_json(z_ana_sp));
    const InitialWavePacket wp = make_initial_wave_packet(at_z.config, at_z.grids, at_z.gas, at_z.fft);
    EXPECT_EQ(wp.info.i_high, last) << z_ana_sp;
    EXPECT_LT(wp.info.i_high, before_absorber);
  }
}

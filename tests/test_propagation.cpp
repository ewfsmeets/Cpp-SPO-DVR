// Tests: kinetic, potential and specular operators, the split-operator step and the averages.
#include <gtest/gtest.h>
#include <omp.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <vector>

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
#include "propagation/propagation.hpp"

using namespace spodvr;
using cplx = std::complex<double>;

namespace {

nlohmann::json reference_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return nlohmann::json::parse(file);
}

// Everything up to the operators, for one input.
struct OperatorSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  PotentialGrid potential;
  FftwFft fft;
  GasPhase gas;
  RotationalBasis basis;
  AngularTransform transform;
  PropagationTables tables;
  ScatteringShape shape;

  explicit OperatorSetup(const nlohmann::json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear, config.input.j_in,
                config.input.mj_in),
        potential(compute_potential_grid(*make_potential(config.input.system_id), config, grids, angular)),
        gas(compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver())),
        basis(angular, config.input.j_in, config.input.mj_in),
        transform(angular, basis, fft),
        tables(make_propagation_tables(config, grids, basis, potential.asymptotic)),
        shape{basis.size(), config.input.nx, config.input.ny, config.input.nr, config.input.nz} {}

  [[nodiscard]] ScatteringWave empty_wave() const {
    return ScatteringWave(shape.n_rot, shape.nx * shape.ny, shape.nr * shape.nz);
  }
};

const OperatorSetup& reference_setup() {
  static const OperatorSetup setup(reference_json());
  return setup;
}

double sum_abs_squared(std::span<const cplx> v) {
  double s = 0.0;
  for (const cplx& c : v) s += std::norm(c);
  return s;
}

void fill_pseudo_random(std::span<cplx> v) {
  for (std::size_t i = 0; i < v.size(); ++i) {
    const double x = static_cast<double>(i);
    v[i] = cplx(std::sin(0.37 * x + 0.1), std::cos(1.31 * x));
  }
}

double max_difference(std::span<const cplx> a, std::span<const cplx> b) {
  double d = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) d = std::max(d, std::abs(a[i] - b[i]));
  return d;
}

}  // namespace

TEST(ScatteringKinetic, ConservesTheNormAndHalfStepsComposeToFullSteps) {
  const OperatorSetup& s = reference_setup();
  const ScatteringKineticOperator kinetic(s.shape, s.tables, s.fft);
  ScatteringWave psi = s.empty_wave();
  fill_pseudo_random(psi.values);
  const double norm = sum_abs_squared(psi.values);

  ScatteringWave full = psi;
  kinetic.apply_zr(full, TimeStep::Full);
  kinetic.apply_xy(full, TimeStep::Full);
  EXPECT_NEAR(sum_abs_squared(full.values) / norm, 1.0, 1e-13);

  ScatteringWave halves = psi;
  for (int i = 0; i < 2; ++i) {
    kinetic.apply_zr(halves, TimeStep::Half);
    kinetic.apply_xy(halves, TimeStep::Half);
  }
  EXPECT_LT(max_difference(full.values, halves.values), 1e-12);
}

TEST(ScatteringKinetic, FreeGaussianSpreadsAnalyticallyInZ) {
  // psi(Z, t) = (2a/pi)^(1/4) / sqrt(1 + i b) exp(-a (Z - Z0)^2 / (1 + i b)), b = 2 a t / M,
  // the same for every (irot, x, y, r) (r and xy kinetic act trivially on a constant).
  const OperatorSetup& s = reference_setup();
  const ScatteringKineticOperator kinetic(s.shape, s.tables, s.fft);
  const double a = 1.0;
  const double z0 = 5.0;
  const double mass = s.config.derived.total_mass;
  const int n_steps = 100;
  const double t = n_steps * s.config.input.time_step;
  const auto gaussian = [&](double z, double time) {
    const cplx denominator(1.0, 2.0 * a * time / mass);
    return std::pow(2.0 * a / M_PI, 0.25) / std::sqrt(denominator) * std::exp(-a * (z - z0) * (z - z0) / denominator);
  };

  ScatteringWave psi = s.empty_wave();
  for (int izr = 0; izr < psi.n_blocks(); ++izr)
    for (cplx& c : psi.block(izr)) c = gaussian(s.grids.z_scattering[izr / s.shape.nr], 0.0);
  for (int i = 0; i < n_steps; ++i) kinetic.apply_zr(psi, TimeStep::Full);

  double worst = 0.0;
  for (int izr = 0; izr < psi.n_blocks(); izr += 37)
    for (const cplx& c : psi.block(izr))
      worst = std::max(worst, std::abs(c - gaussian(s.grids.z_scattering[izr / s.shape.nr], t)));
  EXPECT_LT(worst, 1e-12);
}

TEST(ScatteringKinetic, ParallelPlaneWaveGetsItsKineticPhase) {
  // With parallel momentum k0, the plane wave exp(i (k0 + G) . R) is an eigenfunction of the
  // xy step (shift, FFT, exp(-i T_xy dt), inverse FFT, unshift) with the tabulated energy of G.
  nlohmann::json json = reference_json();
  json["EParallel"] = 0.05;
  json["PhiParallel"] = 30.0;
  const OperatorSetup s(json);
  const DerivedParameters& d = s.config.derived;
  const ScatteringKineticOperator kinetic(s.shape, s.tables, s.fft);
  const int n = 1;   // G = n delta_kx along x
  const int m = -1;  // and m delta_ky along y (FFT index ny - 1)
  ScatteringWave psi = s.empty_wave();
  for (int izr = 0; izr < psi.n_blocks(); ++izr)
    for (int ix = 0; ix < s.shape.nx; ++ix)
      for (int iy = 0; iy < s.shape.ny; ++iy) {
        const double phase = (d.kx + n * d.delta_kx) * ix * d.dx + (d.ky + m * d.delta_ky) * iy * d.dy;
        for (int irot = 0; irot < s.shape.n_rot; ++irot)
          psi.block(izr)[irot + s.shape.n_rot * (iy + s.shape.ny * ix)] = std::exp(cplx(0.0, phase));
      }
  const ScatteringWave initial = psi;
  kinetic.apply_xy(psi, TimeStep::Full);
  const cplx expected_factor = s.tables.xy_kinetic.full[(s.shape.ny + m) + s.shape.ny * n];
  double worst = 0.0;
  for (std::size_t i = 0; i < psi.values.size(); i += 7)
    worst = std::max(worst, std::abs(psi.values[i] - expected_factor * initial.values[i]));
  EXPECT_LT(worst, 1e-13);
}

TEST(PotentialOperator, ConstantPotentialIsAPhaseAndTheInitialChannelMovesToTheSpecularGrid) {
  const OperatorSetup& s = reference_setup();
  const double a = 0.3;  // -V dt
  const std::vector<double> potential(s.potential.values.size(), a);
  const int nzpot = s.grids.nz_potential;
  const int ichn = s.basis.initial_index();
  const PotentialOperator op(s.shape, nzpot, ichn, s.tables, s.transform, potential);

  WavePacket psi{s.empty_wave(), SpecularWave(s.shape.nr, s.grids.z_specular.size())};
  fill_pseudo_random(psi.scattering.values);
  const WavePacket initial = psi;
  op.apply(psi);

  // A point inside the interaction region and outside all absorbers.
  const int iz = nzpot / 2;
  const int ir = 10;
  const int izr = ir + s.shape.nr * iz;
  ASSERT_EQ(s.tables.z_optical[iz], 1.0);
  ASSERT_EQ(s.tables.r_optical[ir], 1.0);
  const int n_rot = s.shape.n_rot;
  const int n_xy = s.shape.nx * s.shape.ny;
  const std::span<const cplx> before = initial.scattering.block(izr);
  const std::span<const cplx> after = psi.scattering.block(izr);
  for (int ixy = 0; ixy < n_xy; ++ixy)
    for (int irot = 0; irot < n_rot; ++irot) {
      cplx expected = std::exp(cplx(0.0, a)) * s.tables.rotation.full[irot + n_rot * ir] * before[irot + n_rot * ixy];
      cplx value = after[irot + n_rot * ixy];
      if (irot == ichn) value += s.tables.specular_to_scattering[ixy] * psi.specular(ir, iz);
      EXPECT_NEAR(std::abs(value - expected), 0.0, 1e-14) << ixy << " " << irot;
    }
  // What remains in the initial channel has no component along the initial (x, y) function.
  cplx overlap(0.0, 0.0);
  for (int ixy = 0; ixy < n_xy; ++ixy) overlap += s.tables.scattering_to_specular[ixy] * after[ichn + n_rot * ixy];
  EXPECT_NEAR(std::abs(overlap), 0.0, 1e-14);

  // Beyond the interaction region the gas-phase potential acts and the specular grid is untouched.
  const int izr_out = ir + s.shape.nr * (nzpot + 5);
  EXPECT_EQ(psi.specular(ir, nzpot + 5), cplx(0.0, 0.0));
  const cplx expected_out = s.tables.rotation.full[0 + n_rot * ir] * s.tables.gas_phase[ir];
  EXPECT_NEAR(std::abs(psi.scattering.block(izr_out)[0] - expected_out * initial.scattering.block(izr_out)[0]), 0.0,
              1e-14);
}

TEST(SplitOperator, AnalysisStepKeepsTheNormAndDoesNotDependOnThreads) {
  // One analysis step (5 time steps) of the initial wave packet: no norm reaches an absorber
  // yet, and the result must be bitwise identical for 1 and 4 threads.
  const OperatorSetup& s = reference_setup();
  std::vector<double> potential = s.potential.values;
  premultiply_with_time_step(potential, s.config.input.time_step);
  const int ichn = s.basis.initial_index();
  const ScatteringKineticOperator kinetic(s.shape, s.tables, s.fft);
  const PotentialOperator potential_op(s.shape, s.grids.nz_potential, ichn, s.tables, s.transform, potential);
  const SpecularOperator specular(s.shape.nr, s.grids.z_specular.size(), s.grids.nz_potential, ichn,
                                  *s.grids.z_specular_optical, s.grids.r_optical, s.tables, s.fft);
  const SplitOperatorPropagator propagator(s.config.derived.steps_per_analysis, kinetic, potential_op, specular);
  const InitialWavePacket initial = make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);

  const auto one_step = [&](int threads) {
    omp_set_num_threads(threads);
    WavePacket psi{s.empty_wave(), initial.specular};
    propagator.propagate_analysis_step(psi);
    return psi;
  };
  const WavePacket serial = one_step(1);
  const WavePacket parallel = one_step(4);
  EXPECT_EQ(serial.scattering.values, parallel.scattering.values);
  EXPECT_EQ(serial.specular.values, parallel.specular.values);

  const Averages before = compute_averages(0.0, WavePacket{s.empty_wave(), initial.specular}, s.config, s.grids,
                                           s.grids.z_specular_optical->first - 1);
  const Averages after = compute_averages(10.0, serial, s.config, s.grids, s.grids.z_specular_optical->first - 1);
  // Only norm inside the absorbers can be lost. In this case that is mainly the tail of the
  // v = 0 state at large r (the periodic FGH grid wraps the decay below r_min around to r_max;
  // the Fortran eigenvectors have the same tail), which the absorber in r removes.
  double in_absorber = 0.0;
  for (int iz = 0; iz < s.grids.z_specular.size(); ++iz)
    for (int ir = 0; ir < s.shape.nr; ++ir)
      if (iz >= s.grids.z_specular_optical->first || ir >= s.grids.r_optical.first)
        in_absorber += std::norm(initial.specular(ir, iz));
  in_absorber *= s.config.input.dz * s.config.input.dr;
  EXPECT_GT(in_absorber, 1e-8);
  EXPECT_NEAR(before.norm, 1.0, 1e-6);
  EXPECT_LE(after.norm, before.norm + 1e-13);
  EXPECT_GE(after.norm, before.norm - in_absorber);
  EXPECT_LT(after.z_mean, before.z_mean);  // moving towards the surface
}

TEST(Averages, NormAndMeanPositionsOfASimpleState) {
  // A single specular point and a single scattering point with known densities.
  const OperatorSetup& s = reference_setup();
  WavePacket psi{s.empty_wave(), SpecularWave(s.shape.nr, s.grids.z_specular.size())};
  const double volume_spec = s.config.input.dz * s.config.input.dr;
  const double volume_scat = volume_spec * s.config.derived.xy_volume_element;
  psi.specular(4, 200) = std::sqrt(0.25 / volume_spec);
  psi.scattering.block(7 + s.shape.nr * 30)[3] = std::sqrt(0.75 / volume_scat);
  const int surface = 210;
  const Averages a = compute_averages(1.0, psi, s.config, s.grids, surface);
  EXPECT_NEAR(a.norm, 1.0, 1e-14);
  EXPECT_NEAR(a.z_mean, 0.25 * s.grids.z_specular[200] + 0.75 * s.grids.z_scattering[30], 1e-13);
  EXPECT_NEAR(a.r_mean, 0.25 * s.grids.r[4] + 0.75 * s.grids.r[7], 1e-13);
  EXPECT_NEAR(a.rest_norm, 200 > s.grids.iz_analysis && 200 <= surface ? 0.25 : 0.0, 1e-14);
}

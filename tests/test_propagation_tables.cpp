// Tests: precomputed propagation factors (kinetic, rotational, absorbing, shift, projection).
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <vector>

#include <nlohmann/json.hpp>

#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "propagation/propagation_tables.hpp"

using namespace spodvr;
using cplx = std::complex<double>;

namespace {

struct TablesSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  RotationalBasis basis;
  std::vector<double> vasy;
  PropagationTables tables;

  explicit TablesSetup(const nlohmann::json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear, config.input.j_in,
                config.input.mj_in),
        basis(angular, config.input.j_in, config.input.mj_in),
        vasy(config.input.nr, 0.01),
        tables(make_propagation_tables(config, grids, basis, vasy)) {}
};

nlohmann::json reference_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return nlohmann::json::parse(file);
}

void expect_consistent_steps(const StepFactors& f) {
  ASSERT_EQ(f.full.size(), f.half.size());
  for (std::size_t i = 0; i < f.full.size(); ++i) {
    EXPECT_NEAR(std::abs(f.full[i]), 1.0, 1e-15) << i;
    EXPECT_NEAR(std::abs(f.half[i] * f.half[i] - f.full[i]), 0.0, 1e-15) << i;
  }
}

}  // namespace

TEST(PropagationTables, SizesOfTheReferenceCase) {
  const TablesSetup s(reference_json());
  const PropagationTables& t = s.tables;
  EXPECT_EQ(t.z_kinetic.full.size(), 176u);
  EXPECT_EQ(t.r_kinetic.full.size(), 56u);
  EXPECT_EQ(t.xy_kinetic.full.size(), 16u);
  EXPECT_EQ(t.specular_kinetic.full.size(), 56u * 280u);
  EXPECT_EQ(t.rotation.full.size(), 15u * 56u);
  EXPECT_EQ(t.z_optical.size(), 176u);
  EXPECT_EQ(t.z_specular_optical.size(), 280u);
  EXPECT_EQ(t.gas_phase.size(), 56u);
}

TEST(PropagationTables, KineticFactorsAreUnitaryAndHalfStepsSquareToFullSteps) {
  const TablesSetup s(reference_json());
  for (const StepFactors* f :
       {&s.tables.z_kinetic, &s.tables.r_kinetic, &s.tables.xy_kinetic, &s.tables.specular_kinetic, &s.tables.rotation})
    expect_consistent_steps(*f);
}

TEST(PropagationTables, KineticFactorsUseTheGridEnergies) {
  const TablesSetup s(reference_json());
  const double dt = s.config.input.time_step;
  const auto ez = s.grids.kz_scattering.kinetic_energies();
  for (std::size_t i = 0; i < ez.size(); ++i) EXPECT_EQ(s.tables.z_kinetic.full[i], exp_minus_i(ez[i] * dt));
  // Specular grid (ir, iz): E_r + E_Z (+ parallel energy, zero here).
  const auto er = s.grids.kr.kinetic_energies();
  const auto esp = s.grids.kz_specular.kinetic_energies();
  EXPECT_EQ(s.tables.parallel_energy, 0.0);
  EXPECT_EQ(s.tables.specular_kinetic.full[3 + 56 * 7], exp_minus_i((er[3] + esp[7]) * dt));
}

TEST(PropagationTables, ParallelEnergiesAtNormalIncidenceOnRectangularCell) {
  // kx = ky = 0 and a 90 degree cell: E(ix, iy) = (kx^2 + ky^2) / 2M with FFT-ordered momenta.
  const TablesSetup s(reference_json());
  const DerivedParameters& d = s.config.derived;
  const int n = 4;
  const double k_fft[n] = {0.0, 1.0, 2.0, -1.0};
  for (int ix = 0; ix < n; ++ix)
    for (int iy = 0; iy < n; ++iy) {
      const double kx = k_fft[ix] * d.delta_kx;
      const double ky = k_fft[iy] * d.delta_ky;
      EXPECT_NEAR(s.tables.xy_kinetic_energies[iy + n * ix], (kx * kx + ky * ky) / (2.0 * d.total_mass), 1e-15);
    }
}

TEST(PropagationTables, RotationalEnergies) {
  const TablesSetup s(reference_json());
  const int nrot = s.basis.size();
  for (int ir = 0; ir < 56; ir += 11)
    for (int irot = 0; irot < nrot; ++irot) {
      const int j = s.basis.j(irot);
      const double r = s.grids.r[ir];
      EXPECT_NEAR(s.tables.rotational_energies[irot + nrot * ir],
                  j * (j + 1) / (2.0 * s.config.derived.reduced_mass * r * r), 1e-16);
    }
}

TEST(PropagationTables, AbsorbingPotentialProfile) {
  const TablesSetup s(reference_json());
  const InputParameters& p = s.config.input;
  const DerivedParameters& d = s.config.derived;
  const double dt = p.time_step;
  const double strongest = std::exp(-d.z_opt_factor * 0.5 * (d.optical_potential_power + 1) * dt);
  double previous = 1.0;
  for (int iz = 0; iz < s.grids.z_scattering.size(); ++iz) {
    const double z = s.grids.z_scattering[iz];
    const double f = s.tables.z_optical[iz];
    if (z < p.z_opt_start || z > p.z_opt_end) {
      EXPECT_EQ(f, 1.0) << iz;
    } else {
      EXPECT_LE(f, previous) << iz;  // monotonically stronger absorption towards the end
      EXPECT_GE(f, strongest * (1.0 - 1e-15)) << iz;
      previous = f;
    }
  }
  EXPECT_LT(previous, 1.0);
}

TEST(PropagationTables, ProjectionPhasesWithParallelMomentum) {
  // With a parallel momentum the shift and projection phases are non-trivial. Forward and
  // backward shifts cancel, and projecting onto the initial (x, y) function after adding it
  // gives back exactly one: sum_xy XYProj * initXY = 1.
  nlohmann::json json = reference_json();
  json["EParallel"] = 0.05;
  json["PhiParallel"] = 30.0;
  const TablesSetup s(json);
  const PropagationTables& t = s.tables;
  EXPECT_GT(t.parallel_energy, 0.0);
  cplx projection(0.0, 0.0);
  for (std::size_t ixy = 0; ixy < t.shift_forward.size(); ++ixy) {
    EXPECT_NEAR(std::abs(t.shift_forward[ixy] * t.shift_backward[ixy] - 1.0), 0.0, 1e-15);
    projection += t.scattering_to_specular[ixy] * t.specular_to_scattering[ixy];
  }
  EXPECT_NEAR(std::abs(projection - 1.0), 0.0, 1e-14);
  EXPECT_NE(t.shift_forward[5], cplx(1.0, 0.0));
  expect_consistent_steps(t.xy_kinetic);
}

TEST(PropagationTables, FftScaleFactors) {
  const TablesSetup s(reference_json());
  EXPECT_EQ(s.tables.fft_scale_z, 0.08 / std::sqrt(2.0 * M_PI));
  EXPECT_EQ(s.tables.fft_scale_r, 0.15 / std::sqrt(2.0 * M_PI));
  EXPECT_EQ(s.tables.fft_scale_x, std::sqrt(s.config.derived.dx / 4));
}

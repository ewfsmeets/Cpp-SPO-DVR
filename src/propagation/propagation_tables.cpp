// Precomputed propagation factors: kinetic and rotational exponentials, absorbing potentials,
// momentum shift and grid projection phases.
// Fortran: GridEnergies.f90, Absorption.f90, MomentumShifting.f90, TwoGrids.f90, KineticOp.f90
#include "propagation/propagation_tables.hpp"

#include <cmath>
#include <stdexcept>

#include "core/constants.hpp"
#include "core/fortran_math.hpp"

namespace spodvr {

namespace {

using cplx = std::complex<double>;

// Full step exp(-i E dt); half step exp(-i E dt / 2) (Fortran: exp(arg/2.0), used for xy and
// the rotation) or the principal square root of the full step (Fortran: sqrt(full), used for
// Z, r and the specular grid).
enum class HalfStep { HalfPhase, SquareRoot };

StepFactors step_factors(std::span<const double> energies, double dt, HalfStep half_step) {
  StepFactors f;
  for (const double e : energies) {
    const double phase = e * dt;
    f.full.push_back(exp_minus_i(phase));
    f.half.push_back(half_step == HalfStep::HalfPhase ? exp_minus_i(phase / 2.0) : std::sqrt(f.full.back()));
  }
  return f;
}

// Momenta along one reciprocal lattice direction in FFT order, shifted by the initial momentum
// k0: k0, k0 + dk, ..., k0 + (n/2) dk, ..., k0 - dk. The offset is accumulated as in Fortran.
// Fortran: GridEnergies.f90 DifChanEnergies (kinx, kiny)
std::vector<double> parallel_momenta(int n, double k0, double dk) {
  std::vector<double> k(n);
  k[0] = k0;
  k[n / 2] = k0 + (n / 2) * dk;
  double offset = dk;
  for (int i = 1; i < n / 2; ++i) {
    k[i] = k0 + offset;
    k[n - i] = k0 - offset;
    offset = offset + dk;
  }
  return k;
}

// Kinetic energy of the parallel motion on the skewed lattice for every (kx, ky) pair,
// ekinxy = (kx^2 + ky^2 - 2 cos(skew) kx ky) / 2M. Fortran: DifChanEnergies (nmoder = 3, 4)
std::vector<double> parallel_kinetic_energies(const Config& config) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  const std::vector<double> kx = parallel_momenta(p.nx, d.kx, d.delta_kx);
  const std::vector<double> ky = parallel_momenta(p.ny, d.ky, d.delta_ky);
  const double ecnvrt = 0.5 / d.total_mass;
  const double prefac = -2.0 * std::cos(d.system.skewing_angle);

  std::vector<double> energies;
  for (int ix = 0; ix < p.nx; ++ix) {
    for (int iy = 0; iy < p.ny; ++iy) {
      const double ekinx = ecnvrt * (kx[ix] * kx[ix]);
      const double ekiny = ecnvrt * (ky[iy] * ky[iy]);
      const double ekxky = ecnvrt * prefac * kx[ix] * ky[iy];
      energies.push_back(ekinx + ekxky + ekiny);
    }
  }
  return energies;
}

// exp(-V_opt dt) for the quadratic (power nabpow) absorbing potential
// V_opt = vabs (n+1)/2 ((q - start) / (end - start))^n on [start, end], 0 elsewhere.
// Fortran: Absorption.f90 SetupAbsorpPot (vimag, vimar, vimas and EXPOpt*)
std::vector<double> optical_factors(const UniformGrid& grid, double start, double end, double strength, int power,
                                    double dt) {
  const double length = end - start;
  std::vector<double> f;
  for (int i = 0; i < grid.size(); ++i) {
    const double q = grid[i];
    double v_opt = 0.0;
    if (q >= start && q <= end)
      v_opt = strength * 0.5 * (power + 1) * fortran_integer_power((q - start) / length, power);
    f.push_back(std::exp(-(v_opt * dt)));
  }
  return f;
}

// Phase k0 . R of the initial parallel momentum at every (x, y) grid point (skewed lattice).
// Fortran: MomentumShifting.f90 / TwoGrids.f90 (conx, conxy)
std::vector<double> parallel_phases(const Config& config) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  const double sin_skew = std::sin(d.system.skewing_angle);
  std::vector<double> phases;
  for (int ix = 0; ix < p.nx; ++ix) {
    const double conx = sin_skew * d.kx * ix * d.dx;
    for (int iy = 0; iy < p.ny; ++iy) phases.push_back(conx + sin_skew * d.ky * iy * d.dy);
  }
  return phases;
}

}  // namespace

PropagationTables make_propagation_tables(const Config& config, const SpatialGrids& grids, const RotationalBasis& basis,
                                          std::span<const double> asymptotic_potential) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  const double dt = p.time_step;
  if (p.nr < 2) throw std::invalid_argument("make_propagation_tables: nR = 1 is not supported");
  PropagationTables t;

  // Kinetic energies. Fortran: SetupGridEnergies
  t.z_kinetic = step_factors(grids.kz_scattering.kinetic_energies(), dt, HalfStep::SquareRoot);
  t.r_kinetic = step_factors(grids.kr.kinetic_energies(), dt, HalfStep::SquareRoot);
  // T_Z and T_r applied together after a 2D FFT: exp(-i T_Z dt) exp(-i T_r dt) per (kr, kZ).
  for (std::size_t iz = 0; iz < t.z_kinetic.full.size(); ++iz)
    for (std::size_t ir = 0; ir < t.r_kinetic.full.size(); ++ir) {
      t.zr_kinetic.full.push_back(t.z_kinetic.full[iz] * t.r_kinetic.full[ir]);
      t.zr_kinetic.half.push_back(t.z_kinetic.half[iz] * t.r_kinetic.half[ir]);
    }

  t.xy_kinetic_energies = parallel_kinetic_energies(config);
  t.parallel_energy = t.xy_kinetic_energies[0];
  t.xy_kinetic = step_factors(t.xy_kinetic_energies, dt, HalfStep::HalfPhase);

  // Specular grid: Z and r kinetic energy plus the (constant) parallel energy. Fortran: ensp + edifin
  std::vector<double> specular_energies;
  for (const double ez : grids.kz_specular.kinetic_energies()) {
    const double ez_total = ez + t.parallel_energy;
    for (const double er : grids.kr.kinetic_energies()) specular_energies.push_back(er + ez_total);
  }
  t.specular_kinetic = step_factors(specular_energies, dt, HalfStep::SquareRoot);

  // Rotational energy J(J+1) / (2 mu r^2) on the r grid. Fortran: RotChanEnergies
  for (int ir = 0; ir < grids.r.size(); ++ir) {
    const double r = grids.r[ir];
    for (int irot = 0; irot < basis.size(); ++irot) {
      const int j = basis.j(irot);
      t.rotational_energies.push_back((j * (j + 1)) / (2.0 * d.reduced_mass * (r * r)));
    }
  }
  t.rotation = step_factors(t.rotational_energies, dt, HalfStep::HalfPhase);

  // Absorbing potentials. Fortran: SetupAbsorpPot
  const int power = d.optical_potential_power;
  t.z_optical = optical_factors(grids.z_scattering, p.z_opt_start, p.z_opt_end, d.z_opt_factor, power, dt);
  t.r_optical = optical_factors(grids.r, p.r_opt_start, p.r_opt_end, d.r_opt_factor, power, dt);
  t.z_specular_optical = optical_factors(grids.z_specular, p.zsp_opt_start, p.zsp_opt_end, d.zsp_opt_factor, power, dt);

  for (const double v : asymptotic_potential) t.gas_phase.push_back(exp_minus_i(v * dt));

  // Momentum shift and grid projection. Fortran: SetupMomentumShift, SetupTwoGrids
  for (const double phase : parallel_phases(config)) {
    const cplx forward = exp_minus_i(phase);
    const cplx backward = std::exp(cplx(0.0, phase));
    t.shift_forward.push_back(forward);
    t.shift_backward.push_back(backward);
    t.specular_to_scattering.push_back(backward * d.unit_cell_norm);
    t.scattering_to_specular.push_back(forward * d.unit_cell_norm * d.xy_volume_element);
  }

  // FFT scale factors. Fortran: SetupScaleFactors
  const double sqrt_two_pi = std::sqrt(2.0 * Constants::pi);
  t.fft_scale_z = p.dz / sqrt_two_pi;
  t.fft_scale_r = p.dr / sqrt_two_pi;
  t.fft_scale_x = std::sqrt(d.dx / p.nx);
  t.fft_scale_y = std::sqrt(d.dy / p.ny);
  return t;
}

}  // namespace spodvr

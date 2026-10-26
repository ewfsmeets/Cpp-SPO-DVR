#pragma once

#include <complex>
#include <span>
#include <vector>

#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"

namespace spodvr {

// exp(-i E dt) for a full and a half time step.
struct StepFactors {
  std::vector<std::complex<double>> full;
  std::vector<std::complex<double>> half;
};

// Everything the split-operator propagator multiplies with, precomputed once (plain arrays).
// Layouts: xy index ixy = iy + ny*ix, specular (ir, iz) = ir + nr*iz, rotation (irot, ir).
// Fortran: GridEnergies.f90 SetupGridEnergies (EXPKinEn), Absorption.f90 SetupAbsorpPot,
// MomentumShifting.f90 SetupMomentumShift, TwoGrids.f90 SetupTwoGrids, GasPhase.f90 XVGasPhase,
// KineticOp.f90 SetupScaleFactors
struct PropagationTables {
  // --- kinetic and rotational energy factors (Fortran: EXPKinEn%fullStep / %halfStep) -------
  StepFactors z_kinetic;          // nz. Fortran: ZScat
  StepFactors r_kinetic;          // nr. Fortran: RScat
  StepFactors zr_kinetic;         // (ir, iz): products of the Z and r factors (2D FFT kinetic step)
  StepFactors xy_kinetic;         // nx*ny. Fortran: XY
  StepFactors specular_kinetic;   // nr*nzsp, includes the parallel energy. Fortran: Spec
  StepFactors rotation;           // nrot*nr. Fortran: RotChan

  std::vector<double> xy_kinetic_energies;  // ekinxy(ix, iy) at ixy (hartree). Fortran: ekinxy
  std::vector<double> rotational_energies;  // J(J+1) / (2 mu r^2) at (irot, ir). Fortran: erot
  double parallel_energy = 0.0;             // kinetic energy parallel to the surface. Fortran: edifin

  // --- absorption: exp(-V_opt dt), real (Fortran: EXPOptZ, EXPOptR, EXPOptZSpec) -------------
  std::vector<double> z_optical;            // nz
  std::vector<double> r_optical;            // nr
  std::vector<double> z_specular_optical;   // nzsp

  // --- gas-phase potential exp(-i vasy(r) dt). Fortran: XVGasPhase ---------------------------
  std::vector<std::complex<double>> gas_phase;  // nr

  // --- momentum shift of the parallel motion (nx*ny). Fortran: shiftForward, shiftBackward ---
  std::vector<std::complex<double>> shift_forward;
  std::vector<std::complex<double>> shift_backward;

  // --- specular <-> scattering grid projection (nx*ny). Fortran: initXY, XYProj -------------
  std::vector<std::complex<double>> specular_to_scattering;
  std::vector<std::complex<double>> scattering_to_specular;

  // --- FFT scale factors: d / sqrt(2 pi) in Z, r; sqrt(d / n) in x, y. Fortran: scalFac* ------
  double fft_scale_z = 0.0;
  double fft_scale_r = 0.0;
  double fft_scale_x = 0.0;
  double fft_scale_y = 0.0;
};

// asymptotic_potential: vasy(r) in hartree (PotentialGrid::asymptotic).
[[nodiscard]] PropagationTables make_propagation_tables(const Config& config, const SpatialGrids& grids,
                                                        const RotationalBasis& basis,
                                                        std::span<const double> asymptotic_potential);

// exp(-i phase), evaluated like the Fortran exp(-ima * phase).
[[nodiscard]] inline std::complex<double> exp_minus_i(double phase) {
  return std::exp(std::complex<double>(0.0, -phase));
}

}  // namespace spodvr

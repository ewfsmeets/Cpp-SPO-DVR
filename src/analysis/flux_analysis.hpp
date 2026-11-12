#pragma once

#include <complex>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "core/aligned_vector.hpp"
#include "core/array_index.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"

// The flux analysis (Fortran FluxAnalyzor.f90). Its source file is compiled with strict floating
// point in every build, also the fast-math one (CMake target spodvr_strict_objects): the results
// must not depend on the slot or thread count, and c_n's complex division must be gfortran's.

namespace spodvr {

// Number of flux analysis energies, INT((EMaxFluxAna - EMinFluxAna) / DeltaEFluxAna) + 1 with the
// energies in hartree: the end point is dropped for most inputs (flux Q3). Shared by
// make_flux_energies and the memory estimate, so the two cannot disagree. 0 if doFluxAnalysis is off.
// Throws InputError (DeltaEFluxAna) above max_flux_energies.
// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor (numScatEnergies)
[[nodiscard]] int flux_energy_count(const InputParameters& input);

// Upper limit of the number of flux energies: far above any count whose coefficients fit in memory
// (each energy holds two complex values per point of the flux surface), and safe for int.
inline constexpr int max_flux_energies = 1'000'000;

// The energies of the flux analysis and the per-energy factors of its time-to-energy transform,
// in atomic units. The Fortran recomputes k_n, b(k_n) and A_n at every analysis step from the same
// operands; here they are computed once (the same values).
// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor, storeTimeIndependentWF, computeFluxProbability
struct FluxEnergies {
  std::vector<double> z_energy;  // E_n = E_min + n dE: the initial state's energy in Z. Fortran: ZScatEnergy
  std::vector<double> total;     // (E_n + parallel energy) + rovibrational energy. Fortran: totalScatEnergy
  std::vector<double> momentum;  // k_n = sqrt(2 M E_n). Fortran: kZ_initial
  std::vector<std::complex<double>> amplitude;         // b(k_n) of the initial wave packet. Fortran: bkz
  std::vector<std::complex<double>> transform_factor;  // A_n = k_n / ((2 pi M) b(k_n)). Fortran: prefactor
  std::vector<std::string> warnings;                   // EMaxFluxAna not reached (flux Q3)

  [[nodiscard]] int size() const { return static_cast<int>(z_energy.size()); }
};

// The flux energies of the input: EMinFluxAna, EMaxFluxAna, DeltaEFluxAna and the wave packet's
// EMinWP, EMaxWP converted with units::precise (flux Q5), M = total mass, b(k) centred at Z0.
// parallel_energy: PropagationTables::parallel_energy (Fortran edifin); initial_state_energy:
// GasPhase::initial.energy (Fortran eviroi); both hartree.
// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor, storeTimeIndependentWF
[[nodiscard]] FluxEnergies make_flux_energies(const Config& config, double parallel_energy,
                                              double initial_state_energy);

// A = k / ((2 pi M) b): the factor of the time-to-energy transform at momentum k (1/bohr), for the
// total mass M (m_e) and the wave packet's momentum amplitude b = b(k). The one complex division of
// the flux analysis: Smith's algorithm, as gfortran's (-fcx-fortran-rules), never conj(b) / |b|^2.
// Fortran: FluxAnalyzor.f90 storeTimeIndependentWF (prefactor)
[[nodiscard]] std::complex<double> flux_transform_factor(double momentum, double mass, std::complex<double> amplitude);

// c_n(t) = (A_n exp(i E_tot,n t)) dt: the weight of the analysis step that ends at `time` in the
// coefficients of energy n; dt (time_step) is the analysis time step. All atomic units.
// Fortran: FluxAnalyzor.f90 storeTimeIndependentWF (coeff)
[[nodiscard]] std::complex<double> flux_step_coefficient(const FluxEnergies& energies, int n, double time,
                                                         double time_step);

// Weights D(ir), ir = 0 .. nr-1, of the Fourier-grid derivative d psi / dr at the flux surface
// r_c = r(ir_flux), with ir_flux 0-based (Fortran irFluxAnalysis - 1) and dr in bohr:
// D(ir) = ((2c) x) / dr with c = pi / nr, d = ir_flux - ir and x = (-1)^d / (2 tan(c d)); D(ir_flux) = 0.
// Fortran: WaveFunction.f90 getWaveFuncDerivAtZr (ddr(:, irFluxAnalysis))
[[nodiscard]] std::vector<double> flux_derivative_weights(int nr, double dr, int ir_flux);

// d psi / dr at the flux surface from one r column (column[ir] = psi(r_ir)) and the weights above:
// the sum in the Fortran order, from 0 over every ir ascending, the zero term D(ir_flux) included.
// The reference of the backend tests (CPU and CUDA), which compare with it by EXPECT_EQ: the backends
// loop ir-outer over the strided wave function, with the same order per element.
// Fortran: WaveFunction.f90 getWaveFuncDerivAtZr
[[nodiscard]] std::complex<double> flux_r_derivative(std::span<const double> weights,
                                                     std::span<const std::complex<double>> column);

// The wave function psi and its r-derivative at the flux surface r = r_c of one analysis step, on
// the planes iz < nz_potential; layout (block, iz), block = (irot, ixy) as in the scattering wave.
// One slot of the FluxAccumulator, which the backend fills (every element: a reused slot still holds
// an earlier step's values). Zero on construction.
// Fortran: WaveFunction.f90 getWaveFunctionAtZr, getWaveFuncDerivAtZr (one iz each)
struct FluxSurfaceWave {
  Index2D index;             // (block, iz)
  ComplexVector psi;         // psi(irot, ixy, ir_c, iz)
  ComplexVector derivative;  // d psi / dr at (irot, ixy, ir_c, iz) (1/bohr)

  FluxSurfaceWave(int n_blocks, int n_planes);

  [[nodiscard]] int n_blocks() const { return static_cast<int>(index.n0); }
  [[nodiscard]] int n_planes() const { return static_cast<int>(index.n1); }
};

// The time-to-energy transforms phi(E_n) of psi and phi'(E_n) of d psi / dr at the flux surface,
// summed over the analysis steps: phi = sum_t c_n(t) psi(t). Layout (block, iz, n), as the Fortran
// arrays and its restart files. Zero on construction.
// Fortran: FluxAnalyzor.f90 phi, phiDeriv
struct FluxCoefficients {
  Index3D index;                 // (block, iz, n)
  ComplexVector phi;             // Fortran: phi
  ComplexVector phi_derivative;  // Fortran: phiDeriv (1/bohr)

  FluxCoefficients(int n_blocks, int n_planes, int n_energies);

  [[nodiscard]] int n_blocks() const { return static_cast<int>(index.n0); }
  [[nodiscard]] int n_planes() const { return static_cast<int>(index.n1); }
  [[nodiscard]] int n_energies() const { return static_cast<int>(index.n2); }
};

// The grid spacings and masses of the flux probability.
// Fortran: FluxAnalyzor.f90 computeFluxProbability (nrotw, nx, ny, delz, delx, dely, skewingAngle,
// rmass, rdmass)
struct FluxGrid {
  int n_rot;             // rotational functions per (x, y) point. Fortran: nrotw
  int nx;
  int ny;
  double dz;             // bohr. Fortran: delz
  double dx;             // bohr. Fortran: delx
  double dy;             // bohr. Fortran: dely
  double skewing_angle;  // rad. Fortran: skewingAngle
  double total_mass;     // M (m_e). Fortran: rmass
  double reduced_mass;   // mu (m_e). Fortran: rdmass
};

// The FluxGrid of a run with n_rot rotational functions (RotationalBasis::size()).
[[nodiscard]] FluxGrid make_flux_grid(const Config& config, const SpatialGrids& grids, int n_rot);

// The reaction probability from the flux through r = r_c, per flux energy.
// Fortran: FluxAnalyzor.f90 computeFluxProbability (reacProbability, reacProbDensity)
struct FluxResults {
  std::vector<double> energy;       // E_n (eV, converted with Constants_eV2Ha: flux Q5)
  std::vector<double> probability;  // P(E_n). Fortran: reacProbability
  // (ix, iy, n) with ix fastest, as the Fortran's reacProbDensity(nx, ny, nE): the transpose of the
  // (x, y) order of psi's blocks, ixy = iy + ny ix.
  Index3D density_index;
  std::vector<double> density;  // P(x, y; E_n) (1/bohr^2). Fortran: reacProbDensity
};

// P(x, y; E_n) = prefac_n sum_iz sum_irot Im(conj(phi) phi'), prefac_n = ((2 pi M) / (k_n mu)) dz,
// and P(E_n) = ((sum_iy sum_ix P(x, y; E_n)) dx) dy sin(skewing angle), in the Fortran's order.
// Fortran: FluxAnalyzor.f90 computeFluxProbability
[[nodiscard]] FluxResults flux_results(const FluxCoefficients& coefficients, const FluxEnergies& energies,
                                       const FluxGrid& grid);

// The time loop's side of the flux analysis: psi and d psi / dr at the flux surface of each
// analysis step go into one of slot_count slots, which are added to the coefficients when all are
// full and on flush(). The caller flushes before it reads the coefficients.
// Fortran: FluxAnalyzor.f90 storeFluxCoeffs, storeTimeIndependentWF
class FluxAccumulator {
 public:
  // n_blocks = n_rot nx ny, n_planes = nz_potential; time_step: the analysis time step (a.u.),
  // the dt of c_n. Throws std::invalid_argument for a slot count below 1.
  FluxAccumulator(int n_blocks, int n_planes, FluxEnergies energies, double time_step, int slot_count);

  // Adds the analysis step that ends at `time` (a.u.): `fill` writes psi and d psi / dr at the flux
  // surface into a free slot, every element (the slot may hold an earlier step's values). Flushes
  // when the slots are full.
  void add_step(double time, const std::function<void(FluxSurfaceWave&)>& fill);

  // Adds the filled slots to the coefficients.
  void flush();

  [[nodiscard]] const FluxEnergies& energies() const { return energies_; }
  [[nodiscard]] int slot_count() const { return static_cast<int>(slots_.size()); }
  [[nodiscard]] int pending_steps() const { return pending_; }  // filled slots not yet added

  // phi and phi' of every step added so far. Throws std::logic_error while steps are pending: the
  // coefficients are complete only after flush().
  [[nodiscard]] const FluxCoefficients& coefficients() const;

 private:
  void add_tile(int iz, int first_block, int tile_size, int n);

  FluxEnergies energies_;
  double time_step_;  // Fortran: dtAnalysis
  FluxCoefficients coefficients_;
  std::vector<FluxSurfaceWave> slots_;
  Index2D slot_coefficient_index_;                       // (n, slot)
  std::vector<std::complex<double>> slot_coefficients_;  // c_n(t) of the step in each slot
  int pending_ = 0;
};

}  // namespace spodvr

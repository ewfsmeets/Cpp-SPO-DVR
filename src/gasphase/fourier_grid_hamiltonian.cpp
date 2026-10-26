// Fourier grid Hamiltonian (FGH) for the H2 vibration. Fortran: FourierGridHamiltonian.f90
#include "gasphase/fourier_grid_hamiltonian.hpp"

#include <cmath>
#include <complex>
#include <stdexcept>

#include "core/constants.hpp"
#include "fft/fft_sizes.hpp"

namespace spodvr {

// Fortran: FourierGridHamiltonian.f90 InitFourierGridHamiltonian, SetUpGrid, SetUpKineticMatrix
FourierGridHamiltonian::FourierGridHamiltonian(int n_approx, double lower_bound, double upper_bound,
                                               double reduced_mass, const Fft& fft)
    : n_(nearest_fft_size(n_approx, 1)), dx_((upper_bound - lower_bound) / n_) {
  const int n = n_;
  grid_.resize(n);
  for (int l = 0; l < n; ++l) grid_[l] = l * dx_ + lower_bound;

  // Kinetic energy in momentum space for the half spectrum, k_l = 2 pi l / (n dx).
  const int n_half = n / 2 + 1;
  std::vector<double> kinetic_diagonal(n_half);
  for (int l = 0; l < n_half; ++l)
    kinetic_diagonal[l] = 2.0 / reduced_mass * std::pow(l * Constants::pi / n / dx_, 2.0);

  // T = c2r( diag(E_k) * r2c(identity) ), column by column. The forward transform has the
  // Fortran "forward" sign (+) and scale 1; the inverse has scale 1/n.
  std::vector<double> identity(n * n, 0.0);
  for (int l = 0; l < n; ++l) identity[l + n * l] = 1.0;
  std::vector<std::complex<double>> transformed(n_half * n);

  const FftLayout real_layout{n, n, 1, n};
  const FftLayout half_layout{n, n, 1, n_half};
  fft.plan_real_to_complex(real_layout, half_layout, FftSign::Plus, 1.0)->execute(identity, transformed);
  for (int column = 0; column < n; ++column)
    for (int l = 0; l < n_half; ++l) transformed[l + n_half * column] *= kinetic_diagonal[l];

  kinetic_.assign(n * n, 0.0);
  const double inverse_scale = 1.0 / (1.0 * n);
  fft.plan_complex_to_real(half_layout, real_layout, FftSign::Minus, inverse_scale)->execute(transformed, kinetic_);
}

namespace {

// Makes the first local maximum of |psi|, seen from the last grid point, positive.
// Fortran: GetEigenvaluesAndEigenvectors (after the sort)
void fix_sign(std::span<double> psi) {
  std::size_t i = psi.size() - 1;
  while (i > 0 && std::abs(psi[i - 1]) >= std::abs(psi[i])) --i;
  if (psi[i] < 0.0)
    for (double& y : psi) y = -1.0 * y;
}

// Normalises psi to sum(psi^2) * dx = 1. Fortran: NormaliseEigenvectors
void normalise(std::span<double> psi, double dx) {
  double norm = 0.0;
  for (const double y : psi) norm = norm + y * y;
  norm = norm * dx;
  norm = 1.0 / std::sqrt(norm);
  for (double& y : psi) y = y * norm;
}

}  // namespace

EigenDecomposition FourierGridHamiltonian::eigenstates(std::span<const double> potential,
                                                       const SymmetricEigensolver& solver) const {
  const int n = n_;
  if (potential.size() != static_cast<std::size_t>(n))
    throw std::invalid_argument("FourierGridHamiltonian: potential size mismatch");

  // H = T + diag(V). Fortran: BuildHamiltonian
  std::vector<double> hamiltonian = kinetic_;
  for (int l = 0; l < n; ++l) hamiltonian[l + n * l] = hamiltonian[l + n * l] + potential[l];

  EigenDecomposition result = solver.solve(hamiltonian, n);
  for (int k = 0; k < n; ++k) {
    const auto psi = std::span<double>(result.vectors).subspan(k * n, n);
    fix_sign(psi);
    normalise(psi, dx_);
  }
  return result;
}

}  // namespace spodvr

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
    : n_(nearest_fft_size(n_approx, 1)), dx_((upper_bound - lower_bound) / static_cast<double>(n_)) {
  const auto n = static_cast<std::size_t>(n_);
  grid_.resize(n);
  for (std::size_t l = 0; l < n; ++l) grid_[l] = static_cast<double>(l) * dx_ + lower_bound;

  // Kinetic energy in momentum space for the half spectrum, k_l = 2 pi l / (n dx).
  const std::size_t n_half = n / 2 + 1;
  std::vector<double> kinetic_diagonal(n_half);
  for (std::size_t l = 0; l < n_half; ++l)
    kinetic_diagonal[l] = 2.0 / reduced_mass *
                          std::pow(static_cast<double>(l) * Constants::pi / static_cast<double>(n_) / dx_, 2.0);

  // T = c2r( diag(E_k) * r2c(identity) ), column by column. The forward transform has the
  // Fortran "forward" sign (+) and scale 1; the inverse has scale 1/n.
  std::vector<double> identity(n * n, 0.0);
  for (std::size_t l = 0; l < n; ++l) identity[l + n * l] = 1.0;
  std::vector<std::complex<double>> transformed(n_half * n);

  const FftLayout real_layout{n_, n_, 1, n_};
  const FftLayout half_layout{n_, n_, 1, static_cast<int>(n_half)};
  fft.plan_real_to_complex(real_layout, half_layout, FftSign::Plus, 1.0)->execute(identity, transformed);
  for (std::size_t m = 0; m < n; ++m)
    for (std::size_t l = 0; l < n_half; ++l) transformed[l + n_half * m] *= kinetic_diagonal[l];

  kinetic_.assign(n * n, 0.0);
  const double inverse_scale = 1.0 / (1.0 * static_cast<double>(n_));
  fft.plan_complex_to_real(half_layout, real_layout, FftSign::Minus, inverse_scale)
      ->execute(transformed, kinetic_);
}

EigenDecomposition FourierGridHamiltonian::eigenstates(std::span<const double> potential,
                                                       const SymmetricEigensolver& solver) const {
  const auto n = static_cast<std::size_t>(n_);
  if (potential.size() != n) throw std::invalid_argument("FourierGridHamiltonian: potential size mismatch");

  // Fortran: BuildHamiltonian
  std::vector<double> hamiltonian = kinetic_;
  for (std::size_t l = 0; l < n; ++l) hamiltonian[l + n * l] = hamiltonian[l + n * l] + potential[l];

  EigenDecomposition result = solver.solve(hamiltonian, n_);

  for (std::size_t k = 0; k < n; ++k) {
    const auto column = std::span<double>(result.vectors).subspan(k * n, n);

    // Sign convention. Fortran: GetEigenvaluesAndEigenvectors (after the sort)
    std::size_t j = n;  // 1-based as in the Fortran loop
    while (j > 1 && !(std::abs(column[j - 2]) < std::abs(column[j - 1]))) --j;
    if (column[j - 1] < 0.0)
      for (double& v : column) v = -1.0 * v;

    // Fortran: NormaliseEigenvectors
    double norm = 0.0;
    for (const double v : column) norm = norm + v * v;
    norm = norm * dx_;
    norm = 1.0 / std::sqrt(norm);
    for (double& v : column) v = v * norm;
  }
  return result;
}

}  // namespace spodvr

#pragma once

#include <span>
#include <vector>

#include "fft/fft.hpp"
#include "linalg/symmetric_eigensolver.hpp"

namespace spodvr {

// Fourier grid Hamiltonian for a 1D coordinate (the H2 bond length r) on its own grid
// x_l = l * dx + lower, dx = (upper - lower) / n, with n = nearest_fft_size(n_approx).
// Fortran: FourierGridHamiltonian.f90
class FourierGridHamiltonian {
 public:
  FourierGridHamiltonian(int n_approx, double lower_bound, double upper_bound, double reduced_mass, const Fft& fft);

  [[nodiscard]] int size() const { return n_; }
  [[nodiscard]] double spacing() const { return dx_; }
  [[nodiscard]] std::span<const double> grid() const { return grid_; }
  // Kinetic energy matrix, column-major n x n. Fortran: FGHType%Kinetic (after SetUpKineticMatrix)
  [[nodiscard]] std::span<const double> kinetic() const { return kinetic_; }

  // Eigenpairs of T + diag(potential), sorted by energy. Each eigenvector is normalised with
  // sum(psi^2) * dx = 1 and its sign is fixed such that the first local maximum of |psi|, seen
  // from the last grid point, is positive.
  // Fortran: UsePotentialInHamiltonian + GetEigenvaluesAndEigenvectors + NormaliseEigenvectors
  [[nodiscard]] EigenDecomposition eigenstates(std::span<const double> potential,
                                               const SymmetricEigensolver& solver) const;

 private:
  int n_;
  double dx_;
  std::vector<double> grid_;
  std::vector<double> kinetic_;
};

}  // namespace spodvr

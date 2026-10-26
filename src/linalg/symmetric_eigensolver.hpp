#pragma once

#include <memory>
#include <span>
#include <vector>

#include "core/numerical_libraries.hpp"

namespace spodvr {

// Eigenpairs of a real symmetric n x n matrix, eigenvalues in ascending order.
// Eigenvector k is column k of `vectors` (column-major: component i at vectors[i + n*k]).
struct EigenDecomposition {
  int n = 0;
  std::vector<double> values;
  std::vector<double> vectors;

  [[nodiscard]] std::span<const double> vector(int k) const {
    return std::span<const double>(vectors).subspan(static_cast<std::size_t>(k) * static_cast<std::size_t>(n),
                                                    static_cast<std::size_t>(n));
  }
};

// Abstract diagonaliser so that other backends (e.g. the Jacobi method of the Fortran code,
// jacobi.f90) can be added without touching the callers.
class SymmetricEigensolver {
 public:
  virtual ~SymmetricEigensolver() = default;
  // `matrix` is column-major n x n; only its lower triangle is guaranteed to be referenced.
  [[nodiscard]] virtual EigenDecomposition solve(std::span<const double> matrix, int n) const = 0;
};

// LAPACK dsyevd (divide and conquer).
class LapackEigensolver final : public SymmetricEigensolver {
 public:
  LapackEigensolver() { configure_numerical_libraries(); }
  [[nodiscard]] EigenDecomposition solve(std::span<const double> matrix, int n) const override;
};

}  // namespace spodvr

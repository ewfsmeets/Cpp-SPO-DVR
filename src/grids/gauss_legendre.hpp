#pragma once

#include <vector>

namespace spodvr {

struct QuadratureRule {
  std::vector<double> roots;
  std::vector<double> weights;
};

// Roots and weights of the associated-Legendre Gauss quadrature of the given order and |mj|,
// scaled to the interval [a, b] (the Fortran code calls it with a = 1, b = -1, which yields
// descending roots and negative weights). Newton iteration converges to
// fortran_machine_epsilon(). Fortran: legendre.f90 GaussAssociatedLegendreQuadr
[[nodiscard]] QuadratureRule gauss_associated_legendre(double a, double b, int order, int mj);

}  // namespace spodvr

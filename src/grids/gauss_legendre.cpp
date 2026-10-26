// Gauss quadrature for associated Legendre functions. Fortran: legendre.f90
#include "grids/gauss_legendre.hpp"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "core/constants.hpp"
#include "core/machine_epsilon.hpp"

namespace spodvr {

namespace {

// log(n!) accumulated as in the Fortran code. Fortran: legendre.f90 SetupLegendreModule
std::vector<double> log_factorials(int n_max) {
  std::vector<double> lf(n_max + 1);
  lf[0] = 0.0;
  for (int n = 1; n <= n_max; ++n)
    lf[n] = std::log(n) + lf[n - 1];
  return lf;
}

}  // namespace

// Fortran: legendre.f90 GaussAssociatedLegendreQuadr. Integer/real mixing of the Fortran
// expressions is reproduced explicitly so that the floating-point operations are identical.
QuadratureRule gauss_associated_legendre(double a, double b, int order, int mj) {
  if (order <= 0) throw std::invalid_argument("gauss_associated_legendre: order must be > 0");
  const int amj = std::abs(mj);
  const double eps = fortran_machine_epsilon();
  const auto lf = log_factorials(2 * order + 2 * amj + 5);
  auto lfac = [&](int n) { return lf[n]; };
  const double ord = order;

  QuadratureRule rule{std::vector<double>(order), std::vector<double>(order)};
  auto root = [&](int i) -> double& { return rule.roots[i - 1]; };
  auto weight = [&](int i) -> double& { return rule.weights[i - 1]; };

  const int m = (order + 1) / 2;
  const double sb = 0.5 * (b + a);
  const double sa = 0.5 * (b - a);
  const double xmj_order = amj / ord;

  double fw = 0.0;
  if (amj != 0) {
    fw = std::pow(2.0, 2 * amj + 1) * (order + amj) * sa *
         std::exp(2.0 * lfac(order + amj - 1) - lfac(order) - lfac(order + 2 * amj));
  } else {
    fw = 2.0 * sa / ord;
  }

  double z = 0.0;
  for (int i = 1; i <= m; ++i) {
    const double di = i;
    if (i == (m - order / 2) * m) {
      z = 0.0;
    } else if (amj == 0) {
      z = std::cos(Constants::pi * (di - 0.25) / (ord + 0.5));
    } else if (i == 1) {
      const double damj = amj;
      z = 1.0 - ((1.0 + damj) * (2.78 / (4.0 + ord * ord) + 0.768 * xmj_order / ord)) /
                    (1.0 + xmj_order * (2.44 + 1.282 * xmj_order));
    } else if (i == 2) {
      const double damj = amj;
      z = z - (1.0 - z) * (1.0 + 0.012 * damj * (1.0 + 0.25 * damj) / ord) *
                  (1.0 + 0.06 * (order - 8) * (1.0 + 0.12 * damj) / ord) *
                  (4.1 + damj) / ((1.0 + damj) * (1.0 + 0.156 * damj));
    } else {
      z = (root(i - 2) - 2.0 * root(i - 1) + sb) / sa;
    }

    double p1 = 0.0;
    double p2 = 0.0;
    double dp = 0.0;
    for (int iteration = 1;; ++iteration) {
      // Order-th polynomial at z by a Jacobi-type recursion.
      p1 = 1.0;
      p2 = 0.0;
      for (int j = 1; j <= order; ++j) {
        const double p3 = p2;
        p2 = p1;
        p1 = (j + amj) *
             ((2 * (j + amj) - 1) * z * p2 - (j + amj - 1) * p3) /
             (j * (j + 2 * amj));
      }
      dp = (ord * z * p1 - (order + amj) * p2) / (z * z - 1.0);
      const double z_previous = z;
      z = z_previous - p1 / dp;
      if (std::abs(z - z_previous) <= eps) break;
      if (iteration >= 500) break;
    }

    root(i) = sb - sa * z;
    weight(i) = fw / (dp * p2);
    root(order + 1 - i) = sb + sa * z;
    weight(order + 1 - i) = weight(i);
  }
  return rule;
}

}  // namespace spodvr

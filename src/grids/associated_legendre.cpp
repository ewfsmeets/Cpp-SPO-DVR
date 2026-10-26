// Normalised associated Legendre functions by upward recursion in J. Fortran: legendre.f90
#include "grids/associated_legendre.hpp"

#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "core/fortran_math.hpp"

namespace spodvr {

std::vector<double> log_factorials(int n_max) {
  std::vector<double> lf(n_max + 1);
  lf[0] = 0.0;
  for (int n = 1; n <= n_max; ++n) lf[n] = std::log(n) + lf[n - 1];
  return lf;
}

// Fortran: legendre.f90 AssociatedLegendre (expression order kept for bitwise results).
std::vector<double> associated_legendre(int j_max, int m, double theta) {
  const int am = std::abs(m);
  if (j_max < 0 || am > j_max) throw std::invalid_argument("associated_legendre: need 0 <= |m| <= j_max");
  const std::vector<double> lf = log_factorials(j_max + am + 1);
  std::vector<double> p(j_max - am + 1);

  // Start of the recursion: P_|m|^|m| up to its J-dependent normalisation.
  double p1 = 1.0;
  if (am > 0) p1 = std::exp(lf[2 * am] / 2 - lf[am]) * fortran_integer_power(std::sin(theta) / 2, am);
  if (m > 0 && am % 2 == 1) p1 = -p1;  // Condon-Shortley phase
  p[0] = std::sqrt(am + 0.5) * p1;

  // Upward recursion (J - |m| + 1) P_{J+1} = (2J + 1) x P_J - (J + |m|) P_{J-1}.
  const double x = std::cos(theta);
  double p2 = 0.0;
  for (int l = am; l < j_max; ++l) {
    const double p3 = p2;
    p2 = p1;
    p1 = ((2 * l + 1) * x * p2 - (l + am) * p3) / (l - am + 1);
    p[l - am + 1] = std::sqrt(l + 1.5) * p1;
  }
  if (am == 0) return p;

  // Remaining normalisation for m != 0: sqrt((J - |m|)! (2|m|)! / (J + |m|)!).
  for (int i = 1; i < static_cast<int>(p.size()); ++i)
    p[i] = std::exp((lf[i] + lf[2 * am] - lf[2 * am + i]) / 2) * p[i];
  return p;
}

}  // namespace spodvr

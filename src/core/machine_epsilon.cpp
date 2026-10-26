// Run-time machine epsilon as computed by the Fortran code. Fortran: constants.f90 Calculate_Constant_EPS
#include "core/machine_epsilon.hpp"

#include <cmath>

namespace spodvr {

namespace {

double compute_machine_epsilon() {
  // volatile forces every intermediate to be rounded to double, as in the Fortran loop.
  volatile double eps = 1.0;
  while (true) {
    volatile double a = 1.0 - eps;
    if (a >= 1.0) break;
    eps = eps / 2.0;
  }
  double result = std::abs(eps);
  if (result <= 0.0) result = 1.0e-12;  // Fortran fallback.
  return result;
}

}  // namespace

double fortran_machine_epsilon() {
  static const double eps = compute_machine_epsilon();
  return eps;
}

}  // namespace spodvr

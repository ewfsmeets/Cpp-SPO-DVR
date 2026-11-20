#pragma once
// Reductions for the comparisons of the tests (CPU, GPU and CRP library tests): a maximum that
// keeps NaN, so that a check on the largest deviation or on a scale fails when a value is NaN.

#include <algorithm>
#include <cmath>
#include <limits>

namespace spodvr::testing {

// max(a, b), or NaN if a or b is NaN. std::max drops NaN (std::max(x, NaN) == x), but a gate must
// fail on NaN: EXPECT_LE(NaN, tolerance) fails. For finite a, b the result is std::max(a, b).
inline double nan_aware_max(double a, double b) {
  if (std::isnan(a) || std::isnan(b)) return std::numeric_limits<double>::quiet_NaN();
  return std::max(a, b);
}

// Half a unit in the last printed digit of a value the Fortran prints as E15.8 (a mantissa
// 0.dddddddd): a C++ value agrees to all printed digits when |c - f| <= this. The 1% margin covers
// the decimal conversion of the printed value.
inline double printed_half_unit(double fortran_value) {
  if (fortran_value == 0.0) return 0.5e-8 * 1e-99;
  const double exponent = std::floor(std::log10(std::abs(fortran_value))) + 1.0;
  return 1.01 * 0.5 * std::pow(10.0, exponent - 8.0);
}

}  // namespace spodvr::testing

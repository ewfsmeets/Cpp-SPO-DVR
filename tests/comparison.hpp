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

}  // namespace spodvr::testing

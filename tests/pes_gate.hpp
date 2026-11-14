#pragma once
// Gate 1b (docs/validation.md): the C++ 6D potential against the subsample of the Fortran .pes
// file in pes_reference.txt (tools/fortran_reference/extract_pes.py). The comparison takes the C++
// values as a list in the order of the reference points, so that a test can show that the gate
// fails on wrong values. Used by test_regression_pes.cpp; negative tests in test_comparison_helpers.cpp.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "comparison.hpp"
#include "input/config.hpp"
#include "potential/potential_grid.hpp"
#include "setup_reference.hpp"

namespace spodvr::testing {

// |a − b| relative to the larger of |a|, |b| and floor: near V = 0 (the LEPS cancellation) a pure
// relative difference is meaningless, so the potential gates use a floor of 1 eV. NaN if a or b is.
inline double scaled_difference(double a, double b, double floor) {
  return std::abs(a - b) / std::max({std::abs(a), std::abs(b), floor});
}

// The VMin / VMax clip of the potential, in hartree; empty when the input does not impose it.
struct PotentialLimits {
  std::optional<double> v_min;  // ImposeVMin: VMin
  std::optional<double> v_max;  // ImposeVMax: VMax
};

// The limits a run imposes on its potential.
inline PotentialLimits potential_limits(const Config& config) {
  PotentialLimits limits;
  if (config.input.impose_v_min) limits.v_min = config.derived.v_min;
  if (config.input.impose_v_max) limits.v_max = config.derived.v_max;
  return limits;
}

// A reference value clipped to the limits. The Fortran writes its .pes file before it imposes them
// (Potential.f90 setupThePES); its potential in memory, like the C++ grid, is the reference clipped.
inline double with_limits(double v, const PotentialLimits& limits) {
  if (limits.v_min && v < *limits.v_min) v = *limits.v_min;
  if (limits.v_max && v > *limits.v_max) v = *limits.v_max;
  return v;
}

// The C++ potential (hartree, before the −dt premultiplication) at the reference points, in their order.
inline std::vector<double> potential_at_points(const PotentialGrid& potential, const test::PesReference& ref) {
  std::vector<double> values;
  values.reserve(ref.points.size());
  for (const test::PesReference::Point& p : ref.points) {
    const std::size_t itp = p.iphi + ref.nphi * p.itheta;
    const std::size_t ixy = p.iy + ref.ny * p.ix;
    const std::size_t izr = p.ir + ref.nr * p.iz;
    values.push_back(potential(itp, ixy, izr));
  }
  return values;
}

// The outcome of gate 1b on a list of C++ values.
struct PesGateResult {
  int clipped = 0;                   // reference values that the limits change (NaN not counted)
  double max_difference = 0.0;       // the largest scaled difference; NaN if a value is NaN
  std::size_t worst = 0;             // the point of max_difference (the first one; with NaN, the first NaN)
  std::vector<std::size_t> failing;  // the points above the tolerance, NaN included
};

// Gate 1b on values (one per reference point, hartree): the scaled difference of each from the
// clipped reference value, with floor (hartree), against tolerance.
inline PesGateResult compare_pes_values(const test::PesReference& ref, const std::vector<double>& values,
                                        const PotentialLimits& limits, double floor, double tolerance) {
  PesGateResult result;
  for (std::size_t i = 0; i < ref.points.size(); ++i) {
    const double expected = with_limits(ref.points[i].v, limits);
    // A NaN reference value is not clipped (NaN != NaN); it fails below.
    if (!std::isnan(ref.points[i].v) && expected != ref.points[i].v) ++result.clipped;
    const double difference = scaled_difference(values.at(i), expected, floor);
    // A new maximum, or the first NaN (true for a NaN difference while the maximum is still finite).
    if (!std::isnan(result.max_difference) && !(difference <= result.max_difference)) result.worst = i;
    result.max_difference = nan_aware_max(result.max_difference, difference);
    if (!(difference <= tolerance)) result.failing.push_back(i);
  }
  return result;
}

// Gate 1b as a test: one failure per failing point (its 0-based grid indices, the C++ and the
// clipped reference value), then the summary line of the regression tests.
inline void expect_pes_gate(const std::string& name, const test::PesReference& ref, const std::vector<double>& values,
                            const PotentialLimits& limits, double floor, double tolerance) {
  ASSERT_EQ(values.size(), ref.points.size());
  const PesGateResult result = compare_pes_values(ref, values, limits, floor, tolerance);
  for (const std::size_t i : result.failing) {
    const test::PesReference::Point& p = ref.points[i];
    const double expected = with_limits(p.v, limits);
    ADD_FAILURE() << "itheta " << p.itheta << " iphi " << p.iphi << " ix " << p.ix << " iy " << p.iy << " iz " << p.iz
                  << " ir " << p.ir << ": " << values[i] << " vs " << expected << ", scaled difference "
                  << scaled_difference(values[i], expected, floor) << " above " << tolerance;
  }
  std::cout << "  " << name << ": " << ref.points.size() << " points (" << result.clipped
            << " clipped to VMin/VMax), max scaled difference " << result.max_difference << '\n';
}

}  // namespace spodvr::testing

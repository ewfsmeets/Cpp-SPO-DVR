#pragma once

#include <cctype>
#include <cmath>
#include <string_view>

namespace spodvr {

// Tolerant real comparisons with an *absolute* tolerance of 1e-12, reproducing the Fortran
// user-defined operators so that accept/reject decisions in input validation match exactly.
// Fortran: MyOperators.f90 (.REQ. .RNE. .RGE. .RGT. .RLE. .RLT.)
namespace fortran_compare {

inline constexpr double tolerance = 1.0e-12;

[[nodiscard]] inline bool req(double a, double b) { return std::abs(a - b) < tolerance; }
[[nodiscard]] inline bool rne(double a, double b) { return !req(a, b); }
[[nodiscard]] inline bool rgt(double a, double b) { return (a - b) > tolerance; }
[[nodiscard]] inline bool rge(double a, double b) {
  const double d = a - b;
  return std::abs(d) < tolerance || d > tolerance;
}
[[nodiscard]] inline bool rle(double a, double b) { return !rgt(a, b); }
[[nodiscard]] inline bool rlt(double a, double b) { return !rge(a, b); }

// Case-insensitive comparison ignoring leading and trailing blanks.
// Fortran: MyOperators.f90 .SEQ. (compareStringsForEQ)
[[nodiscard]] inline bool seq(std::string_view a, std::string_view b) {
  auto trim = [](std::string_view s) {
    const auto first = s.find_first_not_of(' ');
    if (first == std::string_view::npos) return std::string_view{};
    const auto last = s.find_last_not_of(' ');
    return s.substr(first, last - first + 1);
  };
  a = trim(a);
  b = trim(b);
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (std::tolower(static_cast<unsigned char>(a[i])) !=
        std::tolower(static_cast<unsigned char>(b[i])))
      return false;
  }
  return true;
}

}  // namespace fortran_compare
}  // namespace spodvr

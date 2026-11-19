#pragma once

#include <cmath>
#include <string>
#include <string_view>

#include "core/error.hpp"

namespace spodvr::validation {

// An input check: throws InputError naming the parameter when `violated`.
inline void check(bool violated, std::string_view parameter, const std::string& message) {
  if (violated) throw InputError(std::string(parameter), message);
}

// Fortran NINT: round half away from zero.
inline double nint(double x) { return static_cast<double>(std::lround(x)); }

}  // namespace spodvr::validation

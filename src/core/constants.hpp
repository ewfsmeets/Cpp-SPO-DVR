#pragma once

#include <cmath>
#include <complex>

namespace spodvr {

// Mathematical and program constants, values copied verbatim from the Fortran code.
// Fortran: constants.f90 (Constants_PI parameter; the derived values are computed in
// SetupConstantsModule with the same expressions as below).
struct Constants {
  static constexpr double pi = 3.1415926535897932384626433832795028841971;
  static constexpr double two_pi = 2.0 * pi;
  static constexpr double four_pi = 4.0 * pi;
  static constexpr double inv_pi = 1.0 / pi;
  static constexpr double deg_to_rad = pi / 180.0;  // Fortran: Constants_Deg2Rad
  static inline const double sqrt_2 = std::sqrt(2.0);
  static inline const double sqrt_pi = std::sqrt(pi);
  static inline const double sqrt_2pi = std::sqrt(two_pi);
  static inline const double sqrt_4pi = std::sqrt(four_pi);
  static constexpr std::complex<double> imaginary_unit{0.0, 1.0};  // Fortran: Constants_I
  static constexpr int nen95_sentinel = -95;                         // Fortran: Constants_NEN95
};

}  // namespace spodvr

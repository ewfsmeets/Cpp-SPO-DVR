#pragma once

#include <vector>

namespace spodvr {

// log(n!) for n = 0 .. n_max, accumulated term by term as in the Fortran code.
// Fortran: legendre.f90 SetupLegendreModule (LogOfFactorials)
[[nodiscard]] std::vector<double> log_factorials(int n_max);

// Normalised associated Legendre functions P_J^m(cos theta) for J = |m| .. j_max; element
// J - |m| holds P_J^m. Normalisation: integral over x = cos theta in [-1, 1] of P^2 is 1.
// For m > 0 the Condon-Shortley phase (-1)^m is included; P_J^{-m} = P_J^{|m|} without it.
// Fortran: legendre.f90 AssociatedLegendre
[[nodiscard]] std::vector<double> associated_legendre(int j_max, int m, double theta);

}  // namespace spodvr

#pragma once

namespace spodvr {

// The machine precision exactly as the Fortran code determines it at run time: halve eps
// until 1 - eps rounds to 1 (2^-54 for IEEE double with round-to-nearest-even).
// Used as convergence threshold in the Gauss-Legendre root finder and as Constants_EPS
// elsewhere. Fortran: constants.f90 Calculate_Constant_EPS
double fortran_machine_epsilon();

}  // namespace spodvr

#pragma once
// Numerical constants of the CRP library. Fortran: CRP_Constants.f90 and the literals of the
// CRP_*.f90 modules, copied verbatim.

namespace crplib::detail {

// Fortran: CRP_Constants_PI (the same double as spodvr::Constants::pi; kept local so that the
// library depends on nothing else).
inline constexpr double pi = 3.1415926535897932384626433832795028841971;

// sqrt(3) as the Fortran literal (FOUR6II, FOUR3II, ORIGEN2, MURII). Fortran: sq3
inline constexpr double sq3 = 1.732050807568877293500;
// Fortran: usq3 = 1.0 / sq3 (folded at compile time, the same rounding as at run time)
inline constexpr double usq3 = 1.0 / sq3;

// gfortran folds COS/SIN of constant arguments at compile time. The values the reference objects
// contain (CRP_QD_Wrapper.o: the skewing transform; CRP_C3v_AtomSurface.o: cossurf/sinsurf in
// MURII), as hexadecimal literals; tests check that they equal std::cos/std::sin.
inline constexpr double cos_pi_3 = 0x1.0000000000001p-1;   // cos(pi/3), bits 0x3FE0000000000001
inline constexpr double sin_pi_3 = 0x1.bb67ae8584caap-1;   // sin(pi/3), bits 0x3FEBB67AE8584CAA
inline constexpr double cos_pi_2 = 0x1.1a62633145c07p-54;  // cos(pi/2), bits 0x3C91A62633145C07
inline constexpr double sin_pi_2 = 1.0;                    // sin(pi/2)

}  // namespace crplib::detail

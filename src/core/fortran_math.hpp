#pragma once

namespace spodvr {

// x**n for an integer exponent, evaluated as gfortran does it (libgcc __powidf2: binary
// exponentiation by repeated squaring), so results are bitwise identical to the Fortran code.
// std::pow(double, int) would call pow(double, double) and can differ in the last bit.
inline double fortran_integer_power(double x, int n) {
  unsigned int m = n < 0 ? -static_cast<unsigned int>(n) : static_cast<unsigned int>(n);
  double y = (m % 2 != 0) ? x : 1.0;
  while ((m >>= 1) != 0) {
    x = x * x;
    if (m % 2 != 0) y = y * x;
  }
  return n < 0 ? 1.0 / y : y;
}

}  // namespace spodvr

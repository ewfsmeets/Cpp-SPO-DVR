// Flux check 1 (docs/validation_flux_replay/README.md): compares the C++ b(k) and D bit for bit with
// the hex output of flux_f1_check.f90 (stdin). Usage: flux_f1_check < <Fortran output>
// Exit status: 0 every value bitwise equal, 1 otherwise. Built and run by run_checks.sh.
#include <bit>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>

#include "analysis/flux_analysis.hpp"
#include "core/unit_conversion.hpp"
#include "input/wave_packet_params.hpp"

using namespace spodvr;

int main() {
  const double mass = units::legacy.amu_to_me(1.0078 + 1.0078);
  const double low = units::precise.ev_to_hartree(0.2);
  const double high = units::precise.ev_to_hartree(0.6);
  const std::vector<double> d = flux_derivative_weights(56, 0.15, 8);
  int b_total = 0, b_equal = 0, d_total = 0, d_equal = 0;
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string tag;
    in >> tag;
    if (tag == "B") {
      std::string hk, hr, hi;
      in >> hk >> hr >> hi;
      const double k = std::bit_cast<double>(std::stoull(hk, nullptr, 16));
      const auto b = wave_packet_momentum_amplitude(k, low, high, mass, 11.44);
      ++b_total;
      if (std::bit_cast<std::uint64_t>(b.real()) == std::stoull(hr, nullptr, 16) &&
          std::bit_cast<std::uint64_t>(b.imag()) == std::stoull(hi, nullptr, 16))
        ++b_equal;
      else
        std::printf("b differs at k = %.17g\n", k);
    } else if (tag == "D") {
      int ir;
      std::string hd;
      in >> ir >> hd;
      ++d_total;
      if (std::bit_cast<std::uint64_t>(d[ir - 1]) == std::stoull(hd, nullptr, 16))
        ++d_equal;
      else
        std::printf("D differs at ir = %d\n", ir);
    }
  }
  std::printf("b(k): %d of %d bitwise equal; D: %d of %d bitwise equal\n", b_equal, b_total, d_equal, d_total);
  return (b_equal == b_total && d_equal == d_total && b_total > 0 && d_total > 0) ? 0 : 1;
}

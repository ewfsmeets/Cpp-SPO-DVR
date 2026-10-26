#pragma once

#include <span>
#include <vector>

namespace spodvr {

// Momentum grid conjugate to an n-point coordinate grid with spacing d, in FFT order:
// k_0 = 0, k_i = i*dk and k_{n-i} = -k_i for i = 1 .. n/2-1, k_{n/2} = +(n/2)*dk, with
// dk = 2*pi/(n*d). Kinetic energies E_i = 0.5/mass * (i*dk)^2 for the same ordering.
// Fortran: Setup.f90 (zk/ekinz, rkg/ekinr, zksp/ensp/ezsp)
class MomentumGrid {
 public:
  MomentumGrid(int n, double coordinate_spacing, double mass);

  [[nodiscard]] int size() const { return static_cast<int>(k_.size()); }
  [[nodiscard]] double spacing() const { return dk_; }                  // Fortran: delkz, delkr, delksp
  [[nodiscard]] std::span<const double> k() const { return k_; }        // Fortran: zk, rkg, zksp
  [[nodiscard]] std::span<const double> kinetic_energies() const { return ekin_; }  // ekinz, ekinr, ensp
  [[nodiscard]] double max_kinetic_energy() const { return ekin_[k_.size() / 2]; }  // ekzmax, ekrmax

 private:
  double dk_;
  std::vector<double> k_;
  std::vector<double> ekin_;
};

}  // namespace spodvr

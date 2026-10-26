#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace spodvr {

// Equidistant grid x_i = start + i * spacing, i = 0 .. size-1. The Fortran `start + dfloat(i-1)*delta`
// is evaluated in quad precision, so points can differ by 1 ulp (accepted, see CLAUDE.md). Fortran: Setup.f90 zgrid,
// rgrid, xgrid, ygrid
class UniformGrid {
 public:
  UniformGrid(double start, double spacing, int size) : start_(start), spacing_(spacing) {
    points_.resize(size);
    for (int i = 0; i < size; ++i) points_[i] = start + i * spacing;
  }

  [[nodiscard]] double start() const { return start_; }
  [[nodiscard]] double spacing() const { return spacing_; }
  [[nodiscard]] int size() const { return static_cast<int>(points_.size()); }
  [[nodiscard]] double operator[](int i) const { return points_[i]; }
  [[nodiscard]] double back() const { return points_.back(); }
  [[nodiscard]] std::span<const double> points() const { return points_; }

 private:
  double start_;
  double spacing_;
  std::vector<double> points_;
};

}  // namespace spodvr

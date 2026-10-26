#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace spodvr {

// Equidistant grid x_i = start + i * spacing, i = 0 .. size-1 (computed exactly as the Fortran
// `start + dfloat(i-1)*delta`). Fortran: Setup.f90 zgrid, rgrid, xgrid, ygrid
class UniformGrid {
 public:
  UniformGrid(double start, double spacing, int size) : start_(start), spacing_(spacing) {
    points_.resize(static_cast<std::size_t>(size));
    for (std::size_t i = 0; i < points_.size(); ++i)
      points_[i] = start + static_cast<double>(i) * spacing;
  }

  [[nodiscard]] double start() const { return start_; }
  [[nodiscard]] double spacing() const { return spacing_; }
  [[nodiscard]] int size() const { return static_cast<int>(points_.size()); }
  [[nodiscard]] double operator[](int i) const { return points_[static_cast<std::size_t>(i)]; }
  [[nodiscard]] double back() const { return points_.back(); }
  [[nodiscard]] std::span<const double> points() const { return points_; }

 private:
  double start_;
  double spacing_;
  std::vector<double> points_;
};

}  // namespace spodvr

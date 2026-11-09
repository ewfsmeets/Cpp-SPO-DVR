// The CRP PES for the quantum-dynamics code. Fortran: CRP_QD_Wrapper.f90
#include "crplib/qd_wrapper.hpp"

#include <cmath>

#include "crplib/constants.hpp"
#include "crplib/wrapper.hpp"

namespace crplib {

namespace {

// A small shift for making sure the clamped values stay on the grid. Fortran: Delta = 1.0d-12
// (a quad-precision literal rounded to the double 1e-12; checked in the reference object)
constexpr double delta = 1.0e-12;

// The QD codes work in bohr and eV. Fortran: OutputUnits
constexpr Units output_units{LengthUnit::bohr, EnergyUnit::ev};

}  // namespace

// Fortran: Read_PES_From_Disk: SkewingAngle = pi/3 (C3v, C6v) or pi/2 (C4v) and the transformation
// matrix [[1, cos], [0, sin]] (cos/sin folded by gfortran: crplib/constants.hpp)
QdWrapper::QdWrapper(const DataSet& data, EvaluationOptions options)
    : surface_(data, options), value_only_(options.value_only) {
  const bool square = data.symmetry == Symmetry::c4v;
  skewing_angle_ = square ? detail::pi / 2.0 : detail::pi / 3.0;
  cos_skew_ = square ? detail::cos_pi_2 : detail::cos_pi_3;
  sin_skew_ = square ? detail::sin_pi_2 : detail::sin_pi_3;
}

QdParameters QdWrapper::parameters() const {
  return {surface_.lattice_constant(), surface_.z_asymptotic(), skewing_angle_, symmetry_name(surface_.symmetry())};
}

double QdWrapper::clamped_r(double r) const {
  if (r < surface_.r_min() + delta) return surface_.r_min() + delta;
  if (r > surface_.r_max() - delta) return surface_.r_max() - delta;
  return r;
}

double QdWrapper::clamped_z(double z) const {
  if (z < surface_.z_min() + delta) return surface_.z_min() + delta;
  return z;
}

// Fortran: Get_CRP_QD_Wrapper_Potential_6D
double QdWrapper::potential_ev(double u, double v, double z, double r, double theta, double phi) const {
  const double my_r = clamped_r(r);
  const double my_z = clamped_z(z);
  // Skewed (u, v) to Cartesian: MATMUL(TransMatrix, (u, v))
  const double my_x = 1.0 * u + cos_skew_ * v;
  const double my_y = 0.0 * u + sin_skew_ * v;
  // The atom positions (Fortran: rab)
  const Vec3 a{my_x - 0.5 * my_r * std::sin(theta) * std::cos(phi), my_y - 0.5 * my_r * std::sin(theta) * std::sin(phi),
               my_z - 0.5 * my_r * std::cos(theta)};
  const Vec3 b{my_x + 0.5 * my_r * std::sin(theta) * std::cos(phi), my_y + 0.5 * my_r * std::sin(theta) * std::sin(phi),
               my_z + 0.5 * my_r * std::cos(theta)};
  // The library's V in eV (0 where it gives up silently); S3: without the derivatives
  const double value =
      value_only_ ? convert_energy(surface_.value(a, b).value_or(0.00), EnergyUnit::hartree, output_units.energy)
                  : calculate_molecule_surface(surface_, a, b, output_units).v;
  return value - convert_energy(surface_.h2_asymptote(), EnergyUnit::hartree, EnergyUnit::ev);
}

RangeReport QdWrapper::range_report(std::span<const double> z, std::span<const double> r) const {
  RangeReport report;
  report.z_points = static_cast<int>(z.size());
  report.r_points = static_cast<int>(r.size());
  for (double rr : r)
    if (clamped_r(rr) != rr) ++report.r_clamped;
  for (double zz : z) {
    if (clamped_z(zz) != zz) ++report.z_clamped;
    for (double rr : r) {
      ++report.zr_points;
      if (!surface_.covers(clamped_r(rr), clamped_z(zz))) ++report.zr_uncovered;
    }
  }
  return report;
}

}  // namespace crplib

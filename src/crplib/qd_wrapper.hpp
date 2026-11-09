#pragma once
// The CRP PES as the quantum-dynamics codes use it: V in eV at the molecule's centre (skewed
// lattice coordinates u, v), Z, r and orientation, with r and Z clamped into the range of the data,
// shifted so that the free molecule at its last tabulated r has V = Zdep at infinity. V only (as the
// Fortran). Fortran: CRP_QD_Wrapper.f90 (Get_CRP_QD_Wrapper_Potential_6D, Get_CRP_QD_Wrapper_Params).

#include <span>
#include <string>

#include "crplib/dataset.hpp"
#include "crplib/molecule_surface.hpp"

namespace crplib {

// The fixed parameters of the PES. Fortran: Get_CRP_QD_Wrapper_Params
struct QdParameters {
  double lattice_constant;  // bohr, both lattice vectors
  double z_asymptotic;      // bohr. Fortran: Zasy
  double skewing_angle;     // rad, angle between the lattice vectors (pi/3 or pi/2)
  std::string point_group;  // "C3v", "C4v" or "C6v"
};

// How a set of grid points relates to the range of the data (for the setup report).
struct RangeReport {
  int r_points = 0;
  int r_clamped = 0;  // r outside [Rmin, Rmax]: evaluated at the nearest end
  int z_points = 0;
  int z_clamped = 0;  // Z below Zmin: evaluated at Zmin
  long long zr_points = 0;
  long long zr_uncovered = 0;  // (Z, r) where the Fortran returns V = 0 silently (judged at the grid
                               // point's clamped Z and r; the centre's Z varies by rounding only)
};

class QdWrapper {
 public:
  // Fortran: Read_PES_From_Disk (one data set, no multi-PES header). The options select
  // bitwise-neutral speed-ups (only the tests switch them off).
  explicit QdWrapper(const DataSet& data, EvaluationOptions options = {});

  // V in eV. u, v: the centre in skewed lattice coordinates (bohr along the lattice vectors);
  // theta, phi in radians. Fortran: Get_CRP_QD_Wrapper_Potential_6D
  [[nodiscard]] double potential_ev(double u, double v, double z, double r, double theta, double phi) const;

  [[nodiscard]] QdParameters parameters() const;
  [[nodiscard]] RangeReport range_report(std::span<const double> z, std::span<const double> r) const;
  [[nodiscard]] const MoleculeSurface& surface() const { return surface_; }

  // The bond length and centre height actually evaluated. Fortran: MyR, MyZ
  [[nodiscard]] double clamped_r(double r) const;
  [[nodiscard]] double clamped_z(double z) const;

 private:
  MoleculeSurface surface_;
  double skewing_angle_;
  double cos_skew_, sin_skew_;  // as gfortran folded them
  bool value_only_;  // S3
};

}  // namespace crplib

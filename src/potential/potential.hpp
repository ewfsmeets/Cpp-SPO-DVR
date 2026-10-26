#pragma once

#include <string>

namespace spodvr {

// Fixed parameters of a molecule-surface system, supplied by the PES implementation.
// Fortran: H2SurfaceSystem.f90 performAction(GET_PARAMETERS) -> alat, blat, zasymp,
// skewingAngle, pointGroup.
struct SystemParameters {
  double lattice_constant_x;  // bohr. Fortran: alat
  double lattice_constant_y;  // bohr. Fortran: blat
  double z_asymptotic;        // bohr; Z >= z_asymptotic is the gas phase. Fortran: zasymp
  double skewing_angle;       // rad, angle between the surface lattice vectors. Fortran: skewingAngle
  std::string point_group;    // e.g. "C2v". Fortran: pointGroup
};

// Potential energy surface for a diatomic molecule interacting with a surface.
// Coordinates in bohr / rad, energies in eV.
class Potential {
 public:
  virtual ~Potential() = default;

  // V(x, y, Z, r, theta, phi) in eV. iz and ir are the 0-based grid indices of Z and r; only
  // grid-based PES implementations use them (Fortran: ZIndex, RIndex).
  [[nodiscard]] virtual double value(double x, double y, double z, double r, double theta, double phi, int iz,
                                     int ir) const = 0;

  [[nodiscard]] virtual const SystemParameters& system_parameters() const = 0;

  // Name used to select this PES (Fortran: systemID).
  [[nodiscard]] virtual std::string name() const = 0;
};

}  // namespace spodvr

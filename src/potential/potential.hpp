#pragma once

#include <functional>
#include <optional>
#include <string>

namespace spodvr {

class CrpPotential;
class ScmCrpPotential;

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
  // Conventions every caller relies on:
  // - x, y: skewed coordinates (bohr) along a1 = alat (1, 0) and a2 = blat (cos skew, sin skew): the
  //   Cartesian position is x (1, 0) + y (cos skew, sin skew) (symmetry/surface_lattice.hpp), so x, y
  //   are Cartesian only at a skew of 90 degrees.
  // - Z: height of the molecule's centre (bohr, the PES's own zero); r: bond length (bohr).
  // - theta: polar angle of the bond from the surface normal (+Z); phi: its azimuth from a1 towards a2.
  // - The PES places the atoms (all ported PESs at -/+ r/2): Z is the centre of mass only for equal
  //   masses.
  // - The gas-phase V(r) is V(0, 0, Z, r, pi/2, 0) at zasymp (with a switch the blend also uses the
  //   unswitched zasymp), called with iz = nz_potential (potential_grid.cpp gas_phase_potential).
  // - It is also called off the grid (the symmetry self-check, write2DPES) and from several OpenMP
  //   threads at once.
  // - iz and ir are hints for grid-based PESs only (no ported PES uses them). The 6D grid passes its own
  //   indices (iz < nz_potential); the gas phase iz = nz_potential, one past the potential's Z grid; the
  //   symmetry self-check the Z and r grid indices, with x, y and phi off the grid; the 2D cut
  //   (write2DPES) scattering-grid indices, so iz up to nZ - 1, beyond nz_potential; the SCM edge
  //   diagnostic 0, 0.
  // - For a heteronuclear PES: which atom sits at +r/2 along (theta, phi) is the PES's choice (Persson
  //   and the CRP adapter: the second; Dai-Light: the first; for H2 it does not matter).
  [[nodiscard]] virtual double value(double x, double y, double z, double r, double theta, double phi, int iz,
                                     int ir) const = 0;

  [[nodiscard]] virtual const SystemParameters& system_parameters() const = 0;

  // Name used to select this PES (Fortran: systemID).
  [[nodiscard]] virtual std::string name() const = 0;

  // The CRP data set this PES evaluates (CRPlibrary; for CRPDistrtSurface its CRP part), for the
  // set-up report; empty for every other PES. For CRPDistrtSurface this is not the PES the run
  // evaluates: it has neither the correction nor the rescaled lattice.
  [[nodiscard]] virtual std::optional<std::reference_wrapper<const CrpPotential>> crp_potential() const {
    return std::nullopt;
  }
  // The SCM distorted-surface PES (CRPDistrtSurface), for its report and its identity in restart files;
  // empty for every other PES.
  [[nodiscard]] virtual std::optional<std::reference_wrapper<const ScmCrpPotential>> scm_potential() const {
    return std::nullopt;
  }
};

}  // namespace spodvr

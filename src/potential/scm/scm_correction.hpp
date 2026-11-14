#pragma once
// The SCM correction ΔV: the change of the H₂–surface potential when the surface atoms move from
// their ideal lattice positions to those of a distorted slab snapshot, as a sum of H–metal pair
// functions over the atoms of the top two layers. Fortran: DistortedSurfacePES.dir/SCM_Potential.f90
// (SCM_SurfInteract_3Body, Switch_Rydberg_Potential), SCM_Wrapper.f90 (GetSCMPotential,
// ConvertCoMToXYZ), SCM_Types.f90 (Setup_3Body_Rydberg_Par), DistortedSurfaceWrapper.f90
// (GetCRPDistSurfParams, GetCRPDistSurfPotential).
//
// ΔV at the centre of mass (X, Y) (Cartesian), Z, bond length r and orientation (θ, φ), bitwise as
// the Fortran computes it (the comments in scm_correction.cpp give each operation order):
// 1. H positions h₁,₂ = (X, Y, Z) ∓ ½r (sin θ cos φ, sin θ sin φ, cos θ), with the unclamped r and Z.
// 2. With rescaling (EnableSurfExpan), the ideal-frame positions h′: the bond centre's x and y are
//    divided by s = a_T/a₀ and the atoms rebuilt from centre and bond (z changes by rounding);
//    without rescaling h′ = h.
// 3. Per layer, the pair-function parameters p₁…p₈ = slope·t + intercept at t = r clamped to
//    [P17, P18], and C = f(p₈).
// 4. ΔV = +0; for each atom j in file order, for H atom 1 then 2:
//    ΔV ← (ΔV − f(|h′ − R_j^ideal|)) + f(|h − R_j^distorted|), with the parameters of j's layer.
// The pair function f(R), q = R − p₁: C for R > p₈; else −(p₃ + q(p₄ + q(p₅ + q p₆)))·exp(−p₂ q),
// blended towards C for R > p₇ with ρ = ½cos(((R − p₈)π)/(p₈ − p₇)) + ½.
// Every term is kept, also far atoms with f = C: (ΔV − C) + C is not ΔV in floating point.

#include <array>
#include <optional>
#include <vector>

#include "potential/scm/scm_files.hpp"
#include "symmetry/surface_operation.hpp"

namespace spodvr {

// The pair function of one layer at one bond length: p₁…p₈ evaluated at the clamped r, and its
// value beyond the switch. Fortran: SCM_Types.f90 SCM_Parameters_Rydberg_3Body (P(1:8), CacheRyd),
// one element of MySCMPar
struct ScmPairTerms {
  double reference_distance = 0.0;     // p₁: q = R − p₁
  double decay = 0.0;                  // p₂: the exponent is −p₂ q
  std::array<double, 4> polynomial{};  // p₃…p₆: the cubic p₃ + q(p₄ + q(p₅ + q p₆))
  double switch_start = 0.0;           // p₇: f is blended towards C for R > p₇
  double switch_end = 0.0;             // p₈: f = C for R > p₈
  double switched_value = 0.0;         // C = f(p₈) without the switch. Fortran: CacheRyd
};

// The terms of a layer at bond length r (bohr): t = r clamped to [P17, P18] as the Fortran does
// it (IF/ELSE IF; a NaN r stays NaN), p_k = slope_k·t + intercept_k. The Fortran tabulates them at
// the r grid points rgrid(ir) and looks them up by ir; the same expressions give the same doubles,
// so these equal its cached terms only for r = rgrid(ir): the potential must pass the grid's r.
// Fortran: SCM_Types.f90 Setup_3Body_Rydberg_Par
[[nodiscard]] ScmPairTerms scm_pair_terms(const ScmLayerParameters& layer, double r);

// The pair function f(R) (hartree) at an H–metal distance R (bohr). Fortran: SCM_Potential.f90
// Switch_Rydberg_Potential
[[nodiscard]] double scm_pair_energy(const ScmPairTerms& terms, double h_metal_distance);

// |a − b| (bohr), summed as the Fortran's SQRT(SUM((a − b)**2)): √((dx² + dy²) + dz²).
// Fortran: SCM_Potential.f90 SCM_SurfInteract_3Body
[[nodiscard]] double distance(const Vec3& a, const Vec3& b);

// ΔV of one slab snapshot. Built from data in memory (the readers of scm_files.hpp give them);
// evaluation is const and thread-safe. Fortran: SCM_Wrapper.f90 GetSCMPotential with the module
// state MySCMData and MySCMPar
class ScmCorrection {
 public:
  // parameters: params.dat; slab: the snapshot (atoms of layers 1 and 2 only);
  // rescaling_factor: s = a_T/a₀ (the slab's over the CRP's lattice constant) if the slab is
  // rescaled to the CRP lattice (EnableSurfExpan), else empty; s = 1 is not the same as no
  // rescaling, since the ideal-frame round trip rounds. skewing_angle: γ of the CRP data set (rad).
  // Throws std::invalid_argument for an atom outside layers 1–2 or a factor that is not positive.
  // Fortran: SCM_Wrapper.f90 SetSCMCRPParameters (s), DistortedSurfaceWrapper.f90
  // GetCRPDistSurfParams (the matrix of γ)
  ScmCorrection(const ScmParameters& parameters, const ScmSlab& slab, std::optional<double> rescaling_factor,
                double skewing_angle);

  // The Cartesian position (bohr) of skewed lattice coordinates (u, v): (u + cos γ·v, sin γ·v).
  // Fortran: DistortedSurfaceWrapper.f90 GetCRPDistSurfPotential (MATMUL with uvToxyMatrix)
  [[nodiscard]] CartesianVector to_cartesian(double u, double v) const;

  // ΔV in hartree at the Cartesian centre of mass (x, y), height z, bond length r (bohr) and
  // polar and azimuthal angles θ, φ (rad). r is an r grid point, as in the Fortran (whose pair
  // terms are those of the grid index; see scm_pair_terms). Fortran: SCM_Wrapper.f90 GetSCMPotential
  [[nodiscard]] double value_hartree(double x_cartesian, double y_cartesian, double z, double r, double theta,
                                     double phi) const;

 private:
  std::array<ScmLayerParameters, 2> layers_;  // [0] layer 1, [1] layer 2
  std::vector<ScmAtom> atoms_;                // in file order
  std::optional<double> rescaling_factor_;    // s. Fortran: LatExpansion if EnableExpansion
  double cos_skew_ = 0.0;                     // cos γ. Fortran: uvToxyMatrix(1, 2)
  double sin_skew_ = 0.0;                     // sin γ. Fortran: uvToxyMatrix(2, 2)
};

}  // namespace spodvr

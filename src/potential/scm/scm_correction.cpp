// The SCM correction ΔV in the Fortran's operation order. Fortran: DistortedSurfacePES.dir/SCM_Potential.f90
// (SCM_SurfInteract_3Body, Switch_Rydberg_Potential), SCM_Wrapper.f90 (GetSCMPotential, ConvertCoMToXYZ),
// SCM_Types.f90 (Setup_3Body_Rydberg_Par), DistortedSurfaceWrapper.f90 (GetCRPDistSurfPotential: to_cartesian)
#include "potential/scm/scm_correction.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#include "core/constants.hpp"

namespace spodvr {

// SCM has its own π (SCM_Types.f90 Constants_PI = 3.14159265358979323846, made double by
// -fdefault-real-8): the same double as the program's.
static_assert(Constants::pi == 3.14159265358979323846);

namespace {

// The two H atoms of the molecule (bohr). Fortran: SCM_Potential.f90 SCM_SurfInteract_3Body H2Vec(1:3), H2Vec(4:6)
using HydrogenPair = std::array<Vec3, 2>;

// −(p₃ + q(p₄ + q(p₅ + q p₆)))·exp(−p₂ q): f without the switch. The Fortran writes C's product
// the other way round, exp·cubic: the same double.
double unswitched_pair_energy(const ScmPairTerms& terms, double q) {
  const std::array<double, 4>& c = terms.polynomial;
  const double cubic = c[0] + q * (c[1] + q * (c[2] + q * c[3]));
  return -(cubic * std::exp(-(terms.decay * q)));
}

// h₁,₂ = (X, Y, Z) ∓ ½r (sin θ cos φ, sin θ sin φ, cos θ), each offset the Fortran's left-to-right
// product ((0.5·r)·sin θ)·cos φ. Fortran: SCM_Wrapper.f90 ConvertCoMToXYZ
HydrogenPair hydrogen_positions(double x, double y, double z, double r, double theta, double phi) {
  const double sin_theta = std::sin(theta);
  const double cos_theta = std::cos(theta);
  const double sin_phi = std::sin(phi);
  const double cos_phi = std::cos(phi);
  const double half_r = 0.5 * r;
  const Vec3 offset{(half_r * sin_theta) * cos_phi, (half_r * sin_theta) * sin_phi, half_r * cos_theta};
  return {Vec3{x - offset.x, y - offset.y, z - offset.z}, Vec3{x + offset.x, y + offset.y, z + offset.z}};
}

// The H atoms in the frame of the ideal (CRP) lattice: the bond centre's x and y divided by s, the
// bond vector kept. The round trip through centre and bond changes the atoms, z too, by rounding;
// kept, as in the Fortran. Fortran: SCM_Potential.f90 SCM_SurfInteract_3Body (EnableExpansion)
HydrogenPair to_ideal_frame(const HydrogenPair& h, double s) {
  const Vec3 bond{h[1].x - h[0].x, h[1].y - h[0].y, h[1].z - h[0].z};
  const Vec3 centre{((h[0].x + h[1].x) / 2.0) / s, ((h[0].y + h[1].y) / 2.0) / s, (h[0].z + h[1].z) / 2.0};
  return {Vec3{centre.x - 0.5 * bond.x, centre.y - 0.5 * bond.y, centre.z - 0.5 * bond.z},
          Vec3{centre.x + 0.5 * bond.x, centre.y + 0.5 * bond.y, centre.z + 0.5 * bond.z}};
}

}  // namespace

ScmPairTerms scm_pair_terms(const ScmLayerParameters& layer, double r) {
  double t = r;  // Fortran: tR
  if (r < layer.r_clamp_min)
    t = layer.r_clamp_min;
  else if (r > layer.r_clamp_max)
    t = layer.r_clamp_max;
  std::array<double, 8> p{};
  for (int k = 0; k < 8; ++k) p[k] = layer.slope[k] * t + layer.intercept[k];
  ScmPairTerms terms{.reference_distance = p[0],
                     .decay = p[1],
                     .polynomial = {p[2], p[3], p[4], p[5]},
                     .switch_start = p[6],
                     .switch_end = p[7]};
  terms.switched_value = unswitched_pair_energy(terms, terms.switch_end - terms.reference_distance);
  return terms;
}

double scm_pair_energy(const ScmPairTerms& terms, double h_metal_distance) {
  if (h_metal_distance > terms.switch_end) return terms.switched_value;
  const double f = unswitched_pair_energy(terms, h_metal_distance - terms.reference_distance);
  if (!(h_metal_distance > terms.switch_start)) return f;
  // ρ rises from 0 at p₇ to 1 at p₈; its argument is ((R − p₈)·π)/(p₈ − p₇), in this order.
  const double phase =
      ((h_metal_distance - terms.switch_end) * Constants::pi) / (terms.switch_end - terms.switch_start);
  const double rho = 0.5 * std::cos(phase) + 0.5;
  return (1.0 - rho) * f + rho * terms.switched_value;
}

double distance(const Vec3& a, const Vec3& b) {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  const double dz = a.z - b.z;
  // The Fortran's SUM starts from 0, and 0 + dx² is dx² bitwise (a square is never −0): left out.
  return std::sqrt((dx * dx + dy * dy) + dz * dz);
}

ScmCorrection::ScmCorrection(const ScmParameters& parameters, const ScmSlab& slab,
                             std::optional<double> rescaling_factor, double skewing_angle)
    : layers_(parameters.layers),
      atoms_(slab.atoms),
      rescaling_factor_(rescaling_factor),
      cos_skew_(std::cos(skewing_angle)),
      sin_skew_(std::sin(skewing_angle)) {
  for (const ScmAtom& atom : atoms_)
    if (atom.layer < 1 || atom.layer > 2)
      throw std::invalid_argument("ScmCorrection: an atom of layer " + std::to_string(atom.layer) +
                                  "; the correction uses layers 1 and 2 only");
  if (rescaling_factor_ && !(*rescaling_factor_ > 0.0))
    throw std::invalid_argument("ScmCorrection: the rescaling factor must be positive");
}

// gfortran's MATMUL adds to a zero, x = (1·u + 0) + cos γ·v and y = (0 + 0·u) + sin γ·v. Left out:
// that changes only the sign of a zero x or y, which every distance squares away.
CartesianVector ScmCorrection::to_cartesian(double u, double v) const { return {u + cos_skew_ * v, sin_skew_ * v}; }

double ScmCorrection::value_hartree(double x_cartesian, double y_cartesian, double z, double r, double theta,
                                    double phi) const {
  // θ is the polar angle. The Fortran passes θ into GetSCMPotential's argument named phi and φ into
  // theta, and packs them back as (…, phi, theta) for ConvertCoMToXYZ: the two swaps cancel.
  const HydrogenPair h_distorted = hydrogen_positions(x_cartesian, y_cartesian, z, r, theta, phi);
  const HydrogenPair h_ideal = rescaling_factor_ ? to_ideal_frame(h_distorted, *rescaling_factor_) : h_distorted;
  const std::array<ScmPairTerms, 2> terms = {scm_pair_terms(layers_[0], r), scm_pair_terms(layers_[1], r)};
  double delta_v = +0.0;  // Fortran: potE = 0.0 (a slab without atoms gives +0, not −0)
  for (const ScmAtom& atom : atoms_) {
    const ScmPairTerms& pair = terms[atom.layer - 1];
    for (int a = 0; a < 2; ++a) {  // H atom 1, then 2
      const double f_ideal = scm_pair_energy(pair, distance(h_ideal[a], atom.ideal));
      const double f_distorted = scm_pair_energy(pair, distance(h_distorted[a], atom.distorted));
      delta_v = (delta_v - f_ideal) + f_distorted;  // the Fortran order, kept
    }
  }
  return delta_v;
}

}  // namespace spodvr

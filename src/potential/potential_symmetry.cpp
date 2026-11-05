// Surface symmetry of the potential: the grid symmetry of a run and the PES self-check.
// Fortran: CxvSurfaceSymmetry.f90 (the self-check is not in the Fortran)
#include "potential/potential_symmetry.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "symmetry/point_group.hpp"

namespace spodvr {

namespace {

// Fractional part of k * alpha + 0.5 (an additive recurrence: well spread, reproducible).
double sequence(int k, double alpha) { return std::fmod(0.5 + k * alpha, 1.0); }

// The irrational steps of the sequence in each sampled coordinate (golden ratio and relatives).
constexpr double alpha_x = 0.6180339887498949;
constexpr double alpha_y = 0.4142135623730951;
constexpr double alpha_cos_theta = 0.7320508075688772;
constexpr double alpha_phi = 0.2360679774997897;
constexpr double alpha_z = 0.6457513110645907;
constexpr double alpha_r = 0.3166247903554;

double clipped(double v_ev, const Config& config) {
  double v = v_ev;
  if (config.input.impose_v_min) v = std::max(v, config.input.v_min);
  if (config.input.impose_v_max) v = std::min(v, config.input.v_max);
  return v;
}

double wrap(double x, double period) { return x - period * std::floor(x / period); }

}  // namespace

GridSymmetry make_grid_symmetry(const Config& config, const AngularGrid& angular) {
  const InputParameters& p = config.input;
  if (!p.use_symmetry) return GridSymmetry(p.nx, p.ny, angular.n_phi());

  const PointGroup group(config.derived.system.point_group);
  GridSymmetry symmetry(SurfaceLattice(config.derived.system), group, p.nx, p.ny, angular.n_phi());
  if (symmetry.is_trivial()) {
    std::string message = "UseSymmetry: no operation of point group " + group.name() + " is usable";
    if (group.order() > 1)
      message += " with " + std::to_string(angular.n_phi()) + " phi points (a multiple of " +
                 std::to_string(symmetry.phi_points_multiple()) + " is needed; set nPhisOptimal)";
    throw InputError("UseSymmetry", message + ".");
  }
  return symmetry;
}

SymmetryCheck check_potential_symmetry(const Potential& potential, const GridSymmetry& symmetry, const Config& config,
                                       const SpatialGrids& grids, int n_samples) {
  const SurfaceLattice lattice(config.derived.system);
  const double a = lattice.lattice_constant_x();
  const double b = lattice.lattice_constant_y();
  const int nz = grids.nz_potential;
  const int nr = grids.r.size();

  SymmetryCheck check;
  check.n_samples = n_samples;
  check.tolerance = config.input.symmetry_tolerance;
  for (int k = 0; k < n_samples; ++k) {
    const double x = a * sequence(k, alpha_x);
    const double y = b * sequence(k, alpha_y);
    const double theta = std::acos(1.0 - 2.0 * sequence(k, alpha_cos_theta));
    const double phi = Constants::two_pi * sequence(k, alpha_phi);
    const int iz = std::min(nz - 1, static_cast<int>(nz * sequence(k, alpha_z)));
    const int ir = std::min(nr - 1, static_cast<int>(nr * sequence(k, alpha_r)));
    const double z = grids.z_scattering[iz];
    const double r = grids.r[ir];
    const double v = clipped(potential.value(x, y, z, r, theta, phi, iz, ir), config);

    for (const SurfaceOperation& g : symmetry.operations()) {
      if (g == identity_operation()) continue;
      const auto [gx, gy] = lattice.to_skewed(apply_to_cartesian(g, lattice.to_cartesian(x, y)));
      const double gphi = wrap(apply_to_azimuth(g, phi), Constants::two_pi);
      const double gv = clipped(potential.value(wrap(gx, a), wrap(gy, b), z, r, theta, gphi, iz, ir), config);
      const bool first = check.n_comparisons == 0;
      check.n_comparisons = check.n_comparisons + 1;
      const double deviation = std::abs(v - gv);
      if (first || deviation > check.max_deviation) {
        check.max_deviation = deviation;
        check.x = x;
        check.y = y;
        check.z = z;
        check.r = r;
        check.theta = theta;
        check.phi = phi;
        check.operation = g;
      }
    }
  }
  return check;
}

void require_symmetric_potential(const SymmetryCheck& check, const Config& config, const GridSymmetry& symmetry) {
  if (check.passed()) return;
  std::ostringstream message;
  message.precision(6);
  message << "the potential " << config.input.system_id << " is not symmetric under point group "
          << symmetry.point_group() << " to within SymmetryTolerance = " << check.tolerance
          << " eV: |V(q) - V(g q)| = " << check.max_deviation << " eV for g = " << describe(check.operation)
          << " at x = " << check.x << ", y = " << check.y << ", Z = " << check.z << ", r = " << check.r
          << " bohr, theta = " << rad_to_deg(check.theta) << ", phi = " << rad_to_deg(check.phi)
          << " deg. Raise SymmetryTolerance or set UseSymmetry to false.";
  throw InputError("SymmetryTolerance", message.str());
}

}  // namespace spodvr

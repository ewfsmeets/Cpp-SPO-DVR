// Checks of an SCM slab against its CRP data set (SCM D9: the ideal atoms on the CRP lattice), as
// data for the report plus a throwing require_*. No Fortran origin: the Fortran checks none of this.
#include "potential/scm/scm_diagnostics.hpp"

#include <array>
#include <cmath>
#include <format>

#include "core/error.hpp"

namespace spodvr {

namespace {

constexpr const char* position_file_key = "PosFilepath";

// The site of a layer: its offset in skewed lattice units (u, v) and its height (bohr).
struct LayerSite {
  std::array<double, 2> offset{};
  double z = 0.0;
};

// The distance (bohr) of `p` from the nearest site (i + offset_u, j + offset_v)·a at height site.z.
// The nearest site is found by rounding the skewed coordinates, which is exact for deviations far
// below a/2, the only ones that pass the check.
double distance_from_lattice(const Vec3& p, const LayerSite& site, double a, double cos_skew, double sin_skew) {
  const double v = p.y / (a * sin_skew);
  const double u = p.x / a - v * cos_skew;
  const double du = (u - site.offset[0]) - std::round(u - site.offset[0]);
  const double dv = (v - site.offset[1]) - std::round(v - site.offset[1]);
  const double dx = (du + dv * cos_skew) * a;
  const double dy = dv * sin_skew * a;
  const double dz = p.z - site.z;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// True if `deviation` replaces `largest` as the maximum: it is larger, or the first NaN (which is
// then kept, so that the check fails instead of passing).
bool is_new_maximum(double deviation, double largest) {
  if (std::isnan(largest)) return false;
  return std::isnan(deviation) || deviation > largest;
}

}  // namespace

double ScmLatticeCheck::max_deviation() const {
  if (std::isnan(top_layer_deviation) || std::isnan(second_layer_deviation)) return std::nan("");
  return top_layer_deviation > second_layer_deviation ? top_layer_deviation : second_layer_deviation;
}

bool ScmLatticeCheck::warning() const { return !(max_deviation() <= scm_lattice_warning_distance); }

ScmLatticeCheck check_slab_lattice(const ScmSlab& slab, const SystemParameters& crp,
                                   const CrpDataSetSummary& data_set) {
  const double a = crp.lattice_constant_x;
  const double cos_skew = std::cos(crp.skewing_angle);
  const double sin_skew = std::sin(crp.skewing_angle);
  const LayerSite top{{0.0, 0.0}, 0.0};
  const bool check_second = data_set.second_layer_site.has_value() && !data_set.interlayer_spacings.empty();
  const LayerSite second =
      check_second ? LayerSite{*data_set.second_layer_site, -data_set.interlayer_spacings[0]} : LayerSite{};

  ScmLatticeCheck check;
  double worst = -1.0;
  for (int k = 0; k < std::ssize(slab.atoms); ++k) {
    const ScmAtom& atom = slab.atoms[k];
    if (atom.layer == 2 && !check_second) continue;
    const double deviation = distance_from_lattice(atom.ideal, atom.layer == 1 ? top : second, a, cos_skew, sin_skew);
    if (atom.layer == 1) {
      ++check.top_layer_atoms;
      if (is_new_maximum(deviation, check.top_layer_deviation)) check.top_layer_deviation = deviation;
    } else {
      ++check.second_layer_atoms;
      if (is_new_maximum(deviation, check.second_layer_deviation)) check.second_layer_deviation = deviation;
    }
    if (is_new_maximum(deviation, worst)) {
      worst = deviation;
      check.worst_atom = k;
    }
  }
  return check;
}

void require_slab_lattice(const ScmLatticeCheck& check, const ScmSlab& slab, const SystemParameters& crp) {
  if (check.max_deviation() <= scm_lattice_error_distance) return;
  const ScmAtom& atom = slab.atoms.at(check.worst_atom);
  throw InputError(
      position_file_key,
      std::format("configuration {}: the ideal position ({}, {}, {}) bohr of an atom of layer {} lies {} bohr from "
                  "its site of the CRP lattice (a = {} bohr, skew angle {} rad; layer 1 on the lattice at z = 0, "
                  "layer 2 of a C3v data set one interlayer spacing lower under the site (1/3, 1/3)); more than "
                  "{} bohr: the slab's ideal lattice is not that of the data set CRPlibDataSet",
                  slab.configuration, atom.ideal.x, atom.ideal.y, atom.ideal.z, atom.layer, check.max_deviation(),
                  crp.lattice_constant_x, crp.skewing_angle, scm_lattice_error_distance));
}

}  // namespace spodvr

// The SCM set-up diagnostics (SCM D3: ΔV at V_asy; D9: the ideal atoms on the CRP lattice; D10: the
// cell-edge jump; D15: ΔV where the CRP clamps Z or r) and the slab's statistics, as data for the
// report plus throwing require_*. No Fortran origin: the Fortran checks none of this.
#include "potential/scm/scm_diagnostics.hpp"

#include <array>
#include <cmath>
#include <exception>
#include <format>
#include <functional>
#include <optional>
#include <string>
#include <string_view>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "potential/scm/scm_correction.hpp"
#include "potential/scm/scm_crp_potential.hpp"

namespace spodvr {

namespace {

constexpr std::string_view position_file_key = "PosFilepath";

// The (Z, r) of the cell-edge jumps (bohr): the barrier region of H₂/Cu(111) and above it, at the
// SRP48 barrier's bond length and at H₂'s equilibrium one (SCM D10; docs/scm_open_questions.md Q2).
constexpr std::array<double, 3> edge_heights = {2.25, 3.0, 4.0};
constexpr std::array<double, 2> edge_bond_lengths = {1.95, 1.4};
constexpr int edge_points = 12;      // points t = k/12 along each edge
constexpr int edge_phi_points = 12;  // φ = k·15°

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

// Which end of a scan is wanted.
enum class Extreme { largest, lowest };

// True if `value` replaces `best` (`none`: there is no best yet): a larger (lower) value, or the
// first NaN, which then stays. Ties keep the earlier point, so a scan in a fixed order ends on a
// fixed point.
bool replaces(Extreme kind, double value, double best, bool none) {
  if (none) return true;
  if (std::isnan(best)) return false;
  if (std::isnan(value)) return true;
  return kind == Extreme::largest ? value > best : value < best;
}

// Adds one scanned point to `extreme`.
void add_point(ScmExtreme& extreme, Extreme kind, double value, const ScmPoint& at) {
  if (replaces(kind, value, extreme.value, extreme.points == 0)) {
    extreme.value = value;
    extreme.at = at;
  }
  extreme.points = extreme.points + 1;
}

// Adds the result of a later part of a scan to `extreme`, as if its points had been added one by one.
void add_part(ScmExtreme& extreme, Extreme kind, const ScmExtreme& part) {
  if (part.points == 0) return;
  if (replaces(kind, part.value, extreme.value, extreme.points == 0)) {
    extreme.value = part.value;
    extreme.at = part.at;
  }
  extreme.points = extreme.points + part.points;
}

// The extreme over rows 0 … n_rows − 1, scan_row(row, part) scanning one row into `part`. Rows run in
// parallel and are combined in row order, so value and point do not depend on the threads. An
// exception of a row is rethrown after the loop (the first in row order).
ScmExtreme scan_rows(int n_rows, Extreme kind, const std::function<void(int, ScmExtreme&)>& scan_row) {
  std::vector<ScmExtreme> parts(n_rows);
  std::vector<std::exception_ptr> failures(n_rows);
  // dynamic: the cost of a row varies with its position (the CRP's and ΔV's work); each row writes
  // only its own slots, so the schedule does not change the result.
#pragma omp parallel for schedule(dynamic, 1)
  for (int row = 0; row < n_rows; ++row) {
    try {
      scan_row(row, parts[row]);
    } catch (...) {
      failures[row] = std::current_exception();
    }
  }
  for (const std::exception_ptr& failure : failures)
    if (failure) std::rethrow_exception(failure);
  ScmExtreme total;
  for (const ScmExtreme& part : parts) add_part(total, kind, part);
  return total;
}

// The two terms of the SCM potential, V = V_CRP(u/s, v/s) + ΔV (eV).
struct ScmTerms {
  const Potential& scm;
  const Potential& crp;
  double crp_scale;  // s, or 1 without rescaling

  // V at p (iz, ir: the grid indices passed on to the PES).
  [[nodiscard]] double value(const ScmPoint& p, int iz, int ir) const {
    return scm.value(p.u, p.v, p.z, p.r, p.theta, p.phi, iz, ir);
  }
  // ΔV as V − V_CRP(u/s, v/s): zero exactly where ΔV leaves V unchanged.
  [[nodiscard]] double delta_v(const ScmPoint& p, int iz, int ir) const {
    return value(p, iz, ir) - crp.value(p.u / crp_scale, p.v / crp_scale, p.z, p.r, p.theta, p.phi, iz, ir);
  }
};

// The angles a scan runs over (rad), in grid order.
struct ScanAngles {
  std::vector<double> thetas;
  std::vector<double> phis;
};

// Every θ and φ of the angular grid.
ScanAngles all_angles(const AngularGrid& angular) {
  return {{angular.thetas().begin(), angular.thetas().end()}, {angular.phis().begin(), angular.phis().end()}};
}

// Adds ΔV at every (x, y) of the grid and every angle at (z, r_ir), in grid order, to `part`: |ΔV|
// for the largest, ΔV for the lowest.
void scan_xy_and_angles(const ScmTerms& terms, const SpatialGrids& grids, const ScanAngles& angles, double z, int iz,
                        int ir, Extreme kind, ScmExtreme& part) {
  for (int ix = 0; ix < grids.x.size(); ++ix)
    for (int iy = 0; iy < grids.y.size(); ++iy)
      for (const double theta : angles.thetas)
        for (const double phi : angles.phis) {
          const ScmPoint p{grids.x[ix], grids.y[iy], z, grids.r[ir], theta, phi};
          const double delta_v = terms.delta_v(p, iz, ir);
          add_point(part, kind, kind == Extreme::largest ? std::abs(delta_v) : delta_v, p);
        }
}

// The largest |ΔV| on the plane Z = z_asymptotic at the r grid points below ROptStart (outside the r
// absorber); iz = nzpot marks the gas phase, as for V_asy (potential_grid.cpp).
ScmExtreme scan_asymptotic_plane(const ScmTerms& terms, const Config& config, const SpatialGrids& grids,
                                 const AngularGrid& angular) {
  std::vector<int> rows;  // the r indices
  for (int ir = 0; ir < grids.r.size(); ++ir)
    if (fortran_compare::rlt(grids.r[ir], config.input.r_opt_start)) rows.push_back(ir);
  return scan_rows(static_cast<int>(rows.size()), Extreme::largest, [&](int row, ScmExtreme& part) {
    scan_xy_and_angles(terms, grids, all_angles(angular), config.derived.z_asymptotic, grids.nz_potential, rows[row],
                       Extreme::largest, part);
  });
}

// The two lattice points (u, v) that the periodic grid joins at distance t along an edge: (t, 0)
// and (t, a) on the edge v = 0 (edge 0), (0, t) and (a, t) on the edge u = 0 (edge 1).
struct EdgePair {
  double u_inside, v_inside, u_image, v_image;
};

// The pair of edge 0 or 1 at distance t from its start; a: the lattice constant (bohr).
EdgePair edge_pair(int edge, double t, double a) { return edge == 0 ? EdgePair{t, 0.0, t, a} : EdgePair{0.0, t, a, t}; }

// A (Z, r) grid point of the 6D potential: one row of a hole scan.
struct ZrPoint {
  int iz = 0;
  int ir = 0;
};

// The extreme ΔV over every (x, y) and angle at the given (Z, r) grid points (rows, in this order).
ScmExtreme scan_zr_points(const ScmTerms& terms, const SpatialGrids& grids, const ScanAngles& angles,
                          const std::vector<ZrPoint>& rows, Extreme kind) {
  return scan_rows(static_cast<int>(rows.size()), kind, [&](int row, ScmExtreme& part) {
    const ZrPoint& zr = rows[row];
    scan_xy_and_angles(terms, grids, angles, grids.z_scattering[zr.iz], zr.iz, zr.ir, kind, part);
  });
}

// The angles of the scan outside the CRP's r range: every θ, but only the grid's first φ (0). Every
// (x, y) is scanned, so at each θ the lower H atom still passes over the whole cell. All φ would add
// ≈ 85 % to the 6D potential's time, one φ ≈ 9 % (the SCM example grid, 4 threads, 2026-11-14).
ScanAngles outside_r_range_angles(const AngularGrid& angular) {
  ScanAngles angles = all_angles(angular);
  angles.phis.resize(1);
  return angles;
}

// The r grid points where the r absorber acts, as it is applied (specular_operator.cpp
// apply_absorbers, specular_kernels.cu): SpatialGrids::r_optical (ROptStart ≤ r ≤ ROptEnd, tolerant
// comparisons; Fortran: Setup.f90 nrchsp .. nrchmx), none with nR = 1 or when no grid point lies
// between ROptStart and ROptEnd (first > last).
std::optional<IndexRange> r_absorber_points(const SpatialGrids& grids) {
  if (grids.r.size() <= 1 || grids.r_optical.last < grids.r_optical.first) return std::nullopt;
  return grids.r_optical;
}

// The jump across both cell edges at one (Z, r), θ = 90°, with `limit` the accessible limit (eV).
ScmEdgeJump edge_jump(const Potential& scm, double a, double z, double r, double limit) {
  ScmEdgeJump row{z, r};
  const double theta = Constants::pi / 2.0;
  for (int edge = 0; edge < 2; ++edge)
    for (int k = 0; k < edge_points; ++k) {
      const EdgePair pair = edge_pair(edge, a * k / edge_points, a);
      for (int m = 0; m < edge_phi_points; ++m) {
        const double phi = Constants::pi * m / edge_phi_points;
        const double inside = scm.value(pair.u_inside, pair.v_inside, z, r, theta, phi, 0, 0);
        const double image = scm.value(pair.u_image, pair.v_image, z, r, theta, phi, 0, 0);
        const double jump = std::abs(inside - image);
        if (is_new_maximum(jump, row.max_jump)) row.max_jump = jump;
        if (!(inside < limit && image < limit)) continue;
        row.accessible_pairs = row.accessible_pairs + 1;
        if (is_new_maximum(jump, row.max_accessible_jump)) row.max_accessible_jump = jump;
      }
    }
  return row;
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
      std::string(position_file_key),
      std::format("configuration {}: the ideal position ({}, {}, {}) bohr of an atom of layer {} lies {} bohr from "
                  "its site of the CRP lattice (a = {} bohr, skew angle {} rad; layer 1 on the lattice at z = 0, "
                  "layer 2 of a C3v data set one interlayer spacing lower under the site (1/3, 1/3)); more than "
                  "{} bohr: the slab's ideal lattice is not that of the data set CRPlibDataSet",
                  slab.configuration, atom.ideal.x, atom.ideal.y, atom.ideal.z, atom.layer, check.max_deviation(),
                  crp.lattice_constant_x, crp.skewing_angle, scm_lattice_error_distance));
}

ScmSlabStatistics slab_statistics(const ScmSlab& slab) {
  ScmSlabStatistics statistics;
  double sum_of_squares = 0.0;
  Vec3 top_layer_sum;
  for (int k = 0; k < std::ssize(slab.atoms); ++k) {
    const ScmAtom& atom = slab.atoms[k];
    const double displacement = distance(atom.distorted, atom.ideal);
    sum_of_squares = sum_of_squares + displacement * displacement;
    if (is_new_maximum(displacement, statistics.max_displacement)) statistics.max_displacement = displacement;
    if (k == 0 || is_new_maximum(atom.distorted.z, statistics.highest_atom)) statistics.highest_atom = atom.distorted.z;
    if (atom.layer != 1) {
      statistics.second_layer_atoms = statistics.second_layer_atoms + 1;
      continue;
    }
    statistics.top_layer_atoms = statistics.top_layer_atoms + 1;
    top_layer_sum = {top_layer_sum.x + (atom.distorted.x - atom.ideal.x),
                     top_layer_sum.y + (atom.distorted.y - atom.ideal.y),
                     top_layer_sum.z + (atom.distorted.z - atom.ideal.z)};
  }
  const int n_atoms = static_cast<int>(slab.atoms.size());
  if (n_atoms > 0) statistics.rms_displacement = std::sqrt(sum_of_squares / n_atoms);
  const int n_top = statistics.top_layer_atoms;
  if (n_top > 0)
    statistics.top_layer_shift = {top_layer_sum.x / n_top, top_layer_sum.y / n_top, top_layer_sum.z / n_top};
  return statistics;
}

bool ScmAsymptoticCheck::plane_warning() const {
  return largest_on_plane.has_value() && !(largest_on_plane->value <= scm_plane_warning_energy);
}

ScmAsymptoticCheck check_asymptotic_potential(const Potential& scm, const Potential& crp, double crp_scale,
                                              const Config& config, const SpatialGrids& grids,
                                              const AngularGrid& angular, bool with_plane) {
  const ScmTerms terms{scm, crp, crp_scale};
  ScmAsymptoticCheck check;
  check.z_asymptotic = config.derived.z_asymptotic;
  check.point_heights = {check.z_asymptotic};
  const double switch_height = config.derived.z_asymptotic_unswitched;  // V_gas of switch_to_gas_phase
  if (config.input.use_switching && switch_height != check.z_asymptotic) check.point_heights.push_back(switch_height);
  check.plane_r_limit = config.input.r_opt_start;
  ScmExtreme gas_phase;
  // The points of V_gas (potential_grid.cpp gas_phase_potential), with iz = nzpot marking the gas phase.
  for (const double z : check.point_heights)
    for (int ir = 0; ir < grids.r.size(); ++ir) {
      const ScmPoint p{0.0, 0.0, z, grids.r[ir], Constants::pi / 2.0, 0.0};
      const double delta_v = terms.delta_v(p, grids.nz_potential, ir);
      if (!(delta_v == 0.0)) check.changed_points = check.changed_points + 1;
      add_point(check.largest_at_points, Extreme::largest, std::abs(delta_v), p);
      if (z == check.z_asymptotic) add_point(gas_phase, Extreme::lowest, terms.value(p, grids.nz_potential, ir), p);
    }
  check.gas_phase_minimum = gas_phase.value;
  if (with_plane) check.largest_on_plane = scan_asymptotic_plane(terms, config, grids, angular);
  return check;
}

void require_unchanged_asymptotic_potential(const ScmAsymptoticCheck& check, const ScmSlab& slab,
                                            const ScmSlabStatistics& statistics) {
  if (check.changed_points == 0) return;
  const ScmExtreme& worst = check.largest_at_points;
  throw InputError(
      std::string(position_file_key),
      std::format("configuration {}: ΔV changes the gas-phase potential V_asy at {} of its {} points (0, 0, Z, r, "
                  "θ = π/2, φ = 0), by up to {} eV at Z = {} bohr, r = {} bohr: a slab atom reaches the gas phase "
                  "(the highest lies at z = {} bohr). V_asy, the gas-phase states and the switch to the gas phase "
                  "must not depend on the slab (SCM D3); raise z_asymptotic (ZmaxSwitch with UseSwitching) or "
                  "choose another snapshot (PosFileSlabNum) or a smaller SCMCutoff",
                  slab.configuration, check.changed_points, worst.points, worst.value, worst.at.z, worst.at.r,
                  statistics.highest_atom));
}

ScmEdgeJumps check_edge_jumps(const Potential& scm, double accessible_limit) {
  ScmEdgeJumps jumps{accessible_limit, edge_points, edge_phi_points, {}};
  const double a = scm.system_parameters().lattice_constant_x;
  for (const double z : edge_heights)
    for (const double r : edge_bond_lengths) jumps.rows.push_back(edge_jump(scm, a, z, r, accessible_limit));
  return jumps;
}

bool crp_clamps_z(double z, const CrpDataSetSummary& data_set) { return z < data_set.z_min + crp_clamp_margin; }

bool crp_clamps_r(double r, const CrpDataSetSummary& data_set) {
  return r < data_set.r_min + crp_clamp_margin || r > data_set.r_max - crp_clamp_margin;
}

ScmHoleCheck check_holes_below_z_min(const Potential& scm, const Potential& crp, double crp_scale,
                                     const CrpDataSetSummary& data_set, const SpatialGrids& grids,
                                     const AngularGrid& angular) {
  ScmHoleCheck check{data_set.z_min, data_set.r_min, data_set.r_max, grids.z_scattering[0], 0, {}};
  while (check.planes < grids.nz_potential && crp_clamps_z(grids.z_scattering[check.planes], data_set))
    check.planes = check.planes + 1;
  // The r the CRP does not clamp; the others are check_holes_outside_r_range's.
  std::vector<ZrPoint> rows;
  for (int iz = 0; iz < check.planes; ++iz)
    for (int ir = 0; ir < grids.r.size(); ++ir)
      if (!crp_clamps_r(grids.r[ir], data_set)) rows.push_back({iz, ir});
  check.lowest = scan_zr_points(ScmTerms{scm, crp, crp_scale}, grids, all_angles(angular), rows, Extreme::lowest);
  return check;
}

ScmOutsideRRangeCheck check_holes_outside_r_range(const Potential& scm, const Potential& crp, double crp_scale,
                                                  const CrpDataSetSummary& data_set, const SpatialGrids& grids,
                                                  const AngularGrid& angular) {
  const ScanAngles angles = outside_r_range_angles(angular);
  ScmOutsideRRangeCheck check;
  check.r_min = data_set.r_min;
  check.r_max = data_set.r_max;
  check.theta_points = static_cast<int>(angles.thetas.size());
  check.phi_points = static_cast<int>(angles.phis.size());
  // The r the CRP clamps, on every plane of the 6D potential (with check_holes_below_z_min: every
  // clamped point once).
  std::vector<ZrPoint> rows;
  for (int iz = 0; iz < grids.nz_potential; ++iz)
    for (int ir = 0; ir < grids.r.size(); ++ir)
      if (crp_clamps_r(grids.r[ir], data_set)) rows.push_back({iz, ir});
  check.lowest = scan_zr_points(ScmTerms{scm, crp, crp_scale}, grids, angles, rows, Extreme::lowest);
  const std::optional<IndexRange> absorber = r_absorber_points(grids);
  if (!absorber) return check;
  check.r_absorber_points = absorber->size();
  check.r_absorber_first = grids.r[absorber->first];
  check.r_absorber_last = grids.r[absorber->last];
  for (int ir = absorber->first; ir <= absorber->last; ++ir)
    if (check.lowest.points > 0 && grids.r[ir] == check.lowest.at.r) check.in_r_absorber = true;
  return check;
}

ScmDiagnostics diagnose_scm_potential(const ScmCrpPotential& potential, const Config& config, const SpatialGrids& grids,
                                      const AngularGrid& angular) {
  const CrpPotential& crp = potential.crp();
  const double crp_scale = potential.rescaling_factor().value_or(1.0);
  const bool has_6d_potential = config.derived.run_mode != 2;  // none in an analysis-only run
  ScmDiagnostics diagnostics;
  diagnostics.crp_lattice_constant = crp.system_parameters().lattice_constant_x;
  diagnostics.rescaling_factor = potential.rescaling_factor();
  diagnostics.slab = slab_statistics(potential.slab());
  diagnostics.lattice = potential.lattice_check();
  diagnostics.asymptotic =
      check_asymptotic_potential(potential, crp, crp_scale, config, grids, angular, has_6d_potential);
  if (!has_6d_potential) return diagnostics;
  const double accessible_limit =
      diagnostics.asymptotic.gas_phase_minimum + config.input.e_max_wave_packet + scm_accessible_margin;
  diagnostics.edge_jumps = check_edge_jumps(potential, accessible_limit);
  const CrpDataSetSummary& data_set = crp.summary();
  diagnostics.holes = check_holes_below_z_min(potential, crp, crp_scale, data_set, grids, angular);
  diagnostics.outside_r_range = check_holes_outside_r_range(potential, crp, crp_scale, data_set, grids, angular);
  return diagnostics;
}

}  // namespace spodvr

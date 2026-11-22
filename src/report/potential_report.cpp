// Set-up report of the PES: the CRP data set and the SCM block with its diagnostics (SCM D3, D9, D10,
// D15). No Fortran origin: the Fortran prints only the CRP library's messages.
#include "report/potential_report.hpp"

#include <span>
#include <string>
#include <vector>

#include "core/unit_conversion.hpp"
#include "potential/crp_potential.hpp"
#include "potential/scm/scm_crp_potential.hpp"
#include "report/crp_report.hpp"

namespace spodvr {

namespace {

// The CRP data set and the grid points outside its range: the whole r grid and the Z points of the
// 6D potential (none in analysis-only runs).
void report_crp_data_set(std::ostream& out, const Style& style, const CrpPotential& crp, const Config& config,
                         const SpatialGrids& grids) {
  const bool analysis_only = config.derived.run_mode == 2;
  const std::span<const double> z_6d = grids.z_scattering.points().first(analysis_only ? 0 : grids.nz_potential);
  print_crp_data_set(out, style, crp.summary(), crp.range_report(z_6d, grids.r.points()));
}

// "(u, v) = (…, …) bohr, r = … bohr, theta = …, phi = … deg" of a scan's point (ASCII, as the rest of
// the report: the labels are aligned by their length in bytes).
std::string point_text(const ScmPoint& p) {
  return "(u, v) = (" + fixed(p.u, 4) + ", " + fixed(p.v, 4) + ") bohr, r = " + fixed(p.r, 4) +
         " bohr, theta = " + fixed(rad_to_deg(p.theta), 1) + ", phi = " + fixed(rad_to_deg(p.phi), 1) + " deg";
}

// The files, the snapshot, the lattice constants and the displacements.
void add_slab(KeyValueList& list, const ScmInput& input, const ScmSlab& slab, const ScmDiagnostics& d) {
  const ScmSlabHeader& header = slab.header;
  list.add("Parameter file (SCMParFilepath)", input.parameter_file);
  list.add("Position file (PosFilepath)", input.position_file);
  list.add("Slab tag", header.tag);
  list.add("Structure, file version, snapshots", header.structure + ", " + std::to_string(header.version) + ", " +
                                                     std::to_string(header.configuration_count));
  list.add("Snapshot used (PosFileSlabNum)", std::to_string(slab.configuration));
  list.add("Cut-off (SCMCutoff)", fixed(slab.cutoff, 6) + (header.version == 2 ? " (not applied to version 2)" : ""),
           "bohr");
  list.add("Atoms used (layer 1 + layer 2)",
           std::to_string(d.slab.top_layer_atoms) + " + " + std::to_string(d.slab.second_layer_atoms));
  list.add("Lattice constants a_T (slab), a0 (CRP)",
           full_precision(header.lattice_constant) + ", " + full_precision(d.crp_lattice_constant), "bohr");
  list.add("Rescaling (EnableSurfExpan)", d.rescaling_factor ? "on, s = a_T/a0 = " + full_precision(*d.rescaling_factor)
                                                             : "off (the grid has the CRP's a0)");
  list.add("Displacement rms, largest", fixed(d.slab.rms_displacement, 6) + ", " + fixed(d.slab.max_displacement, 6),
           "bohr");
  const Vec3& shift = d.slab.top_layer_shift;
  list.add("Mean displacement of layer 1 (x, y, z)",
           "(" + fixed(shift.x, 6) + ", " + fixed(shift.y, 6) + ", " + fixed(shift.z, 6) + ")", "bohr");
  list.add("Highest atom (distorted z)", fixed(d.slab.highest_atom, 6), "bohr");
}

// Where the lowest point outside the CRP's r range lies: below Rmin, or beyond Rmax inside or outside
// the r absorber (or with no r grid point in the absorber).
std::string r_absorber_text(const ScmOutsideRRangeCheck& check) {
  if (check.lowest.at.r < check.r_min) return "below Rmin";
  if (check.r_absorber_points == 0) return "beyond Rmax; no r grid point lies in the r absorber";
  return std::string("beyond Rmax, ") + (check.in_r_absorber ? "inside" : "outside") +
         " the r absorber (r = " + fixed(check.r_absorber_first, 4) + " .. " + fixed(check.r_absorber_last, 4) +
         " bohr)";
}

// SCM D15, the r region: the lowest V - V_CRP outside the CRP's r range, its point and where it lies.
void add_outside_r_range(KeyValueList& list, const ScmOutsideRRangeCheck& check) {
  const std::string label = "r outside Rmin .. Rmax, every Z plane (SCM D15)";
  const std::string range = "Rmin .. Rmax = " + fixed(check.r_min, 6) + " .. " + fixed(check.r_max, 6) + " bohr: ";
  if (check.lowest.points == 0) {
    list.add(label, range + "no r grid point outside");
    return;
  }
  const ScmExtreme& lowest = check.lowest;
  const std::string scanned = std::to_string(lowest.points) + " points: " + std::to_string(check.theta_points) +
                              " theta x " + std::to_string(check.phi_points) + " phi of the grid";
  list.add(label, range + "lowest V - V_CRP " + scientific(lowest.value, 3) + " eV at Z = " + fixed(lowest.at.z, 4) +
                      " bohr, " + point_text(lowest.at) + ", " + r_absorber_text(check) + "; " + scanned);
}

// SCM D9, D3 and D15 in one line each (D15: one per clamp region).
void add_checks(KeyValueList& list, const ScmDiagnostics& d) {
  list.add("Ideal atoms off the CRP lattice (SCM D9)",
           scientific(d.lattice.max_deviation(), 2) + " bohr at most (" + std::to_string(d.lattice.top_layer_atoms) +
               " + " + std::to_string(d.lattice.second_layer_atoms) + " atoms checked)");
  const ScmAsymptoticCheck& asy = d.asymptotic;
  std::string heights;
  for (const double z : asy.point_heights) heights += (heights.empty() ? "" : ", ") + fixed(z, 6);
  list.add("Gas-phase points with V != V_CRP (SCM D3)", std::to_string(asy.changed_points) + " of " +
                                                            std::to_string(asy.largest_at_points.points) +
                                                            " (Z = " + heights + " bohr)");
  if (asy.largest_on_plane) {
    const ScmExtreme& plane = *asy.largest_on_plane;
    list.add("Largest |V - V_CRP| on Z = " + fixed(asy.z_asymptotic, 4) + ", r < " + fixed(asy.plane_r_limit, 4),
             plane.points == 0 ? "no r point below ROptStart"
                               : scientific(plane.value, 3) + " eV at " + point_text(plane.at) + " (" +
                                     std::to_string(plane.points) + " points)");
  }
  if (!d.holes) return;
  const ScmHoleCheck& holes = *d.holes;
  list.add("Z planes below the CRP's Zmin (SCM D15)",
           holes.planes == 0
               ? "none (ZStart " + fixed(holes.z_start, 6) + " >= Zmin " + fixed(holes.z_min, 6) + " bohr)"
               : std::to_string(holes.planes) +
                     ", lowest V - V_CRP at Rmin <= r <= Rmax: " + scientific(holes.lowest.value, 3) +
                     " eV at Z = " + fixed(holes.lowest.at.z, 4) + " bohr, " + point_text(holes.lowest.at));
  if (d.outside_r_range) add_outside_r_range(list, *d.outside_r_range);
}

// SCM D10: one row per (Z, r).
void print_edge_jumps(std::ostream& out, const Style& style, const ScmEdgeJumps& jumps) {
  out << "\n  Jump of V across the cell edges (SCM D10: the periodic grid repeats one distorted cell) at\n"
      << "  theta = 90 deg, " << jumps.points_per_edge << " points per edge and " << jumps.phi_points
      << " phi; accessible: V below " << fixed(jumps.accessible_limit, 6) << " eV on both sides\n"
      << "  (min V_asy + EMaxWP + " << fixed(scm_accessible_margin, 1) << " eV)\n\n";
  Table table({"Z (bohr)", "r (bohr)", "max jump (eV)", "accessible pairs", "max accessible jump (eV)"});
  const int pairs = 2 * jumps.points_per_edge * jumps.phi_points;
  for (const ScmEdgeJump& row : jumps.rows)
    table.add_row({fixed(row.z, 2), fixed(row.r, 2), fixed(row.max_jump, 6),
                   std::to_string(row.accessible_pairs) + " of " + std::to_string(pairs),
                   row.accessible_pairs > 0 ? fixed(row.max_accessible_jump, 6) : "-"});
  table.print(out, style);
}

// A warning for each check that calls for one (SCM D9, D3, D15), after a blank line.
void print_scm_warnings(std::ostream& out, const Style& style, const ScmDiagnostics& d) {
  std::vector<std::string> warnings;
  if (d.lattice.warning())
    warnings.push_back("an ideal atom lies " + scientific(d.lattice.max_deviation(), 2) +
                       " bohr from its site of the CRP lattice, more than " +
                       scientific(scm_lattice_warning_distance, 0) + " bohr (SCM D9)");
  if (d.asymptotic.plane_warning())
    warnings.push_back("|V - V_CRP| reaches " + scientific(d.asymptotic.largest_on_plane->value, 3) +
                       " eV on the plane Z = z_asymptotic (r < ROptStart), more than " +
                       scientific(scm_plane_warning_energy, 0) + " eV: the gas phase still feels the slab (SCM D3)");
  if (d.holes && d.holes->warning())
    warnings.push_back("the Z grid starts at ZStart = " + fixed(d.holes->z_start, 6) +
                       " bohr, below the CRP's Zmin = " + fixed(d.holes->z_min, 6) +
                       " bohr, where\n  the CRP is evaluated at Zmin but the SCM term at the grid's Z: V - V_CRP "
                       "reaches " +
                       scientific(d.holes->lowest.value, 3) +
                       " eV at Rmin <= r <= Rmax\n  (SCM D15). Use ZStart >= Zmin, or ImposeVMin.");
  if (!warnings.empty()) out << '\n';
  for (const std::string& w : warnings) out << "  " << style.warning("Warning") << ": " << w << '\n';
}

}  // namespace

void print_scm_potential(std::ostream& out, const Style& style, const ScmInput& input, const ScmSlab& slab,
                         const ScmDiagnostics& diagnostics) {
  print_section(out, style, "SCM DISTORTED SURFACE (CRPDistrtSurface)");
  KeyValueList list;
  add_slab(list, input, slab, diagnostics);
  add_checks(list, diagnostics);
  list.print(out, style);
  if (diagnostics.edge_jumps)
    print_edge_jumps(out, style, *diagnostics.edge_jumps);
  else
    out << "\n  Analysis-only run: no 6D potential, so only the gas-phase points are checked.\n";
  print_scm_warnings(out, style, diagnostics);
}

void report_potential(std::ostream& out, const Style& style, const Potential& potential, const Config& config,
                      const SpatialGrids& grids, const AngularGrid& angular) {
  // The CRP data set (CRPlibrary, or the CRP part of CRPDistrtSurface), then the SCM correction.
  if (const auto crp = potential.crp_potential()) report_crp_data_set(out, style, crp->get(), config, grids);
  const auto scm_potential = potential.scm_potential();
  if (!scm_potential) return;
  const ScmCrpPotential& scm = scm_potential->get();
  const ScmDiagnostics diagnostics = diagnose_scm_potential(scm, config, grids, angular);
  print_scm_potential(out, style, config.input.scm, scm.slab(), diagnostics);
  require_unchanged_asymptotic_potential(diagnostics.asymptotic, scm.slab(), diagnostics.slab);
}

}  // namespace spodvr

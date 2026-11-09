// Setup report of the CRP PES: the data set and the grid points outside its range.
#include "report/crp_report.hpp"

#include <string>

namespace spodvr {

namespace {

std::string count_of(long long count, long long total) {
  return std::to_string(count) + " of " + std::to_string(total);
}

void add_range(KeyValueList& list, const CrpDataSetSummary& s, const CrpRangeReport& range) {
  list.add("Grid r points clamped to Rmin, Rmax", count_of(range.r_clamped, range.r_points));
  if (range.z_points == 0) return;
  list.add("Grid Z points clamped to Zmin (6D part)", count_of(range.z_clamped, range.z_points));
  list.add("(Z, r) points without data (V = " + fixed(s.energy_shift, 6) + " eV)",
           count_of(range.zr_uncovered, range.zr_points));
}

}  // namespace

void print_crp_data_set(std::ostream& out, const Style& style, const CrpDataSetSummary& s,
                        const CrpRangeReport& range) {
  print_section(out, style, "CRP DATA SET");
  KeyValueList list;
  list.add("Directory", s.directory);
  if (s.description) list.add("Description", *s.description);
  list.add("Symmetry", s.symmetry);
  list.add("Units of the files", s.file_units);
  list.add("Surface layers", std::to_string(s.surface_layers));
  list.add("Atom-surface cuts", std::to_string(s.atom_cuts));
  list.add("Molecule-surface cuts (used)",
           std::to_string(s.molecule_cuts) + " (" + std::to_string(s.molecule_cuts_used) + ")");
  list.add("The cuts hold", s.cuts_hold_interpolation ? "the interpolation function I (CutsInterpolation)"
                                                      : "V (reduced by the atom-surface PES on reading)");
  list.add("r covered by the cuts (Rmin .. Rmax)", fixed(s.r_min, 6) + " .. " + fixed(s.r_max, 6), "bohr");
  list.add("Z covered by the cuts (from Zmin)", fixed(s.z_min, 6), "bohr");
  list.add("Switch to gas phase (ZsoMin .. ZsoMax)", fixed(s.zso_min, 6) + " .. " + fixed(s.zso_max, 6), "bohr");
  list.add("Z asymptotic of the data (Zasy)", fixed(s.z_asymptotic, 6), "bohr");
  list.add("Energy shift (-VH2 at the last r)", full_precision(s.energy_shift), "eV");
  add_range(list, s, range);
  list.print(out, style);

  out << "\n  r and Z outside the data are evaluated at the nearest end (Rmin + 1e-12, Rmax - 1e-12,\n"
         "  Zmin + 1e-12 bohr), as in the Fortran.\n";
  if (range.zr_uncovered > 0) {
    out << "  " << style.warning("Warning") << ": at " << range.zr_uncovered
        << " (Z, r) points the data set has no data; the Fortran library silently returns V = 0 there\n";
    out << "  (shifted: " << fixed(s.energy_shift, 6) << " eV). Check the r and Z ranges of the cuts.\n";
  }
  for (const std::string& w : s.warnings) out << "  " << style.warning("Warning") << ": " << w << '\n';
}

}  // namespace spodvr

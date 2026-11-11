# Mutations of the CRP library (docs/crp_library.md "Mutation check"): each must make a CrpLib test fail.
# Usage: python3 -I mutations.py <tree> <index>   (applies mutation <index> to <tree>; prints its name)
#        python3 -I mutations.py list
import sys

MUTATIONS = [
    ("bridge images swapped", "src/crplib/hexagonal_interaction.cpp",
     "  s[4] = hexagonal_bridge(bridge, t, phi + pis3, want);\n  s[5] = hexagonal_bridge(bridge, t, phi - pis3, want);",
     "  s[4] = hexagonal_bridge(bridge, t, phi - pis3, want);\n  s[5] = hexagonal_bridge(bridge, t, phi + pis3, want);"),
    ("0.5 instead of the folded cos(pi/3) (QD wrapper skew)", "src/crplib/qd_wrapper.cpp",
     "cos_skew_ = square ? detail::cos_pi_2 : detail::cos_pi_3;",
     "cos_skew_ = square ? detail::cos_pi_2 : 0.5;"),
    ("0.5 instead of the folded cos(pi/3) (C3v MURII)", "src/crplib/hexagonal_atom_surface.cpp",
     "const double xijk = uuc + vuc * cos_pi_3;",
     "const double xijk = uuc + vuc * 0.5;"),
    ("SPLINT slope without its curvature term", "src/crplib/spline.cpp",
     "(y_[khi] - y_[klo]) / h + h * (m_[khi] * (3.0 * b * b - 1.0) - m_[klo] * (3.0 * a * a - 1.0)) / 6.0;",
     "(y_[khi] - y_[klo]) / h;"),
    ("SPLINT instead of DPL for the atom cuts", "src/crplib/hexagonal_atom_surface.cpp",
     "vs[i] = cuts_[i].value_dpl(p.z);",
     "vs[i] = cuts_[i].value_and_slope_splint(p.z).value;"),
    ("x / c instead of x * (1 / c) (Angstrom to bohr)", "src/crplib/units.cpp",
     "  if (from == to) return x;\n  return x * length_conversion(from, to);",
     "  if (from == to) return x;\n  if (from == LengthUnit::angstrom && to == LengthUnit::bohr) return x / conv_angstrom;\n"
     "  return x * length_conversion(from, to);"),
    ("(x * (1/a)) * b instead of x * ((1/a) * b)", "src/crplib/units.cpp",
     "  if (from == to) return x;\n  return x * energy_conversion(from, to);",
     "  if (from == to) return x;\n  if (to != EnergyUnit::hartree) return (x * to_hartree[code(from)]) * hartree_to[code(to)];\n"
     "  return x * energy_conversion(from, to);"),
    ("ORIGEN truncation: floor instead of INT(ori - 1)", "src/crplib/atom_surface.cpp",
     "return ori >= 0.0 ? static_cast<int>(ori) : static_cast<int>(ori - 1.0);",
     "return static_cast<int>(std::floor(ori));"),
    ("C3v t2h and t2f cuts swapped", "src/crplib/hexagonal_interaction.cpp",
     "    s[2] = hexagonal_hollow(hollow, t, phi, want);  // hcp\n    const auto t2h = cuts.subspan(12, 6);",
     "    s[2] = hexagonal_hollow(hollow, t, phi, want);  // hcp\n    const auto t2h = cuts.subspan(18, 6);"),
    ("sign flip in a chain-rule term (dphi/dx of atom a)", "src/crplib/molecule_surface.cpp",
     "d.theta_x * i_theta + d.phi_x * i_phi,",
     "d.theta_x * i_theta - d.phi_x * i_phi,"),
    ("TRANSF regions 4 and 5 swapped", "src/crplib/atom_surface.cpp",
     "    case 4:\n      aux = dvdx;\n      dvdx = -dvdy;\n      dvdy = aux;\n      return;\n"
     "    case 5:\n      aux = dvdx;\n      dvdx = dvdy;\n      dvdy = -aux;",
     "    case 5:\n      aux = dvdx;\n      dvdx = -dvdy;\n      dvdy = aux;\n      return;\n"
     "    case 4:\n      aux = dvdx;\n      dvdx = dvdy;\n      dvdy = -aux;"),
    ("dFso term dropped (atom a)", "src/crplib/molecule_surface.cpp",
     "full_a.z * sw.f + (va.v + vb.v + i_v) * df_a - df_a * gp.v + ",
     "full_a.z * sw.f + "),
    # Speed-up S2 (cluster pruning): real bugs, which the tests must catch ...
    ("S2: reachable atoms in reverse order", "src/crplib/atom_surface.cpp",
     "    if (distance_to_polygon(atom.x, atom.y, cell) <= reach_ + 2.0 * cell_margin_) atoms_.push_back(atom);\n}",
     "    if (distance_to_polygon(atom.x, atom.y, cell) <= reach_ + 2.0 * cell_margin_) atoms_.push_back(atom);\n"
     "  std::reverse(atoms_.begin(), atoms_.end());\n}"),
    ("S2: in-plane selection with half the reach", "src/crplib/atom_surface.cpp",
     "<= reach_ + 2.0 * cell_margin_) atoms_.push_back(atom);",
     "<= 0.5 * reach_) atoms_.push_back(atom);"),
    ("S2: per-call z test with half the reach (C3v, C6v)", "src/crplib/hexagonal_atom_surface.cpp",
     "    if (std::abs(z - atom.z) > z_reach) continue;",
     "    if (std::abs(z - atom.z) > 0.5 * z_reach) continue;"),
    ("S2: no fallback to the full cluster outside the fold cell (C3v, C6v)", "src/crplib/hexagonal_atom_surface.cpp",
     "  if (!in_fold_cell(x, y)) return repulsive_sum(x, y, z, want);\n",
     ""),
    # ... and the safety margins, each of which alone covers the rounding (SAFETY MARGIN: may survive).
    ("S2: no margins at all [SAFETY MARGIN]", "src/crplib/atom_surface.cpp",
     "    : reach_(cutoff > 0.0 ? cutoff * (1.0 + 1.0e-10) : std::numeric_limits<double>::infinity()),\n"
     "      cell_margin_(1.0e-8 * lattice_constant) {",
     "    : reach_(cutoff > 0.0 ? cutoff : std::numeric_limits<double>::infinity()),\n"
     "      cell_margin_(0.0 * lattice_constant) {"),
    ("S2: no relative margin [SAFETY MARGIN]", "src/crplib/atom_surface.cpp",
     "cutoff * (1.0 + 1.0e-10)", "cutoff"),
    ("S2: no cell enlargement [SAFETY MARGIN]", "src/crplib/atom_surface.cpp",
     "cell_margin_(1.0e-8 * lattice_constant)", "cell_margin_(0.0 * lattice_constant)"),
    ("S2: < instead of <= in the in-plane selection [SAFETY MARGIN]", "src/crplib/atom_surface.cpp",
     "<= reach_ + 2.0 * cell_margin_)", "< reach_ + 2.0 * cell_margin_)"),
    # Speed-up S4 (one grid lookup per group of cuts with the same grid)
    ("S4: every cut evaluated in the cell of the first group", "src/crplib/molecule_surface.cpp",
     "values[k] = cuts_[k].evaluate_in_cell(cells[group_of_cut_[k]], want);",
     "values[k] = cuts_[k].evaluate_in_cell(cells[0], want);"),
    ("S4: grids compared with != (merges -0 and +0 knots)", "src/crplib/molecule_surface.cpp",
     "if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;",
     "if (a[i] != b[i]) return false;"),
    ("S4: only the first group located", "src/crplib/molecule_surface.cpp",
     "for (std::size_t g = 0; g < grid_cut_.size(); ++g) {\n    const std::optional<detail::GridCell> cell",
     "for (std::size_t g = 0; g < 1; ++g) {\n    const std::optional<detail::GridCell> cell"),
    # Speed-up S3 (V without the derivatives): things the value needs, wrongly skipped
    ("S3: sin 3 phi of the hollow site skipped for V", "src/crplib/site_interpolation.cpp",
     "const double sin3phi = std::sin(3.00 * phi), cos6phi = std::cos(6.00 * phi);",
     "const double sin3phi = want == Want::value ? 0.0 : std::sin(3.00 * phi), cos6phi = std::cos(6.00 * phi);"),
    ("S3: the C3v sine stars of the atom PES skipped for V", "src/crplib/hexagonal_atom_surface.cpp",
     "const bool star_sines = c3v || gradient;", "const bool star_sines = gradient;"),
    ("S3: V-only blend without the gas-phase term", "src/crplib/molecule_surface.cpp",
     "return MoleculeSurfaceValue{(va + vb + i_v) * sw.f + (1.00 - sw.f) * gp.v, {}, {}};",
     "return MoleculeSurfaceValue{(va + vb + i_v) * sw.f, {}, {}};"),
    # Real data set: missing files of unused extra cuts are skipped (as the Fortran), of used cuts not
    ("a missing file of an unused extra cut rejected", "src/crplib/dataset.cpp",
     '      if (static_cast<int>(k) < used) fail(file, "molecule-surface cut file not found");',
     '      fail(file, "molecule-surface cut file not found");'),
    ("validation accepts a used cut without a file", "src/crplib/dataset.cpp",
     '      if (static_cast<int>(k) < required_molecule_cuts) throw DataError(what + ": molecule-surface cut file not found");\n',
     ""),
]

if sys.argv[1] == "list":
    print(len(MUTATIONS))
    sys.exit(0)
root, index = sys.argv[1], int(sys.argv[2])
name, path, old, new = MUTATIONS[index]
text = open(f"{root}/{path}").read()
if text.count(old) != 1:
    sys.exit(f"mutation {index} ({name}): pattern found {text.count(old)} times")
if "std::floor" in new and "#include <cmath>" not in text:
    text = text.replace("#include ", "#include <cmath>\n#include ", 1)
open(f"{root}/{path}", "w").write(text.replace(old, new))
print(name)

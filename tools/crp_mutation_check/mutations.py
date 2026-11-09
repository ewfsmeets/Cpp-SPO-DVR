# Mutations of the CRP library (docs/crp_library.md "Mutation check"): each must make a CrpLib test fail.
# Usage: python3 -I mutations.py <tree> <index>   (applies mutation <index> to <tree>; prints its name)
#        python3 -I mutations.py list
import sys

MUTATIONS = [
    ("bridge images swapped", "src/crplib/hexagonal_interaction.cpp",
     "  s[4] = hexagonal_bridge(bridge, t, phi + pis3);\n  s[5] = hexagonal_bridge(bridge, t, phi - pis3);",
     "  s[4] = hexagonal_bridge(bridge, t, phi - pis3);\n  s[5] = hexagonal_bridge(bridge, t, phi + pis3);"),
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
     "    s[2] = hexagonal_hollow(hollow, t, phi);  // hcp\n    const auto t2h = cuts.subspan(12, 6);",
     "    s[2] = hexagonal_hollow(hollow, t, phi);  // hcp\n    const auto t2h = cuts.subspan(18, 6);"),
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

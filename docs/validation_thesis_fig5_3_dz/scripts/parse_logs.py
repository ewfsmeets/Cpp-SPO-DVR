"""Parsers of the Fortran stdout (cleaned `script` log) and the C++ screen log.

Values are kept as the printed text (Fortran) so that the printed precision can be derived from
it; C++ values are parsed to float (full precision where the log prints %.16e).
"""
from __future__ import annotations

import re

from common import parse_fortran_real

NUM = r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][+-]?\d+)?"


# ============================================================================ Fortran stdout
def section(lines, title, stop_prefix="INFO FROM"):
    """Lines of the Fortran section whose banner contains `title` (until the next banner)."""
    start = None
    for i, line in enumerate(lines):
        if title in line and (i + 1 < len(lines)) and lines[i + 1].startswith("*" * 20):
            start = i + 2
            break
    if start is None:
        raise KeyError(title)
    out = []
    for line in lines[start:]:
        if line.startswith("*" * 40):
            break
        out.append(line)
    return out


def star_fields(lines, start_marker, end_marker):
    """'* Label:   value' lines between two markers -> [(label, value text)] in order (labels
    repeat in the Fortran input echo, e.g. the optical potentials)."""
    out = []
    active = False
    for line in lines:
        if start_marker in line:
            active = True
            continue
        if active and end_marker in line:
            break
        if not active:
            continue
        s = line.strip()
        if not s.startswith("* "):
            continue
        body = s[2:]
        # 'Label (unit):   value' ; some labels are cut at 58 characters without ':'
        m = re.match(r"^(.*?):\s{2,}(\S.*)$", body) or re.match(r"^(.*?):\s+(\S.*)$", body)
        if m:
            out.append((m.group(1).strip(), m.group(2).strip()))
        else:
            m = re.match(r"^(.*\S)\s+(" + NUM + r")$", body)
            if m:
                out.append((m.group(1).strip(), m.group(2)))
    return out


def equals_fields(lines):
    """'name = value' and 'z(nzchsp)= value' lines -> list of (name, value text) in order."""
    out = []
    for line in lines:
        m = re.match(r"^\s*(.*?)\s*=\s*(" + NUM + r")\s*$", line)
        if m:
            out.append((m.group(1).strip(), m.group(2)))
    return out


def internal_parameters(lines):
    return dict(star_fields(lines, "THE INPUT PARAMETERS SETUP INTERNALLY", "NOTE:"))


def input_read(lines):
    return star_fields(lines, "*** THE INPUT READ ***", "* Validating the input")


def grid_parameters(lines):
    out = {}
    active = False
    for line in lines:
        if "THE GRID PARAMETERS" in line:
            active = True
            continue
        if active and "CONTINUING SETUP II" in line:
            break
        if active:
            m = re.match(r"^\s*\*\s*(.*?)\s*:\s*(" + NUM + r")\s*$", line)
            if m:
                out[m.group(1)] = m.group(2)
    return out


def vasy_table(lines):
    sec = section(lines, "INFO FROM SETTING UP THE ROVIBARIONAL ASYMPTOTIC POTENTIAL")
    rows = []
    z = None
    for line in sec:
        m = re.search(r"potential for Z =\s*(" + NUM + ")", line)
        if m:
            z = m.group(1)
        t = line.split()
        if len(t) == 2 and all(re.fullmatch(NUM, x) for x in t):
            rows.append((t[0], t[1]))
    return z, rows


def rovib_states(lines):
    sec = section(lines, "INFO FROM SETTING UP THE ASYMPTOTIC ROVIBRATIONAL STATES")
    rows = []
    extra = {}
    for line in sec:
        t = line.split()
        if len(t) == 8 and t[0].isdigit() and t[1].isdigit():
            rows.append({"v": int(t[0]), "j": int(t[1]), "e": t[2], "erel": t[3], "psi1": t[4],
                         "psinr": t[5], "norm": t[6], "r": t[7]})
        else:
            m = re.match(r"^\s*(.*?)\s*=\s*(" + NUM + r")\s*$", line)
            if m:
                extra[m.group(1)] = m.group(2)
    return rows, extra


def averages_out(path):
    rows = []
    for line in open(path):
        t = line.split()
        if len(t) == 5 and all(re.fullmatch(NUM, x) for x in t):
            rows.append(t)
    return rows


def energy_blocks(lines):
    """Per-energy blocks of the on-the-fly analysis (Analyzor.f90 smatrix_combine)."""
    blocks = []
    sectn = None
    current_vjm = None
    for line in lines:
        if " total energy =" in line:
            blocks.append({"total_energy": line.split()[-1], "p_v": {}, "edif_v": {}, "erot_v": {},
                           "p_vjm": {}, "p_order": {}})
            sectn = None
            continue
        if not blocks:
            continue
        b = blocks[-1]
        if "trans energy in z =" in line:
            b["trans_energy_z"] = line.split()[-1]
            continue
        if "trans energy =" in line:
            b["trans_energy"] = line.split()[-1]
            continue
        if "total momentum =" in line:
            b["total_momentum"] = line.split()[-1]
            continue
        if "rotation-diffraction summed probabilities" in line:
            sectn = "v"
            continue
        if "diffraction summed probabilities" in line:
            sectn = "vjm"
            continue
        if "probabilities as function of diffraction order" in line:
            sectn = "order"
            continue
        if "sum of the probabilities =" in line:
            b["p_sum"] = line.split()[-1]
            sectn = None
            continue
        m = re.match(r"^\s*average (vibrational|rotational|diffraction)\s+energy \(ev\) =\s*(" + NUM + ")", line)
        if m:
            b["e_" + m.group(1)[:3]] = m.group(2)
            continue
        if "the reaction probability =" in line:
            b["p_reaction"] = line.split()[-1]
            continue
        t = line.split()
        if sectn == "v" and len(t) == 4 and t[0].isdigit():
            v = int(t[0])
            b["p_v"][v], b["edif_v"][v], b["erot_v"][v] = t[1], t[2], t[3]
        elif sectn == "vjm" and len(t) == 5 and t[0].isdigit():
            b["p_vjm"][(int(t[0]), int(t[2]), int(t[3]))] = t[4]
        elif sectn == "order":
            if line.strip().startswith("v, j, mj ="):
                v, j, mj = (int(x) for x in line.split("=")[1].split())
                current_vjm = (v, j, mj)
            elif len(t) == 2 and t[0].isdigit() and current_vjm is not None:
                b["p_order"][current_vjm + (int(t[0]),)] = t[1]
    return blocks


def fort25(path):
    """fort.25 (writeProbs): per energy a header with 'total energy', then per (v, j, mJ, order)
    rows 'n m ReS ImS |S|^2 phase'. v is printed as iv-1. Yields dicts."""
    energy_index = -1
    vjmo = None
    with open(path) as f:
        for line in f:
            if " total energy =" in line:
                energy_index += 1
                continue
            if "v j mj dif.order are:" in line:
                vjmo = tuple(int(x) for x in line.split(":")[1].split())
                continue
            t = line.split()
            if len(t) == 6 and vjmo is not None and re.fullmatch(r"-?\d+", t[0]) and re.fullmatch(r"-?\d+", t[1]):
                yield (energy_index, vjmo[0], vjmo[1], vjmo[2], vjmo[3], int(t[0]), int(t[1]), t[2], t[3], t[4], t[5])


# ============================================================================ C++ log
def cpp_section(lines, title):
    """Lines of a C++ section (banner of '=' lines around the title)."""
    for i, line in enumerate(lines):
        if line.strip() == title and i > 0 and lines[i - 1].startswith("====") and lines[i + 1].startswith("===="):
            out = []
            for l2 in lines[i + 2:]:
                if l2.startswith("=" * 40):
                    break
                out.append(l2)
            return out
    raise KeyError(title)


def cpp_dotted(lines):
    """'  Label ...... value' lines -> {label: value}."""
    out = {}
    for line in lines:
        m = re.match(r"^\s{2,}(\S.*?)\s\.{2,}\s(.*)$", line)
        if m:
            out[m.group(1).strip()] = m.group(2).strip()
    return out


def cpp_table(lines, ncols, first_int=True):
    rows = []
    for line in lines:
        t = line.split()
        if len(t) == ncols and (not first_int or re.fullmatch(r"-?\d+", t[0])):
            try:
                rows.append([float(x) if not re.fullmatch(r"-?\d+", x) else int(x) for x in t])
            except ValueError:
                pass
    return rows


def cpp_rovib(lines):
    sec = cpp_section(lines, "GAS-PHASE ROVIBRATIONAL STATES")
    rows = []
    for line in sec:
        t = line.split()
        if len(t) == 8 and t[0].isdigit() and t[1].isdigit():
            rows.append({"v": int(t[0]), "j": int(t[1]), "e": float(t[2]), "erel": float(t[3]),
                         "psi1": float(t[4]), "psinr": float(t[5]), "norm": float(t[6]), "r": float(t[7])})
    return rows, cpp_dotted(sec)


def cpp_surface_fixed(lines):
    for line in lines:
        m = re.search(r"analysis surface on the specular grid fixed at Z = (" + NUM + r") bohr \(t = (" + NUM + r") a\.u\.\)", line)
        if m:
            return float(m.group(1)), float(m.group(2))
    return None


def fortran_surface_fixed(lines):
    for line in lines:
        m = re.search(r"ZAnalysis on specular channel fixed at Z =\s*(" + NUM + ")", line)
        if m:
            return m.group(1)
    return None


def parse_number(text):
    return parse_fortran_real(text)

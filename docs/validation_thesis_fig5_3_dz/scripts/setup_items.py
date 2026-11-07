#!/usr/bin/env python3
"""Items 1-3 and 5-12: the set-up of the three runs, as a Fortran <-> C++ label map.

Every row: item, quantity, Fortran label and text, C++ label and value (CPU run; the GPU run's
set-up output is compared with the CPU run's as text), kind of comparison, difference,
tolerance, verdict. Kinds:
  int      integers (C++ 0-based indices converted: offset +1)
  text     exact text after normalisation (flags, names)
  printed  Fortran printed value vs C++ value; tolerance = 1.01 x half a unit of the last
           printed digit of the Fortran text (+ the same for the C++ text when the C++ side is
           printed with fewer than 16 digits)
  full     both sides double precision: identical, or rounding if |d| <= 1e-13 max(|x|, floor)
  n/a      no counterpart (reason in the note)
Writes data/items/setup.json and the Markdown fragment tables/setup_label_map.md."""
from __future__ import annotations

import json
import math
import re

import numpy as np

from common import (CHAEV, INPUT_DIR, RUN_DIRS, TABLES, FortranBk, NativeBk, clean_log_lines, load_item,
                    parse_fortran_real, printed_half_unit, save_item, WORK)
from parse_logs import (averages_out, cpp_dotted, cpp_rovib, cpp_section, energy_blocks, equals_fields,
                        grid_parameters, input_read, internal_parameters, rovib_states, section, star_fields,
                        vasy_table)

FULL_REL = 1e-13
ROWS = []


def cpp_text_half_unit(text: str) -> float:
    """Half unit of the last digit of a C++ printed value (0 if printed with >= 16 digits)."""
    t = text.strip()
    m = re.match(r"^[+-]?(\d*)\.(\d+)(?:e([+-]\d+))?$", t)
    if not m:
        return 0.0
    digits = len(m.group(2))
    exp = int(m.group(3)) if m.group(3) else 0
    if m.group(3) and digits >= 16:
        return 0.0
    return 0.5 * 10.0 ** (exp - digits)


def add(item, quantity, f_label, f_text, c_label, c_value, kind, offset=0, note="", c_text=None, floor=0.0,
        gpu_same=True, full_rel=FULL_REL):
    row = {"item": item, "quantity": quantity, "fortran_label": f_label, "fortran": f_text, "cpp_label": c_label,
           "cpp": c_value if c_text is None else c_text, "kind": kind, "note": note, "gpu_equal_to_cpu": gpu_same}
    if kind == "int":
        f = int(f_text)
        c = int(c_value) + offset
        row["delta"] = c - f
        row["verdict"] = "identical" if c == f else "FINDING"
    elif kind == "text":
        row["verdict"] = "identical" if str(f_text).strip() == str(c_value).strip() else "FINDING"
    elif kind == "printed":
        f = parse_fortran_real(f_text)
        c = float(c_value)
        tol = printed_half_unit(f_text) + (1.01 * cpp_text_half_unit(c_text) if c_text else 0.0)
        d = c - f
        row.update(delta=d, tolerance=tol, ratio=abs(d) / tol if tol else (0.0 if d == 0 else float("inf")))
        row["verdict"] = "printed" if abs(d) <= tol else "OUTSIDE"
    elif kind == "full":
        f = float(f_text)
        c = float(c_value)
        d = c - f
        scale = max(abs(f), floor)
        row.update(delta=d, relative=abs(d) / scale if scale else (0.0 if d == 0 else float("inf")),
                   ulps=abs(d) / np.spacing(abs(f)) if f != 0 else (0.0 if d == 0 else float("inf")))
        row["verdict"] = "identical" if d == 0 else ("rounding" if row["relative"] <= full_rel else "OUTSIDE")
    elif kind == "n/a":
        row["verdict"] = "n/a"
    else:
        raise ValueError(kind)
    ROWS.append(row)
    return row


def fl(x):
    return float(x)


def first_num(text):
    return float(re.match(r"^\s*([+-]?[\d.]+(?:e[+-]?\d+)?)", text).group(1))


def nums(text):
    return [float(x) for x in re.findall(r"[+-]?\d+\.?\d*(?:e[+-]?\d+)?", text)]


# ----------------------------------------------------------------------------- DVR.inp
def read_dvr_inp(path):
    out = {}
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = re.match(r'^(\w+)\s*:\s*(?:(INTEGER|REAL|LOGICAL|STRING)\s*:\s*)?(.*?);\s*$', line)
        if m:
            out[m.group(1)] = (m.group(2), m.group(3).strip().strip('"'))
    return out


def json_vs_dvr(js, dvr):
    """input.json <-> DVR.inp field by field (item 1a)."""
    rows = []
    for key, value in js.items():
        if key not in dvr:
            rows.append({"key": key, "json": value, "dvr": None, "equal": False})
            continue
        typ, text = dvr[key]
        if isinstance(value, bool):
            eq = (text.upper() == "TRUE") == value
        elif isinstance(value, (int, float)):
            eq = float(text.replace("D", "E")) == float(value)
        else:
            eq = text == value
        rows.append({"key": key, "json": value, "dvr": text, "type": typ, "equal": eq})
    extra = [k for k in dvr if k not in js]
    return rows, extra


ECHO_TO_KEY = [
    ("SystemID", "SystemID"), ("NewCalc", "NewCalc"), ("Z degree of freedom", "ZDoF"), ("R degree of freedom", "RDoF"),
    ("X degree of freedom", "XDoF"), ("Y degree of freedom", "YDoF"), ("Theta degree of freedom", "ThetaDoF"),
    ("Phi degree of freedom", "PhiDoF"), ("Mass of atom 1 (amu)", "Mass1"), ("Mass of atom 2 (amu)", "Mass2"),
    ("Homo-nuclear molecule", "HomoNuc"), ("Number of points on the scattering grid in Z", "nZ"),
    ("Step size of the scattering Z (bohr)", "dZ"), ("Starting point of the scattering in Z (bohr)", "ZStart"),
    ("Number of points on the specular grid in Z", "nZSpecular"), ("Number of points on the scattering grid in r", "nR"),
    ("Step size of the scattering r (bohr)", "dR"), ("Starting point of the scattering in r (bohr)", "RStart"),
    ("Number of points on the scattering grid in X", "nX"), ("Number of points on the scattering grid in Y", "nY"),
    ("Maximum J value in grid representation.", "JMax"), ("Maximum Mj value in grid representation.", "MJMax"),
    ("Calculation starts at (aut)", "tStart"), ("Calculation ends at (aut)", "tEnd"),
    ("Time step used in the propagation (aut)", "dtStep"), ("Number of time steps taken", "nStep"),
    ("Wave function analysed every (aut)", "dtAna"), ("Averages printed every (aut)", "dtAve"),
    ("Initially using a dynamical specular analysis surface", "DynamZSpAna"), ("Wave function analysed at (bohr)", "Zana"),
    ("Number of energies that are on-the-fly analysed", "NrOfEs"), ("Lowest energy of the initial wave packet", "EMinWP"),
    ("Highest energy of the initial wave packet", "EMaxWP"), ("Average position of the initial wave packet", "Z0"),
    ("Initial rotational quantum number J", "JIn"), ("Initial rotational quantum number Mj", "MJIn"),
    ("Initial rotational quantum number V", "VIn"), ("Using symmetry in setting up the potential", "UseSymmetry"),
    ("Will try to find an optimal number of phis", "nPhisOptimal"),
    ("Optical potential in Z start at (bohr)#1", "ZOptStart"), ("Optical potential in Z ends at (bohr)#1", "ZOptEnd"),
    ("Optical potential in Z absorbs optimally for (eV)", "ZOptEnergy"), ("Optical potential in r start at (bohr)", "ROptStart"),
    ("Optical potential in r ends at (bohr)", "ROptEnd"), ("Optical potential in r absorbs optimally for (eV)", "ROptEnergy"),
    ("Optical potential in Z start at (bohr)#2", "ZSpOptStart"), ("Optical potential in Z ends at (bohr)#2", "ZSpOptEnd"),
    ("Optical potential in Z specular absorbs optimally for (eV", "ZSpOptEnergy"),
    ("Initial energy parallel to the surface (eV)", "EParallel"), ("Azimuthal angle of the incident molecule (degrees)", "PhiParallel"),
    ("S-matrix elements written to file 3", "writeSMatrix"), ("Rot-vib-diff resolved probs written to file 25", "writeProbs"),
    ("Rot-diff resolved probs written to file 10", "writeRDProbs"),
    ("Tolerance for transferring the wavefunction (nproj = 2)", "TolSp"),
]


def echo_vs_input(echo, js):
    """Fortran 'THE INPUT READ' echo vs the input values (shows the Fortran read the intended input)."""
    seen = {}
    labelled = []
    for label, text in echo:
        seen[label] = seen.get(label, 0) + 1
        labelled.append((f"{label}#{seen[label]}" if label.startswith("Optical potential in Z") and "specular" not in label
                         and "absorbs" not in label else label, text))
    d = dict(labelled)
    rows = []
    for label, key in ECHO_TO_KEY:
        text = d[label]
        value = js[key]
        if isinstance(value, bool):
            ok = (text == "T") == value
            tol = None
        elif isinstance(value, str):
            ok = text == value
            tol = None
        elif isinstance(value, int) and re.fullmatch(r"-?\d+", text):
            ok = int(text) == value
            tol = None
        else:
            tol = printed_half_unit(text)
            ok = abs(parse_fortran_real(text) - float(value)) <= tol
        rows.append({"label": label.split("#")[0], "key": key, "fortran": text, "input": value, "ok": ok, "tolerance": tol})
    return rows


# ----------------------------------------------------------------------------- main
def main():
    FL = clean_log_lines("fortran")
    CL = clean_log_lines("cpu")
    GL = clean_log_lines("gpu")
    js = json.loads((INPUT_DIR / "input.json").read_text())
    dvr = read_dvr_inp(INPUT_DIR / "DVR.inp")
    out = {}

    # ---------------- item 1: input
    jrows, extra = json_vs_dvr(js, dvr)
    out["json_vs_dvr"] = {"rows": jrows, "n": len(jrows), "n_equal": sum(r["equal"] for r in jrows),
                          "unequal": [r["key"] for r in jrows if not r["equal"]], "dvr_only": extra}
    erows = echo_vs_input(input_read(FL), js)
    out["fortran_echo_vs_input"] = {"rows": erows, "n": len(erows), "n_ok": sum(r["ok"] for r in erows),
                                    "not_ok": [r["key"] for r in erows if not r["ok"]]}
    cin = cpp_dotted(cpp_section(CL, "INPUT"))
    cgr = cpp_dotted(cpp_section(CL, "GRIDS"))
    # C++ echo vs input.json (provenance: the log is the only record of the input the C++ runs read)
    checks = []

    def chk(label, ok):
        checks.append({"check": label, "ok": bool(ok)})

    chk("Title", cin["Title"] == js["Title"])
    chk("System", cin["System"] == js["SystemID"])
    chk("Run mode", cin["Run mode"] == ("new calculation" if js["NewCalc"] else "?"))
    chk("Masses", nums(cin["Masses"]) == [js["Mass1"], js["Mass2"]])
    chk("Homonuclear", cin["Homonuclear"] == ("yes" if js["HomoNuc"] else "no"))
    chk("Time step", first_num(cin["Time step"]) == js["dtStep"])
    chk("Steps (total / per analysis)", nums(cin["Steps (total / per analysis)"]) == [js["nStep"], js["dtAna"] / js["dtStep"]])
    chk("95% confidence interval", nums(cin["95% confidence interval"]) == [js["EMinWP"], js["EMaxWP"]])
    zs = nums(cgr["Z (scattering)"])
    chk("Z (scattering): points, start, step", [zs[0], zs[1], zs[3]] == [js["nZ"], js["ZStart"], js["dZ"]])
    zsp = nums(cgr["Z (specular)"])
    chk("Z (specular): points, start, step", [zsp[0], zsp[1], zsp[3]] == [js["nZSpecular"], js["ZStart"], js["dZ"]])
    rr = nums(cgr["r"])
    chk("r: points, start, step", [rr[0], rr[1], rr[3]] == [js["nR"], js["RStart"], js["dR"]])
    chk("x, y points", [nums(cgr["x"])[0], nums(cgr["y"])[0]] == [js["nX"], js["nY"]])
    chk("J range", nums(cgr["J range (step)"])[1] == js["JMax"])
    chk("MJMax (basis)", int(cgr["MJMax (basis)"]) == js["MJMax"])
    chk("Analysis at Z", first_num(cgr["Analysis at Z (index)"]) == js["Zana"])
    for k in ("cpu", "gpu"):
        nb = NativeBk(RUN_DIRS[k] / "bk_coefficients.bin")
        p = nb.parameters
        keymap = {"ZStart": "ZStart", "dZ": "dZ", "nZ": "nZ", "nZSpecular": "nZSpecular", "RStart": "RStart", "dR": "dR",
                  "nR": "nR", "nX": "nX", "nY": "nY", "JMax": "JMax", "MJMax": "MJMax", "Mass1": "Mass1", "Mass2": "Mass2",
                  "dtStep": "dtStep", "dtAna": "dtAna", "Zana": "Zana", "Z0": "Z0", "EMinWP": "EMinWP", "EMaxWP": "EMaxWP",
                  "EParallel": "EParallel", "PhiParallel": "PhiParallel", "VIn": "VIn", "JIn": "JIn", "MJIn": "MJIn",
                  "ZOptStart": "ZOptStart", "ZOptEnd": "ZOptEnd", "ROptStart": "ROptStart", "ROptEnd": "ROptEnd",
                  "ZSpOptStart": "ZSpOptStart", "ZSpOptEnd": "ZSpOptEnd"}
        bad = [kk for kk, jk in keymap.items() if float(p[kk]) != float(js[jk])]
        chk(f"{k}: restart parameters in bk_coefficients.bin = input.json (30 keys)", not bad)
    out["cpp_echo_vs_input"] = {"checks": checks, "n": len(checks), "n_ok": sum(c["ok"] for c in checks)}
    out["title_lengths"] = {"json": len(js["Title"]), "dvr": len(dvr["Title"][1])}

    ip = internal_parameters(FL)
    add(1, "lattice constant x", "The surface lattice constant in X (bohr)", ip["The surface lattice constant in X (bohr)"],
        "Lattice constants (x, y)", nums(cin["Lattice constants (x, y)"])[0], "printed", c_text="4.840000")
    add(1, "lattice constant y", "The surface lattice constant in Y (bohr)", ip["The surface lattice constant in Y (bohr)"],
        "Lattice constants (x, y)", nums(cin["Lattice constants (x, y)"])[1], "printed", c_text="4.840000")
    add(1, "Z asymptotic", "The asymptotic value in Z (bohr)", ip["The asymptotic value in Z (bohr)"],
        "Z asymptotic", first_num(cin["Z asymptotic"]), "printed", c_text="15.000000")
    add(1, "skew angle", "The unit cell skewing angle (degrees)", ip["The unit cell skewing angle (degrees)"],
        "Skewing angle", first_num(cin["Skewing angle"]), "printed", c_text="60.000000")
    add(1, "point group", "The unit cell point group", ip["The unit cell point group"], "Point group", cin["Point group"], "text")
    add(1, "r cut-off", "The cutoff value in r (bohr)", ip["The cutoff value in r (bohr)"], "-", None, "n/a",
        note="UseRCutOff is not ported (forced off by the Fortran ProcessInput; fortran_deviations.md)")
    echo = dict(input_read(FL))
    add(1, "projection scheme", "The projection operator scheme used", echo["The projection operator scheme used"],
        "ProjOp (input; only 2 is accepted)", "using two grids all the way" if js["ProjOp"] == 2 else "?", "text")
    add(1, "switching function", "Using a switching function for the potential in Z", echo["Using a switching function for the potential in Z"],
        "Switching to gas phase", "F" if cin["Switching to gas phase"] == "no" else "T", "text")
    add(1, "only/also fine analysis", "Only / Also performing a fine analysis",
        echo["Only performing a fine analysis"] + echo["Also performing a fine analysis"], "Run mode", "FF" if cin["Run mode"] == "new calculation" else "?", "text")
    add(1, "analysed J, v, diffraction order", "Maximum J/V/diffraction order analysed", echo["Maximum J analysed"],
        "(JAna, VAna, MDifAna absent from input.json)", "determined by the program (see info Analysis)"
        if not any(k in js for k in ("JAna", "VAna", "MDifAna")) else "?", "text")
    add(1, "potential limits, mixing", "Imposing a minimum / maximum value; mixing alfa",
        echo["Imposing a minimum value on the potential"] + echo["Imposing a maximum value on the potential"]
        + echo["Potential mixing alfa value"], "(VMin, VMax absent; PESMixingAlfa rejected)", None, "n/a",
        note="C++ has no mixing (NotSupportedError); no limits in either input")
    add(1, "mJ-independent grid", "A general MjIn independant grid will be used", "yes" if any("MjIn independant grid" in l for l in FL) else "no",
        "-", None, "n/a", note="Fortran setupdvr.f90 message; the C++ grid is always MJIn independent")

    # ---------------- item 2: grids
    gp = grid_parameters(FL)
    pot = load_item("item04")
    add(2, "number of phi points", "Number of phi's", gp["Number of phi's"], "Phi points", cgr["Phi points"], "int")
    add(2, "number of theta points", "Number of theta's", gp["Number of theta's"], "Theta points", cgr["Theta points"], "int")
    jr = nums(cgr["J range (step)"])
    add(2, "J max in basis", "JMax in basis set", gp["JMax in basis set"], "J range (step)", jr[1], "int")
    add(2, "J min in basis", "JMin in basis set", gp["JMin in basis set"], "J range (step)", jr[0], "int")
    add(2, "J step in basis", "Delta J in basis set", gp["Delta J in basis set"], "J range (step)", jr[2], "int")
    gtxt = (WORK / "potential" / "build-fast" / "grids.txt").read_text().split("\n")
    phis = []
    for i, line in enumerate(gtxt):
        if line.startswith("phi "):
            n = int(line.split()[1])
            phis = [float(x) for x in gtxt[i + 1:i + 1 + n]]
    add(2, "phi step", "Delta value in phi coordinate", gp["Delta value in phi coordinate"], "phi[1] - phi[0] (helper grids.txt)",
        phis[1] - phis[0], "full", note="Fortran prints 17 significant digits")
    npts = [l for l in FL if "Number of points in (x,y) for which the potential is calculated" in l][0].split("=")[1].strip()
    add(2, "(x, y) points with V computed", "Number of points in (x,y) for which the potential is calculated", npts,
        "(x, y) points computed / total", nums(cpp_dotted(cpp_section(CL, "SURFACE SYMMETRY OF THE POTENTIAL"))["(x, y) points computed / total"])[0], "int")
    add(2, "(x, y) points in .pes vs C++ representatives", ".pes n_xy_computed", str(pot["dims"]["n_xy_computed"]),
        "C++ representatives", pot["cpp_n_representatives"], "int",
        note=f"different fundamental regions: {pot['n_common_representatives']} of 19 points in common")
    for k, lab in (("z1", "Z scattering, first"), ("zn", "Z(nzpot), last point of the potential"), ("r1", "r first"),
                   ("rn", "r last"), ("x1", "x first"), ("xn", "x last"), ("y1", "y first"), ("yn", "y last")):
        e = pot["grid_ends"][k]
        add(2, lab, f".pes header {k}", repr(e["fortran"]), "helper grids.txt", e["cpp"], "full")
    add(2, "theta nodes (15)", ".pes header thetas", "15 values", "helper grids.txt / log theta table",
        None, "n/a", note=f"see full comparison: max |d| {pot['theta']['max_abs_diff']:.2e} rad, "
        f"{pot['theta']['n_identical']}/15 identical, max {pot['theta']['max_ulps']:.0f} ulp")
    out["theta"] = pot["theta"]
    out["phi"] = pot["phi"]
    # theta nodes vs C++ log table (log printed with %.16e = 17 digits)
    th_log = [r[1] for r in [[float(x) for x in l.split()] for l in cpp_section(CL, "GRIDS") if re.match(r"^\s+\d+\s+\S+e[+-]\d+\s+\S+e[+-]\d+$", l)]]
    out["theta_log_equals_helper"] = th_log == [float(x) for x in gtxt[gtxt.index("theta 15") + 1:gtxt.index("theta 15") + 16]]

    # ---------------- item 3: MISC SETUP
    misc = equals_fields(section(FL, "INFO FROM MISCELLANEOUS SETUP"))
    m = {}
    for k, v in misc:
        m.setdefault(k, v)
    dkz = nums(cgr["dkz (scattering / specular)"])
    add(3, "dk_Z scattering", "grid interval size in kz dir. (scat. grid)", m["grid interval size in kz dir. (scat. grid)"],
        "dkz (scattering / specular)", dkz[0], "printed")
    add(3, "dk_Z specular", "grid interval size in kz dir. (spec. grid)", m["grid interval size in kz dir. (spec. grid)"],
        "dkz (scattering / specular)", dkz[1], "printed")
    add(3, "initial k_Z", "initial momentum in z direction", m["initial momentum in z direction"], "Initial kz", first_num(cin["Initial kz"]), "printed")
    add(3, "dk_r", "grid interval size in kr dir.", m["grid interval size in kr dir."], "dkr", first_num(cgr["dkr"]), "printed")
    dkxy = nums(cgr["dkx, dky"])
    add(3, "dk_x", "grid interval size in kx dir.", m["grid interval size in kx dir."], "dkx, dky", dkxy[0], "printed")
    add(3, "dk_y", "grid interval size in ky dir.", m["grid interval size in ky dir."], "dkx, dky", dkxy[1], "printed")
    rbk = cpp_dotted(cpp_section(CL, "ROTATIONAL BASIS AND INITIAL WAVE PACKET"))
    add(3, "initial k_x", "initial momentum in reciprocal x direction", m["initial momentum in reciprocal x direction"],
        "Parallel kinetic energy (hartree)", first_num(rbk["Parallel kinetic energy"]), "printed",
        note="C++ prints the parallel kinetic energy (0 for k_x = k_y = 0)")
    add(3, "initial k_y", "initial momentum in reciprocal y direction", m["initial momentum in reciprocal y direction"],
        "Parallel kinetic energy (hartree)", first_num(rbk["Parallel kinetic energy"]), "printed")
    add(3, "nmoder", "For backward compatibility: nmoder", m["For backward compatibility: nmoder"], "-", None, "n/a",
        note="legacy diagnostic, no C++ counterpart")
    add(3, "analysis line Z", "asymptotic analysis at z", m["asymptotic analysis at z"], "Analysis at Z (index)",
        first_num(cgr["Analysis at Z (index)"]), "printed", c_text="15.000000")
    add(3, "flux analysis", "No flux analysis performed.", "no" if any("No flux analysis performed" in l for l in FL) else "yes",
        "(no FluxAnalysis key)", "no" if "FluxAnalysis" not in js else "yes", "text")
    add(3, "nzpot", "nzpot", m["nzpot"], "Z points in interaction region", cgr["Z points in interaction region"], "int")
    zo = [int(x) for x in nums(cgr["Optical potential Z (indices)"])]
    zso = [int(x) for x in nums(cgr["Optical potential Z specular (indices)"])]
    ro = [int(x) for x in nums(cgr["Optical potential r (indices)"])]
    gz = {}
    for i, line in enumerate(gtxt):
        t = line.split()
        if len(t) == 2 and t[0] in ("z_scattering", "z_specular", "r"):
            gz[t[0]] = [float(x) for x in gtxt[i + 1:i + 1 + int(t[1])]]
    for f_idx, c_idx, f_z, grid, lab in (("nzchsp", zo[0], "z(nzchsp)", "z_scattering", "Optical potential Z (indices)"),
                                          ("nzchmx", zo[1], "z(nzchmx)", "z_scattering", "Optical potential Z (indices)"),
                                          ("nschsp", zso[0], "z(nschsp)", "z_specular", "Optical potential Z specular (indices)"),
                                          ("nschmx", zso[1], "z(nschmx)", "z_specular", "Optical potential Z specular (indices)"),
                                          ("nrchsp", ro[0], "r(nrchsp)", "r", "Optical potential r (indices)"),
                                          ("nrchmx", ro[1], "r(nrchmx)", "r", "Optical potential r (indices)")):
        add(3, f_idx, f_idx, m[f_idx], lab + " (0-based)", c_idx, "int", offset=1)
        add(3, f_z, f_z, m[f_z], f"{grid}[{c_idx}] (helper grids.txt)", gz[grid][c_idx], "printed")
    add(3, "max kinetic energy in r", "max kinetic energy in r (ev)", m["max kinetic energy in r (ev)"],
        "Max. kinetic energy in r", first_num(cgr["Max. kinetic energy in r"]), "printed", c_text="3.654734")
    add(3, "max analysable total energy", "the max total energy ... (eV)", m["can be obtained (eV)"],
        "Max. analysis energy", first_num(cin["Max. analysis energy"]), "printed", c_text="4.274357")
    ci = nums(cin["95% confidence interval"])
    add(3, "ectmin", "ectmin (ev)", m["ectmin (ev)"], "95% confidence interval", ci[0], "printed", c_text="0.400000")
    add(3, "ectmax", "ectmax (ev)", m["ectmax (ev)"], "95% confidence interval", ci[1], "printed", c_text="0.950000")

    # ---------------- item 5: V_asy
    zv, vrows = vasy_table(FL)
    cv = [l.split() for l in (RUN_DIRS["cpu"] / "asymptotic_potential.dat").read_text().splitlines() if not l.startswith("#")]
    gpu_same = (RUN_DIRS["cpu"] / "asymptotic_potential.dat").read_bytes() == (RUN_DIRS["gpu"] / "asymptotic_potential.dat").read_bytes()
    dv, dr = [], []
    for (fr, fv), c in zip(vrows, cv):
        dr.append(float(c[1]) - parse_fortran_real(fr))
        dv.append(abs(float(c[3]) - parse_fortran_real(fv)) / max(abs(parse_fortran_real(fv)), 1.0))
    i = int(np.argmax(dv))
    out["vasy"] = {"n": len(dv), "z": zv, "max_scaled": max(dv), "argmax_ir": i, "argmax_r": float(cv[i][1]),
                   "median_scaled": float(np.median(dv)), "n_identical": int(sum(x == 0 for x in dv)),
                   "r_identical": all(x == 0 for x in dr), "gpu_file_identical": gpu_same, "criterion": 1e-14}
    out["vasy"]["verdict"] = "identical" if max(dv) == 0 else ("rounding" if max(dv) <= 1e-14 else "OUTSIDE")

    # ---------------- item 6: rovibrational states
    frows, fextra = rovib_states(FL)
    crows, cextra = cpp_rovib(CL)
    gas = [l.split() for l in (RUN_DIRS["cpu"] / "gas_phase_states.dat").read_text().splitlines() if not l.startswith("#")]
    psi_full = {(int(t[0]), int(t[1])): (float(t[3]), float(t[-1])) for t in gas}
    e_full = {(int(t[0]), int(t[1])): float(t[2]) for t in gas}
    worst = {}
    sign_flips = []
    n_cmp = 0
    n_out = []
    for f, c in zip(frows, crows):
        assert (f["v"], f["j"]) == (c["v"], c["j"])
        p1, pn = psi_full[(c["v"], c["j"])]
        pairs = {"e": (f["e"], c["e"]), "erel": (f["erel"], c["erel"]), "psi1": (f["psi1"], p1), "psinr": (f["psinr"], pn),
                 "norm": (f["norm"], c["norm"]), "r": (f["r"], c["r"])}
        fp1, fpn = parse_fortran_real(f["psi1"]), parse_fortran_real(f["psinr"])
        if fp1 * p1 < 0 and fpn * pn < 0:
            sign_flips.append((c["v"], c["j"]))
            pairs["psi1"] = (f["psi1"], -p1)
            pairs["psinr"] = (f["psinr"], -pn)
        for q, (ft, cval) in pairs.items():
            tol = printed_half_unit(ft) + (1.01 * 0.5e-12 if q in ("norm", "r") else 0.0)
            ratio = abs(cval - parse_fortran_real(ft)) / tol
            n_cmp += 1
            if ratio > 1:
                n_out.append({"v": c["v"], "j": c["j"], "q": q, "fortran": ft, "cpp": cval, "ratio": ratio})
            if ratio > worst.get(q, {"ratio": -1})["ratio"]:
                worst[q] = {"ratio": ratio, "v": c["v"], "j": c["j"], "fortran": ft, "cpp": cval}
    out["rovib"] = {"n_states": len(frows), "n_values": n_cmp, "worst_ratio_by_quantity": worst,
                    "outside": n_out, "sign_flipped_states": sign_flips}
    ce = {"energy initial v,j state": "energy", "potential energy": "potential energy", "rotational energy": "rotational energy",
          "vibrational energy": "vibrational energy", "rovibrational energy": "rovibrational kinetic energy", "expectation value r": "<r>"}
    for fk, ck in ce.items():
        add(6, "initial state: " + fk, fk, fextra[fk], ck, first_num(cextra[ck]), "printed")
    fb = FortranBk(RUN_DIRS["fortran"] / "fort.21")
    add(6, "initial rovibrational energy (hartree)", "fort.21 record 2: eviroi", repr(fb.e_initial),
        "gas_phase_states.dat E(0,0)", e_full[(0, 0)], "full", full_rel=1e-12,
        note="criterion 1e-12 relative (gate 2: LAPACK vs Jacobi)")

    # ---------------- items 7-9: bases
    dif = section(FL, "INFO FROM SETTING UP THE DIFFRACTION BASIS")
    rot = section(FL, "INFO FROM SETTING UP THE ROTATION BASIS")
    rd = section(FL, "INFO FROM SETTING UP THE ROTATION-DIFFRACTION BASIS")
    ana = cpp_dotted(cpp_section(CL, "ASYMPTOTIC ANALYSIS"))
    cpp_warn = [l.strip() for l in cpp_section(CL, "ASYMPTOTIC ANALYSIS") if l.strip().startswith("Warning")]
    f_warn = [l.strip() for l in dif if l.strip().startswith("WARNING")]
    add(7, "warning: basis too small", f_warn[0] if f_warn else "-", "yes" if f_warn else "no",
        cpp_warn[0] if cpp_warn else "-", "yes" if cpp_warn else "no", "text")

    def grab(lines, pattern):
        for l in lines:
            mm = re.search(pattern, l)
            if mm:
                return mm.group(1)
        raise KeyError(pattern)

    ndif = nums(ana["Diffraction channels analysed (ndifa)"])
    add(7, "diffraction channels in basis", "ndftow", grab(dif, r"(\d+) diffraction channels in basis"),
        "Diffraction channels analysed (ndifa): '144 of N'", ndif[1], "int")
    add(7, "diffraction channels analysed", "ndifa", grab(dif, r"ndifa =\s*(\d+)"), "Diffraction channels analysed (ndifa)", ndif[0], "int")
    add(7, "highest diffraction order", "maxord", grab(dif, r"maxord =\s*(\d+)"), "Highest diffraction order analysed (mdifa)",
        ana["Highest diffraction order analysed (mdifa)"], "int")
    nb = NativeBk(RUN_DIRS["cpu"] / "bk_coefficients.bin")
    first_rot = nb.labels["J"] * 1000 + nb.labels["mJ"]
    sel = first_rot == first_rot[0]
    orders = nb.labels["order"][sel]
    counts = np.bincount(orders)
    add(7, "largest number of channels of one order", "maxnum", grab(dif, r"maxnum =\s*(\d+)"),
        "derived: max count per order of the 144 channels (bk_coefficients.bin labels)", int(counts.max()), "int")
    add(8, "rotation channels in basis", "nrotw", grab(rot, r"(\d+) rotation channels in basis"),
        "Rotational basis functions", rbk["Rotational basis functions"], "int")
    add(8, "rotational states analysed", "nana", grab(rot, r"nana =\s*(\d+)"), "Rotational states analysed (nana)",
        ana["Rotational states analysed (nana)"], "int")
    add(8, "initial rotational state", "ichn", grab(rot, r"ichn =\s*(\d+)"), "Initial channel (index: J, mJ) (0-based)",
        rbk["Initial channel (index: J, mJ)"].split(":")[0], "int", offset=1)
    add(8, "initial rotational state (fort.21)", "fort.21 record 3: ichn", str(fb.ichn), "Initial channel (0-based)",
        rbk["Initial channel (index: J, mJ)"].split(":")[0], "int", offset=1)
    add(9, "rotation-diffraction channels in basis", "nstatw", grab(rd, r"(\d+) rotation-diffraction channels in basis"),
        "derived: rotational basis x diffraction channels", int(rbk["Rotational basis functions"]) * int(ndif[1]), "int")
    add(9, "rotation-diffraction states analysed", "nstata", grab(rd, r"nstata =\s*(\d+)"), "Rotation-diffraction states (nstata)",
        ana["Rotation-diffraction states (nstata)"], "int")
    add(9, "initial rot-dif channel", "inchna", grab(rd, r"inchna =\s*(\d+)"), "Initial channel (index) (0-based)",
        ana["Initial channel (index)"].split(":")[0], "int", offset=1)
    # labels of all 9504 states and channel energies at full precision
    lab_eq = {k: bool(np.array_equal(fb.labels[k], nb.labels[k])) for k in ("n", "m", "J", "mJ")}
    nbg = NativeBk(RUN_DIRS["gpu"] / "bk_coefficients.bin")
    lab_gpu = all(np.array_equal(nb.labels[k], nbg.labels[k]) for k in nb.labels)
    de = nb.energies.reshape(-1) - fb.enymo
    rel = np.abs(de) / np.abs(fb.enymo)
    out["state_labels"] = {"n_states": int(nb.nstata), "nviba": int(nb.nviba), "equal": lab_eq,
                           "gpu_equal_to_cpu": lab_gpu,
                           "rotation_pointer_consistent": bool(np.array_equal(
                               np.unique(fb.labels["PntToRot"]).size, int(grab(rot, r"nana =\s*(\d+)"))))}
    k = int(np.argmax(np.abs(de)))
    out["channel_energies"] = {"n": int(de.size), "n_identical": int(np.sum(de == 0)), "max_abs_hartree": float(np.max(np.abs(de))),
                               "max_relative": float(rel.max()), "median_abs": float(np.median(np.abs(de))),
                               "argmax": {"v": int(k // nb.nstata), "is": int(k % nb.nstata),
                                          "J": int(nb.labels["J"][k % nb.nstata]), "n": int(nb.labels["n"][k % nb.nstata]),
                                          "m": int(nb.labels["m"][k % nb.nstata]), "E_hartree": float(fb.enymo[k])},
                               "gpu_identical_to_cpu": bool(np.array_equal(nb.energies, nbg.energies))}
    # open channels at each analysis energy are compared in items 14-16 and 19

    # ---------------- item 10: optical potentials
    op = star_fields(FL, "INFO FROM SETTING UP THE OPTICAL POTENTIALS", "INFO FROM SETTING UP THE GRID ENERGIES")
    caps = nums(cin["CAP strengths (Z, r, Zsp)"])
    for (lab, text), c, ctext, key in zip(op, caps, cin["CAP strengths (Z, r, Zsp)"].replace(" hartree", "").split(", "),
                                         ("ZOptFactor (effective)", "ROptFactor (effective)", "ZSpOptFactor (effective)")):
        add(10, "CAP coefficient: " + lab.replace("Optical potential in ", "").replace(" absorbs with coefficient", ""),
            lab, text, "CAP strengths (log, %.8e)", c, "printed", c_text=ctext)
        add(10, "CAP coefficient (full precision C++)", lab, text, f"bk_coefficients.bin: {key}", nb.parameters[key], "printed")

    # ---------------- item 11: initial wave packet
    wf = section(FL, "INFO FROM SETTING UP THE WAVE FUNCTION")
    wfx = {}
    for l in wf:
        mm = re.match(r"^\s*(?:\*|-)\s*(.*?)\s*(?:=|\s{3,})\s*(" + r"[+-]?\d*\.?\d+(?:E[+-]\d+)?" + r")\s*$", l)
        if mm:
            wfx[mm.group(1)] = mm.group(2)
    add(11, "norm of the initial packet", "Norm initial wave packet (Z,r)", wfx["Norm initial wave packet (Z,r)"],
        "Norm of the initial wave packet", first_num(rbk["Norm of the initial wave packet"]), "printed")
    add(11, "|psi(k)| at k = 0", "at k = 0", wfx["at k = 0"], "|psi(k)| at k = 0", first_num(rbk["|psi(k)| at k = 0"]), "printed")
    add(11, "|psi(k)| at the grid edge", "at the end/beginning of the grid", wfx["at the end/beginning of the grid"],
        "|psi(k)| at the grid edge", first_num(rbk["|psi(k)| at the grid edge"]), "printed")
    add(11, "|psi(k)| maximum", "maximum value", wfx["maximum value"], "|psi(k)| maximum", first_num(rbk["|psi(k)| maximum"]), "printed")
    lo_hi = nums(rbk["Allowed region (indices)"])
    add(11, "iLow", "iLow", grab(wf, r"iLow =\s*(\d+)"), "Allowed region (indices) (0-based)", lo_hi[0], "int", offset=1)
    add(11, "iHigh", "iHigh", grab(wf, r"iHigh =\s*(\d+)"), "Allowed region (indices) (0-based)", lo_hi[1], "int", offset=1)
    add(11, "norm between iLow and iHigh", "Norm of the initial wave packet between iLow and iHigh",
        grab(wf, r"between iLow and iHigh\s*=\s*(\S+)"), "Norm inside the allowed region",
        first_num(rbk["Norm inside the allowed region"]), "printed")
    add(11, "initial k_Z (fort.21, full)", "fort.21 record 2: rkz", repr(fb.kz), "Initial kz (log %.16e)", first_num(cin["Initial kz"]), "full")
    add(11, "packet alpha (fort.21, full)", "fort.21 record 2: alfwv", repr(fb.alpha), "Wave packet alpha (log %.16e)",
        first_num(cin["Wave packet alpha"]), "full")
    add(11, "translational energy (fort.21, full)", "fort.21 record 3: enrgy", repr(fb.e_translational),
        "Translational energy (log %.16e)", first_num(cin["Translational energy"]), "full")
    add(11, "Z0 (fort.21)", "fort.21 record 2: z0", repr(fb.z0), "input Z0", js["Z0"], "full")
    add(11, "analysis line (fort.21)", "fort.21 record 2: ZAnalysis", repr(fb.z_analysis), "Analysis at Z", first_num(cgr["Analysis at Z (index)"]), "full")
    add(11, "analysis time step (fort.21)", "fort.21 record 2: deltat", repr(fb.deltat), "dtAna (bk header)", nb.parameters["dtAna"], "full")

    # ---------------- item 12: analysis set-up
    an = section(FL, "INFO FROM SETTING UP THE ANALYSIS")
    add(12, "tolerance for fixing the specular analysis line", "Tolerance for fixing the dynamical analysis surface",
        grab(an, r"=\s*(\S+E[+-]\d+)"), "tolerance for fixing it", first_num(rbk["tolerance for fixing it"]), "printed")
    add(12, "jana", "Maximum j in analysis: jana", grab(an, r"jana =\s*(\d+)"), "Highest j analysed (jana)", ana["Highest j analysed (jana)"], "int")
    add(12, "nviba", "Maximum v in analysis: nviba", grab(an, r"nviba =\s*(\d+)"), "Vibrational levels analysed (nviba)",
        ana["Vibrational levels analysed (nviba)"], "int")
    add(12, "mdifa", "Maximum diffraction order in analysis: mdifa", grab(an, r"mdifa =\s*(\d+)"),
        "Highest diffraction order analysed (mdifa)", ana["Highest diffraction order analysed (mdifa)"], "int")
    add(12, "number of analysis energies", "per-energy blocks in stdout", str(len(energy_blocks(FL))),
        "Analysis energies: number", nums(ana["number, step"])[0], "int")
    add(12, "specular analysis line start", "iHigh (Z = zsp[iHigh-1])", grab(wf, r"iHigh =\s*(\d+)"),
        "Analysis surface on the specular grid: starts at Z", gz["z_specular"].index(first_num(rbk["Analysis surface on the specular grid"].split("=")[1])),
        "int", offset=1)
    # energy grid of the analysis: all 111 energies
    blocks = energy_blocks(FL)
    summ = [l.split("\t") for l in (RUN_DIRS["cpu"] / "scattering_summary.tsv").read_text().splitlines() if not l.startswith("#")]
    head = summ[0]
    cols = {h: i for i, h in enumerate(head)}
    egrid = {"E_coll": [], "E_total": [], "E_z": [], "momentum": []}
    mtot = first_num(cin["Total mass"])
    for b, rr_ in zip(blocks, summ[1:]):
        for fkey, ckey in (("trans_energy", "E_coll"), ("total_energy", "E_total"), ("trans_energy_z", "E_z")):
            ft = b[fkey]
            c = float(rr_[cols[ckey]])
            egrid[ckey].append(abs(c - parse_fortran_real(ft)) / printed_half_unit(ft))
        k = math.sqrt(2.0 * mtot * float(rr_[cols["E_coll"]]) / CHAEV)
        egrid["momentum"].append(abs(k - parse_fortran_real(b["total_momentum"])) / printed_half_unit(b["total_momentum"]))
    out["energy_grid"] = {k: {"n": len(v), "worst_ratio": max(v), "n_within": int(sum(x <= 1 for x in v))} for k, v in egrid.items()}
    out["energy_grid"]["momentum"]["note"] = "derived on the C++ side: k = sqrt(2 M E_coll), M = total mass of the C++ log, legacy eV"
    gsum = (RUN_DIRS["gpu"] / "scattering_summary.tsv").read_text().splitlines()
    out["energy_grid"]["gpu_E_columns_equal"] = all(a.split("\t")[:3] == b.split("\t")[:3] for a, b in zip(
        (RUN_DIRS["cpu"] / "scattering_summary.tsv").read_text().splitlines(), gsum))

    # ---------------- GPU vs CPU set-up output: identical text up to the propagation section
    def setup_part(lines):
        start = next(i for i, l in enumerate(lines) if l.strip() == "INPUT")
        stop = next(i for i, l in enumerate(lines) if l.strip() == "PROPAGATION")
        keep = []
        skip = False
        for l in lines[start:stop]:
            if l.strip() == "GPU MEMORY ESTIMATE":
                skip = True
            elif skip and l.strip() in ("GAS-PHASE POTENTIAL",):
                skip = False
            if not skip and not l.strip().startswith("Data files written to"):
                keep.append(l)
        return keep
    cpu_setup, gpu_setup = setup_part(CL), setup_part(GL)
    diff_lines = [(a, b) for a, b in zip(cpu_setup, gpu_setup) if a != b]
    out["gpu_vs_cpu_setup_log"] = {"cpu_lines": len(cpu_setup), "gpu_lines": len(gpu_setup), "differing": len(diff_lines),
                                   "first_differences": diff_lines[:3],
                                   "excluded": "GPU MEMORY ESTIMATE block (GPU run only), 'Data files written to' line"}
    for r in ROWS:
        r["gpu_equal_to_cpu"] = out["gpu_vs_cpu_setup_log"]["differing"] == 0
    # item 11: the |psi(k)| tails outside the printed criterion are reproduced by item11_psi_k.py
    psi_k = load_item("item11_psi_k")
    for r in ROWS:
        key = {"|psi(k)| at k = 0": "zero", "|psi(k)| at the grid edge": "edge"}.get(r["quantity"])
        if key and r["verdict"] == "OUTSIDE" and psi_k.get(key, {}).get("within_fortran_half_unit_with_fortran_alpha_kz"):
            r["verdict"] = "explained"
            r["note"] = ("reproduced: the correctly rounded DFT of the samples with the Fortran's alpha and k_Z "
                         f"({psi_k['alpha_ulps']:.0f} and {psi_k['kz_ulps']:.0f} ulp from the C++ values, quad-precision site) "
                         "gives the Fortran value within its half unit (item11_psi_k.py)")
    out["rows"] = ROWS
    by_item = {}
    for r in ROWS:
        by_item.setdefault(r["item"], []).append(r["verdict"])
    out["verdicts_by_item"] = {str(k): {v: vs.count(v) for v in sorted(set(vs))} for k, vs in by_item.items()}
    save_item("setup", out)
    write_table(ROWS)
    print("json vs DVR.inp:", out["json_vs_dvr"]["n_equal"], "/", out["json_vs_dvr"]["n"], "unequal:", out["json_vs_dvr"]["unequal"], "dvr-only:", extra)
    print("Fortran echo vs input:", out["fortran_echo_vs_input"]["n_ok"], "/", out["fortran_echo_vs_input"]["n"], out["fortran_echo_vs_input"]["not_ok"])
    print("C++ echo vs input:", out["cpp_echo_vs_input"]["n_ok"], "/", out["cpp_echo_vs_input"]["n"])
    print("verdicts by item:", out["verdicts_by_item"])
    for r in ROWS:
        if r["verdict"] in ("OUTSIDE", "FINDING"):
            print("  ", r["item"], r["quantity"], r["fortran"], r["cpp"], r.get("ratio", r.get("relative")))
    print("vasy:", out["vasy"])
    print("rovib:", {k: out["rovib"][k] for k in ("n_states", "n_values", "sign_flipped_states")}, "outside:", len(out["rovib"]["outside"]))
    for q, w in out["rovib"]["worst_ratio_by_quantity"].items():
        print("   ", q, round(w["ratio"], 3), w["v"], w["j"])
    print("labels:", out["state_labels"], "\nchannel energies:", out["channel_energies"])
    print("energy grid:", out["energy_grid"])
    print("GPU vs CPU set-up log:", out["gpu_vs_cpu_setup_log"]["differing"], "differing of", out["gpu_vs_cpu_setup_log"]["cpu_lines"])


def write_table(rows):
    TABLES.mkdir(parents=True, exist_ok=True)
    lines = ["| # | Quantity | Fortran (label: value) | C++ (label: value) | Kind | Δ (C++ − F) | Tolerance | Verdict |",
             "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["kind"] == "printed":
            d, tol = f"{r['delta']:.2e}", f"{r['tolerance']:.2e}"
        elif r["kind"] == "full":
            d, tol = f"{r['delta']:.2e} ({r['ulps']:.0f} ulp)" if r["delta"] != 0 else "0", "1e-13 rel." if "1e-12" not in r["note"] else "1e-12 rel."
        elif r["kind"] == "int":
            d, tol = str(r["delta"]), "exact"
        else:
            d, tol = "", "exact" if r["kind"] == "text" else ""
        c = r["cpp"] if r["cpp"] is not None else "—"
        if isinstance(c, float):
            c = f"{c:.17g}"
        note = f" ({r['note']})" if r["note"] else ""
        cells = [str(r["item"]), r["quantity"], f"{r['fortran_label']}: `{r['fortran']}`", f"{r['cpp_label']}: `{c}`",
                 r["kind"], d, tol, r["verdict"] + note]
        lines.append("| " + " | ".join(x.replace("|", "\\|") for x in cells) + " |")
    (TABLES / "setup_label_map.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

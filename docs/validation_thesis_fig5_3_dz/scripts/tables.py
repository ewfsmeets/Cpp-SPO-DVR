#!/usr/bin/env python3
"""Markdown table fragments of the report, generated from scripts/results.json (tables/*.md).
fill_tables.py copies them into README.md, stages.md and appendix.md between the markers
<!-- BEGIN table:NAME --> and <!-- END table:NAME -->.
Markers in the cells: (p) printed precision (Fortran printed value vs C++ full precision),
(f) full precision on both sides."""
from __future__ import annotations

import json
import re

from common import RESULTS, TABLES, fmt

R = json.loads(RESULTS.read_text())
P1, P2 = " (p)", " (f)"


def table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    for r in rows:
        out.append("| " + " | ".join(re.sub(r"(?<!\\)\|", r"\\|", str(c)) for c in r) + " |")
    return "\n".join(out) + "\n"


def write(name, text):
    TABLES.mkdir(exist_ok=True)
    (TABLES / f"{name}.md").write_text(text)


def link(n):
    return f"[{n}](stages.md#item-{n.lower()})"


def verdict_counts(item):
    v = R["setup"]["verdicts_by_item"][str(item)]
    return ", ".join(f"{k} {n}" for k, n in v.items())


def main():
    s, i4, i13, i14 = R["setup"], R["item04"], R["item13"], R["item14_16"]
    i13p, i17, i18a, i18b, i19 = R["item13_psi"], R["item17"], R["item18a"], R["item18b"], R["item19"]
    po = i18b["propagation_only"]
    q = lambda k: i14["item14"][f"{k}_vs_fortran"]
    bk = i13["bk"]

    def pr(tally):   # "n_within+n_explained/n"
        return f"{tally['n_within_printed']}/{tally['n']} within ½ unit, {tally['n_explained']} explained"

    rows = [
        (link("1"), "Input read, decisions taken", "identical / printed",
         f"{s['fortran_echo_vs_input']['n_ok']}/{s['fortran_echo_vs_input']['n']} echoed values = input{P1}; "
         f"input.json ≡ DVR.inp {s['json_vs_dvr']['n_equal']}/{s['json_vs_dvr']['n']} keys (Title by design)", "same", "shared code",
         "exact; ½ printed unit"),
        (link("2"), "Grids (Z, r, x, y, k, θ, φ, J)", "identical / rounding",
         f"grids, φ, counts identical; θ nodes ≤ {i4['theta']['max_ulps']:.0f} ulp{P2}", "same", "shared code", "exact; few ulp"),
        (link("3"), "MISC SETUP lines", "identical / printed", verdict_counts(3), "same", "shared code", "exact; ½ printed unit"),
        (link("4"), "6D potential, full grid", "rounding",
         f"{fmt(i4['same_point']['max'])} at the same point, {fmt(i4['used']['max'])} as used{P2}", "same (V computed on the host)",
         "identical (bitwise)", "1e-13 · max(\\|V\\|, 1 eV)"),
        (link("5"), "Gas-phase V_asy(r)", "rounding", f"{fmt(s['vasy']['max_scaled'])}{P2}", "same", "identical", "1e-14 · max(\\|V\\|, 1 eV)"),
        (link("6"), "Rovibrational states", "printed / rounding",
         f"{s['rovib']['n_values']} values ≤ {max(w['ratio'] for w in s['rovib']['worst_ratio_by_quantity'].values()):.2f} ½ unit{P1}; "
         f"E(0,0) {fmt([r for r in s['rows'] if r['quantity'].startswith('initial rovibrational energy')][0]['relative'])} rel.{P2}",
         "same", "identical", "½ printed unit; 1e-12 rel."),
        (link("7"), "Diffraction basis (incl. warning)", "identical", verdict_counts(7), "same", "shared code", "exact"),
        (link("8"), "Rotation channels", "identical", verdict_counts(8), "same", "shared code", "exact"),
        (link("9"), "Rotation-diffraction channels", "identical / rounding",
         f"labels of {s['state_labels']['n_states']} states identical; channel energies ≤ {fmt(s['channel_energies']['max_relative'])} rel.{P2}",
         "same", "identical", "exact; 1e-13 rel."),
        (link("10"), "Optical potentials", "printed", verdict_counts(10), "same", "shared code", "½ printed unit"),
        (link("11"), "Initial wave packet", "printed / rounding / explained", verdict_counts(11), "same", "shared code",
         "½ printed unit; 1e-13 rel."),
        (link("12"), "Analysis set-up, energy grid", "identical / printed",
         f"{verdict_counts(12)}; {s['energy_grid']['E_coll']['n_within']}/{s['energy_grid']['E_coll']['n']} energies{P1}", "same",
         "identical", "exact; ½ printed unit"),
        (link("13a"), "Averages every 500 a.u.", "printed",
         f"{i13['averages']['cpu_vs_fortran']['n_values'] - i13['averages']['cpu_vs_fortran']['n_outside']}/{i13['averages']['cpu_vs_fortran']['n_values']} "
         f"within ½ unit{P1}", f"{i13['averages']['gpu_vs_fortran']['n_values'] - i13['averages']['gpu_vs_fortran']['n_outside']}/"
         f"{i13['averages']['gpu_vs_fortran']['n_values']}{P1}", f"norm {fmt(i13['averages']['gpu_vs_cpu']['norm']['max'])} rel.{P2}",
         "½ printed unit; 1e-10"),
        (link("13b"), "Specular analysis line fixed", "identical",
         f"same step (t = {i13['analysis_line']['time_of_that_record']['fortran']:.0f} a.u.)", "same step", "same step", "exact"),
        (link("13c"), "BK coefficients, 875 steps", "rounding · open (channel-local excess)", f"{fmt(bk['cpu-fortran']['max_over_scale'])}{P2}",
         f"{fmt(bk['gpu-fortran']['max_over_scale'])}{P2}", f"{fmt(bk['gpu-cpu']['max_over_scale'])}{P2}", "1e-10 · max\\|c\\|"),
        (link("13d"), "Final wave function", "rounding",
         f"‖Δψ‖/‖ψ_T‖ {fmt(i13p['cpu-fortran']['rel_l2'])}, {fmt(i13p['cpu-fortran']['scattering_abs_rel_initial_norm'])} of the initial norm{P2}",
         f"{fmt(i13p['gpu-fortran']['scattering_abs_rel_initial_norm'])}{P2}", f"{fmt(i13p['gpu-cpu']['scattering_abs_rel_initial_norm'])}{P2}",
         "reported (no threshold)"),
        (link("14"), "P_reaction, ΣP, P(v), P(v, j, mJ)", "printed / explained / rounding",
         f"P_reaction {pr(q('cpu')['P_reaction'])}{P1}; full precision ≤ {fmt(po['cpu-fortran']['summary']['P_reaction']['max'])}{P2}",
         f"≤ {fmt(po['gpu-fortran']['summary']['P_reaction']['max'])}{P2}", f"≤ {fmt(i14['item14']['gpu_vs_cpu']['P_reaction']['max'])}{P2}",
         "½ printed unit, floor 1e-11 for P_reaction; 1e-10"),
        (link("15"), "Diffraction-order probabilities", "printed / explained",
         f"{pr(i14['item15']['cpu_vs_fortran'])}{P1}", f"{pr(i14['item15']['gpu_vs_fortran'])}{P1}",
         f"{fmt(i14['item15']['gpu_vs_cpu']['max'])}{P2}", "½ printed unit; 1e-10"),
        (link("16"), "S-matrix elements", "printed / explained / rounding · open (channel-local excess)",
         f"\\|S\\|² {i14['item16']['cpu_vs_fortran']['|S|^2']['n_within_printed']}/{i14['item16']['cpu_vs_fortran']['|S|^2']['n']} within ½ unit{P1}; "
         f"\\|ΔS\\| ≤ {fmt(po['cpu-fortran']['s_matrix']['S_abs']['max'])}{P2}", f"\\|ΔS\\| ≤ {fmt(po['gpu-fortran']['s_matrix']['S_abs']['max'])}{P2}",
         f"\\|ΔS\\| ≤ {fmt(i14['item16']['gpu_vs_cpu']['S_abs']['max'])}{P2}", "½ printed unit; 1e-10"),
        (link("17"), "Thesis Fig. 5.3 (dotted)", "agrees",
         f"{i17['cpu']['n_within']}/{i17['n_points_in_range']} points; max \\|ΔP\\| {i17['cpu']['max_abs_dP']:.3f}",
         f"{i17['gpu']['n_within']}/{i17['n_points_in_range']}", "—", "digitisation + 0.01"),
        (link("18a"), "Fine analysis (onlyAnalysis) of each BK file", "rounding",
         f"P_reaction ≤ {fmt(i18a['cpu-fortran']['P_reaction']['max'])}{P2}", f"≤ {fmt(i18a['gpu-fortran']['P_reaction']['max'])}{P2}",
         f"≤ {fmt(i18a['gpu-cpu']['P_reaction']['max'])}{P2}", "1e-10"),
        (link("18b"), "Replay through one analysis code", "rounding / printed",
         f"propagation only: \\|ΔP\\| ≤ {fmt(po['cpu-fortran']['s_matrix']['P_abs']['max'])}{P2}; analysis only: all printed P within ½ unit{P1}",
         f"\\|ΔP\\| ≤ {fmt(po['gpu-fortran']['s_matrix']['P_abs']['max'])}{P2}", f"\\|ΔP\\| ≤ {fmt(po['gpu-cpu']['s_matrix']['P_abs']['max'])}{P2}",
         "1e-10; ½ printed unit"),
        (link("19"), "Physics checks", "same in all codes",
         f"thresholds, bounds, P(mJ) = P(−mJ) ≤ {fmt(max(v['mJ_symmetry']['max_relative'] for v in i19['symmetry'].values()))} rel.",
         "same", "same", "property shared by all three"),
    ]
    write("verdict", table(["#", "Item", "Verdict", "C++ CPU − Fortran", "C++ GPU − Fortran", "C++ GPU − C++ CPU", "Criterion"], rows))

    # provenance
    p0 = R["phase0"]
    b, r = p0["builds"], p0["runs"]
    prov = [
        ("C++ CPU", f"`{b['cpu']['dir']}/spodvr` ({b['cpu']['spodvr_mtime'][:16].replace('T', ' ')})",
         f"FAST_MATH {b['cpu']['cache']['SPODVR_FAST_MATH']}, ARCH {b['cpu']['cache']['SPODVR_ARCH']}, oneMKL", r["cpu"]["start"][:16],
         r["cpu"]["time"]["wall"], "`runs/thesis_fig5_3_dz_cpu_fast/`"),
        ("C++ GPU", f"`{b['gpu']['dir']}/spodvr` ({b['gpu']['spodvr_mtime'][:16].replace('T', ' ')})",
         f"FAST_MATH {b['gpu']['cache']['SPODVR_FAST_MATH']}, ARCH {b['gpu']['cache']['SPODVR_ARCH']}, CUDA, cuFFT", r["gpu"]["start"][:16],
         r["gpu"]["time"]["wall"], "`runs/thesis_fig5_3_dz_cuda_fast/`"),
        ("Fortran 5.19", f"`{b['fortran']['exe']}` ({b['fortran']['mtime'][:16].replace('T', ' ')})", b["fortran"]["flags"],
         r["fortran"]["start"][:16], r["fortran"]["time"]["wall"], "`runs/thesis_fig5_3_dz_fortran/`"),
    ]
    write("provenance", table(["Run", "Executable (built)", "Build", "Started", "Wall", "Output"], prov))

    # parameters (below Figure 1)
    js = R["input"]
    par = [("PES", "Dai–Zhang LEPS (`H2Cu111DaiLight`), lattice constant 4.84 a₀, 60° cell, C6v symmetry used"),
           ("Initial state", f"H₂ v = {js['VIn']}, j = {js['JIn']}, mJ = {js['MJIn']}, normal incidence"),
           ("Z grid", f"{js['nZ']} points from {js['ZStart']} a₀, ΔZ {js['dZ']} a₀; specular grid {js['nZSpecular']} points"),
           ("r grid", f"{js['nR']} points from {js['RStart']} a₀, Δr {js['dR']} a₀"),
           ("x, y grid", f"{js['nX']} × {js['nY']}"),
           ("Rotation", f"j_max = m_j,max = {js['JMax']}: 435 functions; 15 θ, 60 φ points"),
           ("Time", f"Δt {js['dtStep']} a.u., {js['nStep']} steps to T = {js['tEnd']:.0f} a.u.; analysis every {js['dtAna']} a.u."),
           ("Wave packet", f"{js['EMinWP']}–{js['EMaxWP']} eV, Z₀ = {js['Z0']} a₀; {js['NrOfEs']} analysis energies"),
           ("Analysis line", f"Z_ana = {js['Zana']} a₀ (dynamic specular line)"),
           ("Absorbers", f"Z {js['ZOptStart']}–{js['ZOptEnd']} ({js['ZOptEnergy']} eV), Z_spec {js['ZSpOptStart']}–{js['ZSpOptEnd']} "
                         f"({js['ZSpOptEnergy']} eV), r {js['ROptStart']}–{js['ROptEnd']} ({js['ROptEnergy']} eV)"),
           ("Input", "`inputs/thesis_fig5_3_dz/input.json` (C++), `DVR.inp` (Fortran); choices: `inputs/thesis_fig5_3_dz/README.md`")]
    write("parameters", table(["Parameter", "Value"], par))

    # ---- per-item tables for stages.md
    rows = [(x["label"], x["key"], f"`{x['fortran']}`", x["input"], "yes" if x["ok"] else "**no**") for x in s["fortran_echo_vs_input"]["rows"]]
    write("item01_echo", table(["Fortran echo", "Key", "Fortran printed", "Input", "Within"], rows))
    write("item01", table(["Check", "Result"], [
        ("`input.json` ≡ `DVR.inp`, field by field", f"{s['json_vs_dvr']['n_equal']}/{s['json_vs_dvr']['n']} keys equal; differs: "
                                                   f"{', '.join(s['json_vs_dvr']['unequal'])} ({s['title_lengths']['dvr']} vs {s['title_lengths']['json']} characters)"),
        ("Fortran input echo = input values", f"{s['fortran_echo_vs_input']['n_ok']}/{s['fortran_echo_vs_input']['n']} within ½ printed unit"),
        ("C++ log echo and BK-file restart parameters = `input.json`", f"{s['cpp_echo_vs_input']['n_ok']}/{s['cpp_echo_vs_input']['n']} checks"),
        ("Set-up rows (label map)", verdict_counts(1)),
        ("GPU set-up log = CPU set-up log", f"{s['gpu_vs_cpu_setup_log']['cpu_lines'] - s['gpu_vs_cpu_setup_log']['differing']}/"
                                            f"{s['gpu_vs_cpu_setup_log']['cpu_lines']} lines identical")]))
    ge = i4["grid_ends"]
    write("item02", table(["Quantity", "Fortran", "C++", "Result"], [
        ("Z, r, x, y end points (`.pes` header, full precision)", ", ".join(f"{k} {v['fortran']:.6g}" for k, v in ge.items()),
         "same", f"{sum(v['diff'] == 0 for v in ge.values())}/{len(ge)} bitwise identical"),
        ("θ nodes (15, full precision)", "`.pes`", "helper grids.txt", f"{i4['theta']['n_identical']}/15 identical, max {fmt(i4['theta']['max_abs_diff'])} rad "
                                                                     f"({i4['theta']['max_ulps']:.0f} ulp)"),
        ("φ nodes (60)", "`.pes`", "helper grids.txt", f"{i4['phi']['n_identical']}/60 identical"),
        ("nφ, nθ, J range, (x, y) points computed", "60, 15, 0..28 step 2, 19", "same", "identical"),
        ("Momentum grids dk_Z, dk_r, dk_x, dk_y", "6 decimals", "%.16e", verdict_counts(3).split(",")[0] + " (item 3)")]))
    misc = [r for r in s["rows"] if r["item"] == 3]
    write("item03", table(["Fortran line", "C++", "Kind", "Verdict"], [
        (f"`{r['fortran_label']} = {r['fortran']}`", f"{r['cpp_label']}: `{r['cpp']}`", r["kind"], r["verdict"]) for r in misc]))
    write("item04", table(["Comparison", "n", "median", "p95", "max", "where (max)"], [
        ("Same physical point: Fortran stored value vs C++ PES", i4["same_point"]["n"], fmt(i4["same_point"]["median"]),
         fmt(i4["same_point"]["p95"]), fmt(i4["same_point"]["max"]),
         "Z {Z:.2f}, r {r:.2f}, V {V_eV:.3f} eV".format(**i4["same_point"]["argmax"])),
        ("… restricted to V < 2 eV", i4["same_point_v_below_2eV"]["n"], fmt(i4["same_point_v_below_2eV"]["median"]),
         fmt(i4["same_point_v_below_2eV"]["p95"]), fmt(i4["same_point_v_below_2eV"]["max"]), ""),
        ("As used: Fortran expanded (symmetryIndex) vs C++ gathered (own maps)", i4["used"]["n"], fmt(i4["used"]["median"]),
         fmt(i4["used"]["p95"]), fmt(i4["used"]["max"]), "Z {Z:.2f}, r {r:.2f}, V {V_eV:.3f} eV".format(**i4["used"]["argmax"])),
        ("PES asymmetry, C++ maps: V(q) vs V(rep(q))", i4["asym_cpp_maps"]["n"], fmt(i4["asym_cpp_maps"]["median"]),
         fmt(i4["asym_cpp_maps"]["p95"]), fmt(i4["asym_cpp_maps"]["max"]), ""),
        ("PES asymmetry, Fortran maps", i4["asym_fortran_maps"]["n"], fmt(i4["asym_fortran_maps"]["median"]),
         fmt(i4["asym_fortran_maps"]["p95"]), fmt(i4["asym_fortran_maps"]["max"]), "")]))
    v = s["vasy"]
    write("item05", table(["Quantity", "Result"], [("points", v["n"]), ("r grid", "identical" if v["r_identical"] else "differs"),
                                                  ("max \\|ΔV\\| / max(\\|V\\|, 1 eV)", f"{fmt(v['max_scaled'])} at r = {v['argmax_r']} a₀"),
                                                  ("bitwise identical values", f"{v['n_identical']}/{v['n']}"),
                                                  ("GPU run file = CPU run file", "yes" if v["gpu_file_identical"] else "no")]))
    w = s["rovib"]["worst_ratio_by_quantity"]
    write("item06", table(["Quantity (58 states)", "Fortran format", "Worst \\|Δ\\| / ½ unit", "at (v, j)"], [
        (q_, fm, f"{w[q_]['ratio']:.3f}", f"({w[q_]['v']}, {w[q_]['j']})") for q_, fm in
        (("e", "F 6 decimals (eV)"), ("erel", "F 6 decimals"), ("psi1", "D, 7 digits"), ("psinr", "D, 7 digits"),
         ("norm", "F 7 decimals"), ("r", "F 7 decimals"))]))
    rows = [(r["quantity"], f"`{r['fortran']}`", r["cpp"], r["verdict"]) for r in s["rows"] if r["item"] in (7, 8, 9)]
    write("item07_09", table(["Quantity", "Fortran", "C++", "Verdict"], rows))
    ce = s["channel_energies"]
    write("item09_energies", table(["Quantity", "Result"], [
        ("state labels (J, mJ, n, m) of all states", f"{s['state_labels']['n_states']} × identical: {s['state_labels']['equal']}"),
        ("channel energies (fort.21 enymo vs bk_coefficients.bin)", f"{ce['n']} values, max \\|Δ\\| {fmt(ce['max_abs_hartree'])} hartree "
                                                                   f"({fmt(ce['max_relative'])} rel.), median {fmt(ce['median_abs'])}"),
        ("GPU = CPU", "identical" if ce["gpu_identical_to_cpu"] else "differs")]))
    rows = [(r["quantity"], f"`{r['fortran']}`", f"{r['cpp']:.10g}" if isinstance(r["cpp"], float) else r["cpp"],
             f"{r.get('ratio', 0):.2f}" if r["kind"] == "printed" else (f"{r.get('ulps', 0):.0f} ulp" if r["kind"] == "full" else ""),
             r["verdict"]) for r in s["rows"] if r["item"] in (10, 11, 12)]
    write("item10_12", table(["Quantity", "Fortran", "C++", "\\|Δ\\|/½ unit or ulp", "Verdict"], rows))
    pk = R["item11_psi_k"]
    write("item11_psik", table(["\\|ψ(k)\\|", "Fortran (FFTW)", "C++ (oneMKL)", "exact DFT, C++ α, k_Z", "exact DFT, Fortran α, k_Z"], [
        (k, pk[k]["fortran_text"], f"{pk[k]['cpp_onemkl']:.10e}", f"{pk[k]['exact_fsum']:.10e}",
         f"{pk[k]['exact_fsum'] + pk[k]['effect_of_fortran_alpha_kz']:.10e}") for k in ("zero", "edge")]))
    av = i13["averages"]
    write("item13a", table(["Column", "C++ CPU vs Fortran (worst / ½ unit)", "C++ GPU vs Fortran", "GPU − CPU (full)"], [
        (c, f"{av['cpu_vs_fortran']['worst_ratio_by_column'][c]:.3f}", f"{av['gpu_vs_fortran']['worst_ratio_by_column'][c]:.3f}",
         f"{fmt(av['gpu_vs_cpu'][c]['max'])} {av['gpu_vs_cpu'][c]['metric'][:3]}.") for c in ("norm", "Z", "r", "rest_norm")]))
    al = i13["analysis_line"]
    write("item13b", table(["t (a.u.)", "\\|c_init\\| Fortran", "\\|c_init\\| C++ CPU", "\\|c_init\\| C++ GPU"],
                           [(f"{t:.0f}", f"{a:.3e}", f"{b_:.3e}", f"{c:.3e}") for t, a, b_, c in
                            zip(al["initial_channel_abs_around_fix"]["t"], al["initial_channel_abs_around_fix"]["fortran"],
                                al["initial_channel_abs_around_fix"]["cpu"], al["initial_channel_abs_around_fix"]["gpu"])]))
    av_ = bk["amplitude_vs_time"]
    write("item13c", table(["Pair", "max \\|Δc\\|", "/ max\\|c\\|", "at", "per-channel rel. (median / max)", "‖Δc‖/‖c‖ in window (median)",
                            "amplitude part κ(t), window quartile medians"], [
        (k, fmt(bk[k]["max_abs"]), fmt(bk[k]["max_over_scale"]), "t {t:.0f}, J {J} (n,m) ({n},{m})".format(**bk[k]["argmax"]),
         f"{fmt(bk[k]['per_channel_relative']['median'])} / {fmt(bk[k]['per_channel_relative']['max'])}",
         fmt(bk[k]["l2"]["rel_in_window"]["median"]),
         " / ".join(fmt(x) for x in av_[k]["window_quartile_medians"]))
        for k in ("cpu-fortran", "gpu-fortran", "gpu-cpu")]))
    fam = bk["specular_inelastic_family"]
    rows = []
    for jj in ("2", "4", "6", "8", "10"):
        rows.append((jj, fmt(fam["cpu-fortran"]["per_J"][jj]["max_abs_c"]), fmt(fam["cpu-fortran"]["per_J"][jj]["max_abs_dc"]),
                     fmt(fam["cpu-fortran"]["per_J"][jj]["rel"]), fmt(fam["gpu-cpu"]["per_J"][jj]["max_abs_dc"]),
                     fmt(fam["gpu-cpu"]["per_J"][jj]["rel"])))
    write("item13c_family", table(["j (v=0, (n,m)=(0,0), m_J=0)", "max \\|c\\|", "C++ − F: max \\|Δc\\|", "relative",
                                   "GPU − CPU: max \\|Δc\\|", "relative"], rows))
    dS = i18b["amplitude_decomposition"]["cpu-fortran"]
    dG = i18b["amplitude_decomposition"]["gpu-cpu"]
    rel_all = i18b["S_relative_normweighted"]
    fj = dS["specular_inelastic_family"]["per_J"]
    gj = dG["specular_inelastic_family"]["per_J"]
    write("composition", table(["Part of the S-matrix difference (propagation only; relative to \\|S\\|)", "C++ CPU − Fortran", "C++ GPU − C++ CPU"], [
        ("all channels, ‖ΔS‖/‖S‖ (norm-weighted over all channels and energies)", fmt(rel_all["cpu-fortran"]), fmt(rel_all["gpu-cpu"])),
        ("large channels (\\|S\\|² > 1e-2), \\|ΔS\\|/\\|S\\| median", fmt(dS["rel_dS_large_channels_median"]), fmt(dG["rel_dS_large_channels_median"])),
        ("uniform amplitude part κ(E), one factor for every channel at an energy", f"{fmt(dS['alpha_min'])} … {fmt(dS['alpha_max'])}",
         f"{fmt(dG['alpha_min'])} … {fmt(dG['alpha_max'])}"),
        ("phase part of that factor", f"≤ {fmt(dS['beta_abs_max'])}", f"≤ {fmt(dG['beta_abs_max'])}"),
        ("what remains per channel after removing it: share of the norm ‖ΔS‖ (of ‖ΔS‖²)",
         f"{100 * dS['specular_inelastic_family']['residual_norm_fraction_of_total_dS']:.0f}% "
         f"({100 * dS['specular_inelastic_family']['residual_norm2_fraction_of_total']:.0f}%)",
         f"{100 * dG['specular_inelastic_family']['residual_norm_fraction_of_total_dS']:.0f}% "
         f"({100 * dG['specular_inelastic_family']['residual_norm2_fraction_of_total']:.0f}%)"),
        ("… diffuse, outside the specular inelastic channels", fmt(dS["specular_inelastic_family"]["residual_rel_normweighted_outside_family"]),
         fmt(dG["specular_inelastic_family"]["residual_rel_normweighted_outside_family"])),
        ("… specular inelastic channels (v=0, (n,m)=(0,0), m_J=0), j = 2 / 4 / 6 / 8 / 10",
         " / ".join(fmt(fj[j]["rel_residual_normweighted"]) for j in ("2", "4", "6", "8", "10")),
         " / ".join(fmt(gj[j]["rel_residual_normweighted"]) for j in ("2", "4", "6", "8", "10"))),
        ("… the same, single elements \\|ΔS\\|/\\|S\\|: largest per j",
         " / ".join(fmt(dS["specular_inelastic_family"]["per_element_relative"][j]["max"]) for j in ("2", "4", "6", "8", "10")),
         " / ".join(fmt(dG["specular_inelastic_family"]["per_element_relative"][j]["max"]) for j in ("2", "4", "6", "8", "10")))]))
    pr_ = R["item18b"]["propagation_only"]
    write("prediction", table(["Pre-registered prediction", "Outcome", "Met?"], [
        ("set-up integers, decisions, channel lists identical", "identical (items 1–3, 7–9, 12)", "yes"),
        ("potential at identical points ≤ 1e-13", f"{fmt(i4['same_point']['max'])}", "yes"),
        ("C++↔Fortran BK ≈ 1e-12, at most 1e-11 of max\\|c\\|", f"{fmt(bk['cpu-fortran']['max_over_scale'])}",
         "within the bound, above the central estimate"),
        ("growth ∝ √t to ∝ t, no jump", "no jump at the fixing step; the growth law is not identifiable (window spans a factor 2.2 in t)", "partly"),
        ("C++↔Fortran \\|ΔP_reaction\\| ≲ 1e-11", f"{fmt(pr_['cpu-fortran']['summary']['P_reaction']['max'])} (CPU), "
                                                  f"{fmt(pr_['gpu-fortran']['summary']['P_reaction']['max'])} (GPU)", "**no**, narrowly"),
        ("C++↔Fortran within ~10× of CPU↔GPU", f"{rel_all['cpu-fortran'] / rel_all['gpu-cpu']:.0f}× (‖ΔS‖/‖S‖), "
                                               f"{R['ladder']['gap_ratios_of_medians'][0]:.0f}× (P_reaction medians)", "**no**, narrowly"),
        ("… also per channel", f"specular inelastic channels {R['derived']['family_S_cpuF_over_gpuC_min']:.0f}–"
                               f"{R['derived']['family_S_cpuF_over_gpuC_max']:.0f}× (S), "
                               f"{R['derived']['family_cpuF_over_gpuC_min']:.0f}–{R['derived']['family_cpuF_over_gpuC_max']:.0f}× (BK)", "**no**"),
        ("printed set-up reals within printed precision", "all but the two \\|ψ(k)\\| tails (337 and 346 half units, explained, item 11)", "partly"),
        ("printed comparisons within ½ unit except small values", f"{R['derived']['explained_printed_total_cpu']} values outside (CPU), "
         "all explained by a rule adopted after seeing them", "partly"),
        ("analysis line fixed at the same step", "same step, t = 2220 a.u.", "yes")]))
    write("item13d", table(["Pair", "‖Δψ‖/‖ψ_T‖ scattering grid", "… of the initial norm", "specular grid", "global phase (rad)"], [
        (k, fmt(i13p[k]["rel_l2"]), fmt(i13p[k]["scattering_abs_rel_initial_norm"]), fmt(i13p[k]["specular_rel_l2"]),
         fmt(i13p[k]["global_phase_rad"])) for k in ("cpu-fortran", "gpu-fortran", "gpu-cpu")]))
    rows = []
    for name in ("P_reaction", "P_sum", "P_v", "P_vjm", "E_avg"):
        c_, g_ = q("cpu")[name], q("gpu")[name]
        where = ", ".join(f"{kk} {vv}" for kk, vv in c_["worst"]["where"].items())
        rows.append((name, c_["n"], f"{c_['n_within_printed']} / {c_['n_explained']} / {c_['n_outside']}",
                     f"{g_['n_within_printed']} / {g_['n_explained']} / {g_['n_outside']}", f"{c_['worst']['ratio']:.2f}", where))
    write("item14", table(["Quantity", "n", "C++ CPU: within / explained / outside", "C++ GPU", "worst ratio (CPU)", "at"], rows))
    write("item14_full", table(["Pair (full precision)", "P_reaction max", "median", "P(v,j,mJ) max", "source"], [
        ("C++ CPU − Fortran", fmt(po["cpu-fortran"]["summary"]["P_reaction"]["max"]), fmt(po["cpu-fortran"]["summary"]["P_reaction"]["median"]),
         fmt(po["cpu-fortran"]["rovib"]["max"]), "replay (item 18b)"),
        ("C++ GPU − Fortran", fmt(po["gpu-fortran"]["summary"]["P_reaction"]["max"]), fmt(po["gpu-fortran"]["summary"]["P_reaction"]["median"]),
         fmt(po["gpu-fortran"]["rovib"]["max"]), "replay (item 18b)"),
        ("C++ GPU − C++ CPU", fmt(i14["item14"]["gpu_vs_cpu"]["P_reaction"]["max"]), fmt(i14["item14"]["gpu_vs_cpu"]["P_reaction"]["median"]),
         fmt(i14["item14"]["gpu_vs_cpu"]["P_vjm"]["max"]), "run tables")]))
    t15 = i14["item15"]
    write("item15", table(["Pair", "n", "within ½ unit", "explained", "outside", "worst ratio", "at"], [
        (k, t15[f"{k}_vs_fortran"]["n"], t15[f"{k}_vs_fortran"]["n_within_printed"], t15[f"{k}_vs_fortran"]["n_explained"],
         t15[f"{k}_vs_fortran"]["n_outside"], f"{t15[f'{k}_vs_fortran']['worst']['ratio']:.2f}",
         ", ".join(f"{kk} {vv}" for kk, vv in t15[f"{k}_vs_fortran"]["worst"]["where"].items())) for k in ("cpu", "gpu")]))
    t16 = i14["item16"]["cpu_vs_fortran"]
    loc = lambda w: ", ".join(f"{kk} {(round(vv, 3) if isinstance(vv, float) else vv)}" for kk, vv in w.items())
    write("item16", table(["Column (C++ CPU vs fort.25)", "n", "within ½ unit", "explained", "outside", "worst ratio", "at"], [
        (c, t16[c]["n"], t16[c]["n_within_printed"], t16[c]["n_explained"], t16[c]["n_outside"], f"{t16[c]['worst_ratio']:.2f}",
         loc(t16[c]["worst"])) for c in ("Re S", "Im S", "|S|^2")]
        + [("phase (\\|S\\|² > 1e-8)", t16["phase"]["n"], t16["phase"]["n_within_printed"], 0, 0,
            f"{t16['phase']['worst_ratio']:.2f}", loc(t16["phase"]["worst"]))]))
    write("item17", table(["Run", "points in range", "within tolerance", "max \\|ΔP\\|", "at E (eV)", "mean ΔP (run − thesis)",
                           "best shift"], [(k, i17["n_points_in_range"], i17[k]["n_within"], f"{i17[k]['max_abs_dP']:.4f}",
                                            f"{i17[k]['at_E']:.3f}", f"{i17[k]['mean_dP']:+.4f}", f"{1000 * i17[k]['best_shift_eV']:+.1f} meV")
                                           for k in ("fortran", "cpu", "gpu")]))
    write("item18a", table(["Pair", "P_reaction max", "median", "P(v,j,mJ) max", "\\|ΔS\\| max", "within 1e-10"], [
        (k, fmt(i18a[k]["P_reaction"]["max"]), fmt(i18a[k]["P_reaction"]["median"]), fmt(i18a[k]["P_vjm"]["max"]),
         fmt(i18a[k]["S"]["max"]), "yes" if i18a[k]["within_1e-10"] else "no") for k in ("cpu-fortran", "gpu-fortran", "gpu-cpu")]))
    ao = i18b["analysis_only"]
    write("item18b", table(["Analysis only: Fortran printed vs replay(fort.21)", "n", "within ½ unit"], [
        (k, v_["n"], v_.get("n_within", "")) for k, v_ in ao.items()
        if isinstance(v_, dict) and "n" in v_ and "corrected" not in k and "bound" not in k]))
    sc = i18b["self_check"]
    dec = i18b["amplitude_decomposition"]
    write("item18b_prop", table(["Pair (propagation only, same analysis binary)", "\\|ΔS\\| max", "\\|ΔP\\| max", "P_reaction max",
                                 "amplitude κ (median)", "phase β (max)"], [
        (k, fmt(po[k]["s_matrix"]["S_abs"]["max"]), fmt(po[k]["s_matrix"]["P_abs"]["max"]), fmt(po[k]["summary"]["P_reaction"]["max"]),
         fmt(dec[k]["alpha_median"]) if k in dec else "—", fmt(dec[k]["beta_abs_max"]) if k in dec else "—")
        for k in ("cpu-fortran", "gpu-fortran", "gpu-cpu")]))
    write("item18b_self", table(["Run", "bitwise", "run − replay \\|ΔS\\| max", "\\|ΔP\\| max", "\\|ΔP_reaction\\| max"], [
        (k, "yes" if sc[k]["bitwise"] else "no", fmt(sc[k]["s_matrix"]["S_abs"]["max"]), fmt(sc[k]["s_matrix"]["P_abs"]["max"]),
         fmt(sc[k]["summary"]["P_reaction"]["max"])) for k in ("cpu", "gpu")]))
    sy = i19["symmetry"]
    th = i19["thresholds"]
    write("item19", table(["Check", "Fortran", "C++ CPU", "C++ GPU"], [
        ("open channels (E − E_ch > 0)", f"{th['n_open_fortran_by_energy_rule']}", f"{th['n_open_cpp_by_energy_rule']} (listed = rule: {th['cpu_listed_equals_rule']})",
         f"listed = rule: {th['gpu_listed_equals_rule']}"),
        ("smallest threshold margin (eV)", fmt(th["smallest_margin_eV"]["fortran"], 3), fmt(th["smallest_margin_eV"]["cpp"], 3), "same"),
        ("min P, max ΣP", f"{fmt(i19['bounds']['fortran']['min_P'])}, {i19['bounds']['fortran']['max_sum_P']}",
         f"{fmt(i19['bounds']['cpu']['min_P'])}, {i19['bounds']['cpu']['max_sum_P']:.7f}",
         f"{fmt(i19['bounds']['gpu']['min_P'])}, {i19['bounds']['gpu']['max_sum_P']:.7f}"),
        ("P(v,j,mJ) = P(v,j,−mJ), max rel. difference", fmt(sy["fortran"]["mJ_symmetry"]["max_relative"]),
         fmt(sy["cpu"]["mJ_symmetry"]["max_relative"]), fmt(sy["gpu"]["mJ_symmetry"]["max_relative"])),
        ("C6v operations other than (n, m) ↔ (m, n): median rel. difference (C++ CPU)", "",
         f"{100 * i19['per_operation_summary']['broken_median_min']:.1f}–{100 * i19['per_operation_summary']['broken_median_max']:.1f}%", ""),
        ("mirror (n, m) ↔ (m, n), max rel. difference (Fortran: printed values)", fmt(sy["fortran"]["diagonal_mirror"]["max_relative"]),
         fmt(sy["cpu"]["diagonal_mirror"]["max_relative"]), fmt(sy["gpu"]["diagonal_mirror"]["max_relative"])),
        ("C6v orbits of (n, m), rel. spread median / max", *(f"{sy[k]['diffraction_orbits']['median_relative_spread']:.3f} / "
                                                             f"{sy[k]['diffraction_orbits']['max_relative_spread']:.4f}" for k in ("fortran", "cpu", "gpu"))),
        ("P_reaction at E ≤ 0.50 eV (sub-threshold floor)", f"{fmt(i19['subthreshold']['fortran_same_range'][0])}–{fmt(i19['subthreshold']['fortran_same_range'][1])}",
         f"{fmt(i19['subthreshold']['P_reaction_min'])}–{fmt(i19['subthreshold']['P_reaction_max'])}", "same to 1e-11"),
        ("norm left on the grids at T", i19["norm_left_at_T"]["fortran_text"], fmt(i19["norm_left_at_T"]["cpu"], 6), fmt(i19["norm_left_at_T"]["gpu"], 6))]))

    # ---- appendix tables
    rows = [(o["component"], f"{o['E_coll']:.3f}", o["v"], o["j"], o["mJ"], f"({o['n']}, {o['m']})", f"`{o['fortran']}`", f"{o['replay']:.9e}",
             fmt(o["abs_diff"]), fmt(o["half_unit"]), fmt(o["bound"]), "yes" if o["explained"] else "**no**")
            for o in ao["analysis_side_S_outliers"]]
    write("appendix_s_outliers", table(["S", "E (eV)", "v", "j", "mJ", "(n, m)", "fort.25", "replay(fort.21)", "\\|Δ\\|", "½ unit", "bound",
                                        "explained"], rows))
    lines = (RESULTS.parent.parent / "data" / "reaction_probability.tsv").read_text().splitlines()[2:]
    rep = (RESULTS.parent.parent / "data" / "replay_reaction_probability.tsv").read_text().splitlines()[2:]
    rows = []
    for l, r_ in zip(lines, rep):
        a, b_ = l.split("\t"), r_.split("\t")
        rows.append((a[0], f"`{a[1]}`", f"{float(a[4]):.10e}", f"{float(a[5]):.10e}", b_[4], b_[5], a[8]))
    write("appendix_p_reaction", table(["E (eV)", "Fortran (printed)", "C++ CPU", "C++ GPU", "\\|CPU − F\\| (replay, full)",
                                        "\\|GPU − F\\| (replay, full)", "\\|GPU − CPU\\| (full)"], rows))
    print("tables written:", len(list(TABLES.glob("*.md"))))


if __name__ == "__main__":
    main()

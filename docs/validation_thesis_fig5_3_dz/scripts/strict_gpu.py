#!/usr/bin/env python3
"""Strict-build GPU run (criteria.md, 2026-11-08): the same input with build-cuda (no fast math on
the host, no FMA in the kernels) in runs/thesis_fig5_3_dz_cuda/, compared with the estimators of
items 13a, 13c, 14 and 18b, with the strict run in the place of the fast GPU run:

  validity      exit status, steps, build options, source age, analysis line, set-up files
  averages      averages.dat vs the Fortran averages.out, half a printed unit (item 13a)
  BK            bk_compare of item 13c on (Fortran, fast CPU, strict) and (Fortran, fast GPU, strict)
  replay        replay(strict BK file) with the same build-fast helper as the other runs: s_diff and
                decompose of item 18b against replay(fort.21), replay(fast CPU), replay(fast GPU)
  self-check    the run's tables vs replay(strict BK file) with a helper linked against build-cuda
  end to end    P_reaction vs the Fortran's printed values (item 14)

Needs `make_all.sh --strict` first (helper, replays, cleaned log). Writes data/items/strict_gpu.json and,
for Figures 1–3, data/strict_reaction_probability.tsv and data/strict_bk_error_vs_time.tsv."""
from __future__ import annotations

import datetime as dt
import re

import numpy as np
import pandas as pd

from common import (DATA, REPO, RUN_DIRS, RUNS, WORK, FortranBk, NativeBk, clean_log_lines, load_item, parse_fortran_real,
                    printed_half_unit, save_item)
from item13_propagation import averages_vs_fortran, bk_compare
from item18b_replay import REPLAY, decompose, files_identical, load, s_diff, summary_diff
from parse_logs import averages_out, cpp_surface_fixed, energy_blocks
from provenance import cmake_cache, cpp_timing, mtime, script_header, sh, time_v

STRICT = RUNS / "thesis_fig5_3_dz_cuda"
STRICT_LOG = RUNS / "thesis_fig5_3_dz_cuda.log"
STRICT_CLEAN = WORK / "logs" / "gpu_strict.txt"
BUILD = REPO / "build-cuda"
INPUT = "inputs/thesis_fig5_3_dz/input.json"
REPLAY_FAST = REPLAY / "build-fast"
STRICT_RP_TSV = DATA / "strict_reaction_probability.tsv"
STRICT_BK_TSV = DATA / "strict_bk_error_vs_time.tsv"
NUMBER = re.compile(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?")


def numeric_diff(a, b):
    """Largest absolute and relative difference of the floating-point numbers of two text files."""
    x = np.array([float(v) for v in NUMBER.findall(a.read_text())])
    y = np.array([float(v) for v in NUMBER.findall(b.read_text())])
    assert x.size == y.size
    d = np.abs(x - y)
    return {"identical": bool(a.read_bytes() == b.read_bytes()), "n_numbers": int(x.size),
            "max_abs": float(d.max()), "max_abs_over_max_value": float(d.max() / np.abs(y).max())}


def echo_comparison(a, b):
    """INPUT echoes of two runs, line by line: identical; numbers only (same text, numbers differ:
    derived values such as E_trans, alpha, k_Z, the theta nodes and weights), split at 1e-13 relative
    (the set-up criterion of criteria.md); diagnostic (the PES symmetry self-check lines, a maximum
    over operations); or text differing (an input difference)."""
    out = {"n_lines": len(a), "identical": 0, "numbers_within_1e-13": 0, "numbers_above_1e-13": [], "numbers_max_rel": 0.0,
           "diagnostic": [], "text_differing": []}
    assert len(a) == len(b)
    for la, lb in zip(a, b):
        if la == lb:
            out["identical"] += 1
        elif "Self-check" in la:
            out["diagnostic"].append([la.strip(), lb.strip()])
        elif NUMBER.sub("#", la) == NUMBER.sub("#", lb):
            x = np.array([float(v) for v in NUMBER.findall(la)])
            y = np.array([float(v) for v in NUMBER.findall(lb)])
            rel = float(np.max(np.abs(x - y) / np.maximum(np.abs(y), 1e-300)))
            out["numbers_max_rel"] = max(out["numbers_max_rel"], rel)
            if rel <= 1e-13:
                out["numbers_within_1e-13"] += 1
            else:
                out["numbers_above_1e-13"].append({"strict": la.strip(), "fast": lb.strip(), "rel": rel})
        else:
            out["text_differing"].append([la.strip(), lb.strip()])
    return out


def input_echo(lines):
    """The INPUT section of a C++ log (every input value as the run read it)."""
    i = lines.index("  INPUT")
    j = lines.index("  MEMORY ESTIMATE")
    return lines[i:j]


def bk_header_diff():
    """Stored restart parameters and channel energies of the strict vs the fast GPU BK file."""
    a, b = NativeBk(STRICT / "bk_coefficients.bin"), NativeBk(RUN_DIRS["gpu"] / "bk_coefficients.bin")
    diff = {k: {"strict": a.parameters[k], "fast": b.parameters[k], "rel": abs(a.parameters[k] - b.parameters[k]) / abs(b.parameters[k])}
            for k in a.parameters if a.parameters[k] != b.parameters[k]}
    d = np.abs(a.energies - b.energies)
    return {"n_parameters": len(a.parameters), "differing_parameters": diff,
            "channel_energies_max_rel": float(np.max(d / np.abs(b.energies)))}


def validity(lines):
    built = mtime(BUILD / "spodvr")
    started = dt.datetime.fromisoformat(load_item("phase0")["runs"]["gpu"]["start"].replace(" ", "T"))
    input_changed = sh("git", "log", "-1", "--format=%cI", "--", INPUT)
    return {"script": script_header(STRICT_LOG), "time": time_v(lines), "program_timing": cpp_timing(lines),
            "build": {"dir": "build-cuda", "cmake": cmake_cache(BUILD), "spodvr_mtime": built},
            "source_commits_after_build": sh("git", "log", "--oneline", f"--since={built}", "--", "src", "CMakeLists.txt",
                                             "cmake"),
            "input_last_commit": input_changed,
            "input_unchanged_since_fast_runs": bool(dt.datetime.fromisoformat(input_changed) < started.astimezone()),
            "input_working_tree_clean": sh("git", "status", "--porcelain", "--", INPUT) == "",
            "input_echo_vs_fast_gpu": echo_comparison(input_echo(lines), input_echo(clean_log_lines("gpu"))),
            "analysis_line": dict(zip(("Z", "t"), cpp_surface_fixed(lines))),
            "bk_header_vs_fast_gpu": bk_header_diff(),
            "setup_files_vs_fast_gpu": {f: numeric_diff(STRICT / f, RUN_DIRS["gpu"] / f)
                                        for f in ("asymptotic_potential.dat", "gas_phase_states.dat")}}


def averages():
    f = averages_out(RUN_DIRS["fortran"] / "averages.out")
    rows = [l.split() for l in (STRICT / "averages.dat").read_text().splitlines() if not l.startswith("#")]
    fast = [l.split() for l in (RUN_DIRS["gpu"] / "averages.dat").read_text().splitlines() if not l.startswith("#")]
    a = np.array(rows, dtype=float)
    b = np.array(fast, dtype=float)
    return {"strict_vs_fortran": averages_vs_fortran(f, rows),
            "strict_vs_fast_gpu_norm_rel_max": float(np.max(np.abs(a[:, 1] - b[:, 1]) / np.abs(b[:, 1]))),
            "final_norm": {"strict": float(a[-1, 1]), "fast_gpu": float(b[-1, 1])}}


def bk_summary(r, key):
    """The numbers of one pair of an item-13c result that the strict-run comparison reports."""
    p = r[key]
    amp = r["amplitude_vs_time"][key]
    fam = r["specular_inelastic_family"][key]["per_J"]
    return {"identical": p["identical"], "max_over_scale": p["max_over_scale"],
            "per_channel_relative": {k: p["per_channel_relative"][k] for k in ("median", "p95", "max")},
            "rel_l2_in_window_median": p["l2"]["rel_in_window"]["median"],
            "amplitude_window_median": amp["alpha_window_median"],
            "amplitude_window_quartile_medians": amp["window_quartile_medians"],
            "amplitude_weighted_slope_per_au": amp["fit_weighted_slope_per_au"],
            "family_rel_by_J": {j: v["rel"] for j, v in fam.items()}}


def bk_part(t_fix):
    F = FortranBk(RUN_DIRS["fortran"] / "fort.21")
    readers = {"cpu": NativeBk(RUN_DIRS["cpu"] / "bk_coefficients.bin"),
               "gpu": NativeBk(RUN_DIRS["gpu"] / "bk_coefficients.bin"),
               "strict": NativeBk(STRICT / "bk_coefficients.bin")}
    t = readers["strict"].times()
    assert all(r.n_records == F.n_records for r in readers.values())
    assert all(np.array_equal(r.times(), t) for r in readers.values())
    c = {k: np.stack([r.record(i) for i in range(r.n_records)]) for k, r in readers.items()}
    base = {"fortran": np.stack([F.record(i) for i in range(F.n_records)]), "times": t,
            "labels": readers["strict"].labels, "nstata": readers["strict"].nstata,
            "times_fortran_minus_cpp_max": float(np.max(np.abs(F.times() - t)))}
    # bk_compare's pairs are (cpu, fortran), (gpu, fortran), (gpu, cpu): the strict run takes the
    # place of "gpu", and "cpu" is the fast CPU run (first call) or the fast GPU run (second call)
    s_cpu, s_gpu = {}, {}
    with_cpu = bk_compare(base | {"cpu": c["cpu"], "gpu": c["strict"]}, t_fix, tsv_path=None, series_out=s_cpu)
    with_gpu = bk_compare(base | {"cpu": c["gpu"], "gpu": c["strict"]}, t_fix, tsv_path=None, series_out=s_gpu)
    write_bk_series(s_cpu, s_gpu, with_cpu["scale"])
    return {"strict-fortran": bk_summary(with_cpu, "gpu-fortran"),
            "strict-cpu_fast": bk_summary(with_cpu, "gpu-cpu"),
            "strict-gpu_fast": bk_summary(with_gpu, "gpu-cpu"),
            "cpu_fast-fortran": bk_summary(with_cpu, "cpu-fortran"),     # check: equals item 13c
            "gpu_fast-fortran": bk_summary(with_gpu, "cpu-fortran")}


def write_bk_series(s_cpu, s_gpu, scale):
    """Figure 3 time series of the strict run, in the format of data/bk_error_vs_time.tsv."""
    cols = {"t": s_cpu["t"], "norm_c_t": s_cpu["norm_c_t"],
            "strict-fortran": s_cpu["gpu-fortran"], "strict-fortran_alpha": s_cpu["gpu-fortran_alpha"],
            "strict-gpu_fast": s_gpu["gpu-cpu"], "strict-gpu_fast_alpha": s_gpu["gpu-cpu_alpha"]}
    lines = ["# Strict-build GPU run: BK coefficient differences per analysis step, max over channels |c_a - c_b| / scale "
             f"(scale = max(max|c_Fortran|, 1e-4) = {scale:.6e}); *_alpha: amplitude part (item 13c)", "\t".join(cols)]
    for i in range(len(cols["t"])):
        lines.append("\t".join(f"{v[i]:.6e}" for v in cols.values()))
    STRICT_BK_TSV.write_text("\n".join(lines) + "\n")


def write_reaction_tsv():
    """Per-energy P_reaction of the strict run (Figures 1 and 2): the run's own table, the replay of its
    BK file and its propagation-only difference from the Fortran (same build-fast helper as
    data/replay_reaction_probability.tsv), and the largest |dP(v, j, mJ)| per energy."""
    run = load(STRICT, "scattering_summary.tsv")
    rs, rf = (load(REPLAY_FAST / k, "scattering_summary.tsv") for k in ("gpu_strict", "fortran"))
    vs, vf = (load(REPLAY_FAST / k, "probabilities_rovib.tsv") for k in ("gpu_strict", "fortran"))
    vjm = pd.Series(np.abs(vs["P"].to_numpy() - vf["P"].to_numpy())).groupby(vf["E_coll"].to_numpy()).max().to_numpy()
    E, P_run, P_rs, P_rf = (x.to_numpy() for x in (run["E_coll"], run["P_reaction"], rs["P_reaction"], rf["P_reaction"]))
    assert np.allclose(E, rf["E_coll"].to_numpy(), rtol=0, atol=1e-12) and vjm.size == E.size   # grids agree to 1 ulp
    assert all(np.array_equal(vs[k].to_numpy(), vf[k].to_numpy()) for k in ("v", "j", "mJ"))
    lines = ["# Strict-build GPU run: P_reaction of the run, of the replay of its BK file and of the Fortran's BK file "
             "(helper `replay`, build-fast); abs_strict_minus_F = |replay strict - replay Fortran| (propagation only).",
             "E_coll\tP_run_strict\tP_replay_strict\tP_replay_fortran\tabs_strict_minus_F\tmax_dPvjm_strict_F"]
    for i in range(E.size):
        lines.append(f"{E[i]:.4f}\t{P_run[i]:.16e}\t{P_rs[i]:.16e}\t{P_rf[i]:.16e}\t{abs(P_rs[i] - P_rf[i]):.3e}\t{vjm[i]:.3e}")
    STRICT_RP_TSV.write_text("\n".join(lines) + "\n")


def decomposition_summary(d):
    fam = d["specular_inelastic_family"]
    return {"kappa_median": d["alpha_median"], "kappa_min": d["alpha_min"], "kappa_max": d["alpha_max"],
            "beta_abs_max": d["beta_abs_max"], "kappa_corr_inv_sqrtE": d["kappa_corr_inv_sqrtE"],
            "rel_dS_large_channels_median": d["rel_dS_large_channels_median"],
            "residual_rel_normweighted_outside_family": fam["residual_rel_normweighted_outside_family"],
            "family_rel_normweighted_by_J": {j: v["rel_residual_normweighted"] for j, v in fam["per_J"].items()},
            "family_per_element_max_by_J": {j: v["max"] for j, v in fam["per_element_relative"].items()},
            "residual_by_size": fam["residual_by_size"]}


def replay_part():
    names = {"strict": REPLAY_FAST / "gpu_strict", "fortran": REPLAY_FAST / "fortran", "cpu_fast": REPLAY_FAST / "cpu",
             "gpu_fast": REPLAY_FAST / "gpu"}
    sm = {k: load(p, "s_matrix.tsv") for k, p in names.items()}
    su = {k: load(p, "scattering_summary.tsv") for k, p in names.items()}
    out = {}
    for ref in ("fortran", "cpu_fast", "gpu_fast"):
        sa = sm["strict"]["Re_S"].to_numpy() + 1j * sm["strict"]["Im_S"].to_numpy()
        sb = sm[ref]["Re_S"].to_numpy() + 1j * sm[ref]["Im_S"].to_numpy()
        out[f"strict-{ref}"] = {"S_abs_max": s_diff(sm["strict"], sm[ref])["S_abs"]["max"],
                                "P_reaction": summary_diff(su["strict"], su[ref])["P_reaction"],
                                "S_relative_normweighted": float(np.sqrt(np.sum(np.abs(sa - sb) ** 2) / np.sum(np.abs(sb) ** 2))),
                                "decomposition": decomposition_summary(decompose(sm["strict"], sm[ref]))}
    # the fast runs with the same summary, for the side-by-side table
    for a in ("cpu_fast", "gpu_fast"):
        out[f"{a}-fortran"] = {"decomposition": decomposition_summary(decompose(sm[a], sm["fortran"]))}
    return out


def self_check():
    own = REPLAY / "build-cuda" / "gpu_strict"
    res = {"bitwise": files_identical(STRICT, own)}
    if not res["bitwise"]:
        res["S_abs_max"] = s_diff(load(STRICT, "s_matrix.tsv"), load(own, "s_matrix.tsv"))["S_abs"]["max"]
    return res


def end_to_end():
    blocks = energy_blocks(clean_log_lines("fortran"))
    run = load(STRICT, "scattering_summary.tsv")["P_reaction"].to_numpy()
    texts = [b["p_reaction"] for b in blocks]
    ref = np.array([parse_fortran_real(x) for x in texts])
    half = np.array([printed_half_unit(x) for x in texts])
    replay_f = load(REPLAY_FAST / "fortran", "scattering_summary.tsv")["P_reaction"].to_numpy()
    d = np.abs(run - ref)
    within = d <= np.maximum(half, 1e-11)
    # item 14's rule: outside its half unit, a value is explained when the Fortran's own BK coefficients
    # through the C++ analysis reproduce the printed value and the run is within 1e-10 of that replay
    explained = ~within & (np.abs(replay_f - ref) <= half) & (np.abs(run - replay_f) <= 1e-10)
    return {"n": int(d.size), "n_within_half_unit": int(np.sum(d <= half)),
            "n_within_half_unit_or_floor": int(within.sum()), "n_explained": int(explained.sum()),
            "n_outside_unexplained": int(np.sum(~within & ~explained)),
            "max_abs_vs_printed": float(d.max()),
            "full_precision_vs_replay_fortran_max": float(np.max(np.abs(run - replay_f)))}


def verdicts(res):
    """The pre-registered rules of criteria.md (2026-11-08) applied to the measured values."""
    kappa = abs(res["replay"]["strict-fortran"]["decomposition"]["kappa_median"])
    fam = res["replay"]["strict-fortran"]["decomposition"]["family_rel_normweighted_by_J"][10]
    v = res["validity"]
    return {"valid": bool(v["script"]["exit_code"] == 0 and v["program_timing"].get("steps") == 7000
                          and v["build"]["cmake"].get("SPODVR_FAST_MATH") == "OFF" and v["source_commits_after_build"] == ""
                          and v["input_unchanged_since_fast_runs"] and v["input_working_tree_clean"]
                          and not v["input_echo_vs_fast_gpu"]["text_differing"] and v["analysis_line"]["t"] == 2220.0),
            "kappa": {"value": kappa, "verdict": "fast math" if kappa <= 1.5e-12 else "not fast math" if kappa >= 4.5e-12 else "partly"},
            "family_j10": {"value": fam, "verdict": "fast math involved" if fam <= 8.4e-10 else "not fast math" if fam >= 4.2e-9
                           else "partly"},
            "fast_math_effect_S_rel": res["replay"]["strict-gpu_fast"]["S_relative_normweighted"],
            "end_to_end_pass": bool(res["averages"]["strict_vs_fortran"]["n_outside"] == 0
                                    and res["end_to_end"]["n_outside_unexplained"] == 0),
            "max_abs_dP_reaction_vs_fortran": res["replay"]["strict-fortran"]["P_reaction"]["max"],
            "thresholds": {"kappa_fast_math": 1.5e-12, "kappa_not_fast_math": 4.5e-12, "family_fast_math": 8.4e-10,
                           "family_not_fast_math": 4.2e-9, "fast_math_effect_predicted": 1e-12}}


def ratios(res):
    """Ratios quoted in the report: fast CPU run (item 18b) over strict GPU run, against the Fortran."""
    i18 = load_item("item18b")
    rp = res["replay"]
    fast_kappa = rp["cpu_fast-fortran"]["decomposition"]["kappa_median"]
    strict = rp["strict-fortran"]
    fam = lambda k: rp[k]["decomposition"]["family_rel_normweighted_by_J"][10]
    return {"kappa_fast_cpu_over_strict": fast_kappa / strict["decomposition"]["kappa_median"],
            "S_rel_fast_cpu_over_strict": i18["S_relative_normweighted"]["cpu-fortran"] / strict["S_relative_normweighted"],
            "P_reaction_fast_cpu_over_strict": i18["propagation_only"]["cpu-fortran"]["summary"]["P_reaction"]["max"]
            / strict["P_reaction"]["max"],
            "family_j10_fast_cpu_over_strict": fam("cpu_fast-fortran") / fam("strict-fortran"),
            "family_j10_strict_over_outside": fam("strict-fortran")
            / strict["decomposition"]["residual_rel_normweighted_outside_family"],
            "S_rel_strict_F_over_fast_gpu_cpu": strict["S_relative_normweighted"] / i18["S_relative_normweighted"]["gpu-cpu"],
            "fast_math_effect_over_prediction": rp["strict-gpu_fast"]["S_relative_normweighted"] / 1e-12}


def main():
    lines = STRICT_CLEAN.read_text(encoding="utf-8").splitlines()
    res = {"validity": validity(lines), "averages": averages()}
    res["bk"] = bk_part(res["validity"]["analysis_line"]["t"])
    res["replay"] = replay_part()
    write_reaction_tsv()
    res["self_check"] = self_check()
    res["end_to_end"] = end_to_end()
    fast_wall = load_item("phase0")["runs"]["gpu"]["time"]["wall_s"]
    res["timing"] = {"wall_s": res["validity"]["time"]["wall_s"], "fast_gpu_wall_s": fast_wall,
                     "strict_over_fast": res["validity"]["time"]["wall_s"] / fast_wall}
    res["verdicts"] = verdicts(res)
    res["ratios"] = ratios(res)
    save_item("strict_gpu", res)
    v = res["validity"]
    print("exit", v["script"]["exit_code"], v["time"]["exit_status"], "steps", v["program_timing"].get("steps"),
          "FAST_MATH", v["build"]["cmake"].get("SPODVR_FAST_MATH"), "src commits after build:", repr(v["source_commits_after_build"]),
          "input unchanged:", v["input_unchanged_since_fast_runs"], v["input_working_tree_clean"],
          "line", v["analysis_line"])
    print("INPUT echo vs fast GPU:", v["input_echo_vs_fast_gpu"])
    print("BK header vs fast GPU:", v["bk_header_vs_fast_gpu"])
    print("set-up vs fast GPU:", v["setup_files_vs_fast_gpu"])
    print("averages:", res["averages"]["strict_vs_fortran"]["n_outside"], "outside of", res["averages"]["strict_vs_fortran"]["n_values"],
          "worst", round(res["averages"]["strict_vs_fortran"]["worst"]["ratio"], 3))
    for k, b in res["bk"].items():
        print(f"BK {k}: max/scale {b['max_over_scale']:.2e} rel l2 {b['rel_l2_in_window_median']:.2e} "
              f"amplitude {b['amplitude_window_median']:.3e} family {({j: f'{x:.1e}' for j, x in b['family_rel_by_J'].items()})}")
    for k, r in res["replay"].items():
        d = r["decomposition"]
        extra = (f"|dS| {r['S_abs_max']:.2e} dP_reaction {r['P_reaction']['max']:.2e} S_rel {r['S_relative_normweighted']:.2e} "
                 if "S_abs_max" in r else "")
        print(f"replay {k}: {extra}kappa {d['kappa_median']:.3e} [{d['kappa_min']:.2e}, {d['kappa_max']:.2e}] outside family "
              f"{d['residual_rel_normweighted_outside_family']:.2e} family "
              f"{({j: f'{x:.1e}' for j, x in d['family_rel_normweighted_by_J'].items()})}")
    print("self-check:", res["self_check"])
    print("end to end:", res["end_to_end"])
    print("timing:", res["timing"])
    print("verdicts:", res["verdicts"])
    print("ratios:", {k: round(v, 3) for k, v in res["ratios"].items()})


if __name__ == "__main__":
    main()

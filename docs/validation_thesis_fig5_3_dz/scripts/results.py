#!/usr/bin/env python3
"""Merge the per-item extracts (data/items/*.json) and the run input into scripts/results.json,
the single source of the numbers in the report, plus a few derived values used in the prose."""
from __future__ import annotations

import json
import math

from common import INPUT_DIR, ITEMS, RESULTS, to_jsonable


def main():
    res = {p.stem: json.loads(p.read_text()) for p in sorted(ITEMS.glob("*.json"))}
    res["input"] = json.loads((INPUT_DIR / "input.json").read_text())
    p0 = res["phase0"]["runs"]
    wall = {k: p0[k]["time"]["wall_s"] for k in p0}
    i14 = res["item14_16"]
    rb = res["item18b"]["propagation_only"]
    res["derived"] = {
        "speedup_cpu_vs_fortran": wall["fortran"] / wall["cpu"],
        "speedup_gpu_vs_fortran": wall["fortran"] / wall["gpu"],
        "speedup_gpu_vs_cpu": wall["cpu"] / wall["gpu"],
        "wall_hours": {k: v / 3600 for k, v in wall.items()},
        "propagation_share_min": min(
            [p0["fortran"]["program_timing"]["Propagating the wavefunction"] / wall["fortran"],
             p0["cpu"]["program_timing"]["propagation_s"] / wall["cpu"], p0["gpu"]["program_timing"]["propagation_s"] / wall["gpu"]]),
        # with the strict GPU run (added 2026-11-08), rounded down to 0.1% for "≥ x%" (Figure 4)
        "propagation_share_min_with_strict_floor": math.floor(1000 * min(
            [p0["fortran"]["program_timing"]["Propagating the wavefunction"] / wall["fortran"],
             p0["cpu"]["program_timing"]["propagation_s"] / wall["cpu"], p0["gpu"]["program_timing"]["propagation_s"] / wall["gpu"],
             res["strict_gpu"]["validity"]["program_timing"]["propagation_s"] / res["strict_gpu"]["validity"]["time"]["wall_s"]])) / 1000,
        "p_reaction_full_precision_max": {
            "cpu-fortran": rb["cpu-fortran"]["summary"]["P_reaction"]["max"],
            "gpu-fortran": rb["gpu-fortran"]["summary"]["P_reaction"]["max"],
            "gpu-cpu": i14["item14"]["gpu_vs_cpu"]["P_reaction"]["max"]},
        "p_reaction_full_precision_median": {
            "cpu-fortran": rb["cpu-fortran"]["summary"]["P_reaction"]["median"],
            "gpu-fortran": rb["gpu-fortran"]["summary"]["P_reaction"]["median"],
            "gpu-cpu": i14["item14"]["gpu_vs_cpu"]["P_reaction"]["median"]},
        "amplitude_per_kinetic_block_measured": res["item13"]["bk"]["amplitude_vs_time"]["cpu-fortran"]["fit_slope_per_au"] * 20.0 / 9.0,
        "norm_per_kinetic_block_measured": 2 * res["item13"]["bk"]["amplitude_vs_time"]["cpu-fortran"]["fit_slope_per_au"] * 20.0 / 9.0,
        "amplitude_per_kinetic_block_fft_scales": res["item13_norm_bias"]["predicted_amplitude_ratio_minus_1"]["scattering_grid"]
        / res["item13_norm_bias"]["n_kinetic_blocks"],
        "fft_scale_share_of_measured_rate": (res["item13_norm_bias"]["predicted_amplitude_ratio_minus_1"]["scattering_grid"]
                                             / res["item13_norm_bias"]["n_kinetic_blocks"])
        / (res["item13"]["bk"]["amplitude_vs_time"]["cpu-fortran"]["fit_slope_per_au"] * 20.0 / 9.0),
        "explained_printed_counts_cpu": {k: v["n_explained"] for k, v in i14["item14"]["cpu_vs_fortran"].items()}
        | {"P_order": i14["item15"]["cpu_vs_fortran"]["n_explained"],
           "Re S": i14["item16"]["cpu_vs_fortran"]["Re S"]["n_explained"],
           "Im S": i14["item16"]["cpu_vs_fortran"]["Im S"]["n_explained"],
           "|S|^2": i14["item16"]["cpu_vs_fortran"]["|S|^2"]["n_explained"]},
        "n_rotational_functions": int([r for r in res["setup"]["rows"] if r["quantity"] == "rotation channels in basis"][0]["fortran"]),
        "n_scattering_grid_values": int([r for r in res["setup"]["rows"] if r["quantity"] == "rotation channels in basis"][0]["fortran"])
        * res["input"]["nZ"] * res["input"]["nR"] * res["input"]["nX"] * res["input"]["nY"],
    }
    fam = res["item13"]["bk"]["specular_inelastic_family"]
    ratios = [fam["cpu-fortran"]["per_J"][j]["max_abs_dc"] / fam["gpu-cpu"]["per_J"][j]["max_abs_dc"] for j in fam["cpu-fortran"]["per_J"]]
    res["derived"]["family_cpuF_over_gpuC_min"] = min(ratios)
    res["derived"]["family_cpuF_over_gpuC_max"] = max(ratios)
    res["derived"]["gate4_margin_min"] = 1e-10 / max(res["item13"]["bk"][k]["max_over_scale"] for k in ("cpu-fortran", "gpu-fortran", "gpu-cpu"))
    dcf = res["item18b"]["amplitude_decomposition"]["cpu-fortran"]["specular_inelastic_family"]["per_J"]
    dgc = res["item18b"]["amplitude_decomposition"]["gpu-cpu"]["specular_inelastic_family"]["per_J"]
    sr = [dcf[j]["rel_residual_normweighted"] / dgc[j]["rel_residual_normweighted"] for j in dcf]
    res["derived"]["family_S_cpuF_over_gpuC_min"] = min(sr)
    res["derived"]["family_S_cpuF_over_gpuC_max"] = max(sr)
    sn = res["item18b"]["S_relative_normweighted"]
    res["derived"]["S_rel_cpuF_over_gpuC"] = sn["cpu-fortran"] / sn["gpu-cpu"]
    res["derived"]["explained_S_components_cpu"] = (i14["item16"]["cpu_vs_fortran"]["Re S"]["n_explained"]
                                                    + i14["item16"]["cpu_vs_fortran"]["Im S"]["n_explained"])
    res["derived"]["explained_printed_total_cpu"] = sum(res["derived"]["explained_printed_counts_cpu"].values())
    RESULTS.write_text(json.dumps(to_jsonable(res), indent=1, sort_keys=True) + "\n")
    print(f"results.json: {len(res)} sections, {RESULTS.stat().st_size // 1024} KiB")


if __name__ == "__main__":
    main()

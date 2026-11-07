#!/usr/bin/env python3
"""Figure 4: execution time of the three runs as run (horizontal bars, set-up + rest vs
propagation stacked, speed-up against the Fortran at the bar ends), and the timing table
fragment tables/timing.md. Reads data/items/phase0.json."""
from __future__ import annotations

from common import FIGURES, TABLES, load_item
from style import PROGRAM, INK_2, MUTED, plt, save, setup


def main():
    setup()
    p0 = load_item("phase0")
    R, H, B = p0["runs"], p0["hardware"], p0["builds"]
    wall = {k: R[k]["time"]["wall_s"] for k in R}
    prop = {"fortran": R["fortran"]["program_timing"]["Propagating the wavefunction"],
            "cpu": R["cpu"]["program_timing"]["propagation_s"], "gpu": R["gpu"]["program_timing"]["propagation_s"]}
    order = ["gpu", "cpu", "fortran"]
    fig, ax = plt.subplots(figsize=(8.8, 3.9))
    for i, k in enumerate(order):
        st = PROGRAM[k]
        ax.barh(i, prop[k] / 3600, color=st["color"], height=0.55, label=None)
        ax.barh(i, (wall[k] - prop[k]) / 3600, left=prop[k] / 3600, color=st["color"], alpha=0.35, height=0.55)
        sp = wall["fortran"] / wall[k]
        txt = f"{wall[k] / 3600:.2f} h" + ("" if k == "fortran" else f"   ×{sp:.1f} vs Fortran")
        ax.text(wall[k] / 3600 + 0.2, i, txt, va="center", fontsize=9.5, color="#0b0b0b")
    ax.set_yticks(range(3), [PROGRAM[k]["short"] + ("\n(RTX 4070 Ti)" if k == "gpu" else "\n(8 threads)") for k in order])
    share = min(prop[k] / wall[k] for k in wall)
    ax.set_xlabel(f"Wall-clock time (hours); the propagation is ≥ {100 * share:.1f}% of each run (light end: set-up, analysis, output)")
    ax.set_xlim(0, 21.5)
    ax.grid(axis="y", visible=False)
    gpu = H["gpu"].split(",")
    box = (f"CPU  {H['cpu_model']}\n     {H['cores_per_socket']} cores / {H['logical_cpus']} threads, {H['mem_total_gib']:.1f} GiB RAM\n"
           f"GPU  {gpu[0].strip()}, {gpu[1].strip()}\nOS   {H['os']} under WSL2\n"
           f"all runs: 8 OpenMP threads, one after another,\nas run (not like-for-like, see table)")
    ax.text(0.985, 0.04, box, transform=ax.transAxes, ha="right", va="bottom", multialignment="left", family="DejaVu Sans Mono", fontsize=7.8,
            bbox=dict(boxstyle="round,pad=0.4", fc="#f6f6f4", ec=MUTED, lw=0.6))
    ax.set_title(f"Execution time of the 7000-step run (7000 × 2.5 a.u.)", loc="left")
    save(fig, FIGURES / "fig4_execution_time.png")

    # table fragment
    ffts = {"fortran": "FFTW 3.3.10, ESTIMATE plans", "cpu": "oneMKL 2026.2 (FFTW3 interface)", "gpu": "cuFFT 12.3 (CUDA 13.3)"}
    comp = {"fortran": "gfortran 15.2 -O3 -fopenmp -fdefault-real-8", "cpu": "g++ 15.2, fast-math, -march=native, LTO",
            "gpu": "g++ 15.2 + nvcc 13.3, fast-math (FMA kernels), -march=native"}
    lines = ["| Run | Wall | Propagation | s/step | User CPU | CPU busy | Max RSS | Compiler, flags | FFT library | Threads |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for k in ("fortran", "cpu", "gpu"):
        t = R[k]["time"]
        lines.append(f"| {PROGRAM[k]['short']} | {t['wall']} ({wall[k]:.0f} s) | {prop[k]:.0f} s | {prop[k] / 7000:.2f} | "
                     f"{t['user_s']:.0f} s | {t['cpu_percent']:.0f}% | {t['max_rss_gib']:.1f} GiB | {comp[k]} | {ffts[k]} | 8 |")
    TABLES.mkdir(exist_ok=True)
    (TABLES / "timing.md").write_text("\n".join(lines) + "\n")
    print({k: round(wall["fortran"] / wall[k], 2) for k in wall}, {k: round(prop[k] / 7000, 3) for k in prop})


if __name__ == "__main__":
    main()

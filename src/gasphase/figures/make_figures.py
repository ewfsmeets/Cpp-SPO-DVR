"""The gas-phase figure of src/gasphase/README.md: V_asy(r) with its J = 0 levels and the available energy (left),
and the tails of u_00 and u_10 with the sign rule (right), for the README's worked example inputs/hardwall_small.

Data: a potential-only run of inputs/hardwall_small (its asymptotic_potential.dat, gas_phase_states.dat and
run_info.json). E_avail = EMaxWP + (E(VIn, JIn) - E(0, 0)) + EParallel, as the README defines it, from run_info.json.
Writes gas_phase.png next to this script.

Reproduce (from the repository root; the run takes well under a second, the figure about 3 s, about 20 s the
first time matplotlib builds its font cache):

    python3 -I -c "import json; d = json.load(open('inputs/hardwall_small/input.json')); d['Potential'] = True;
        json.dump(d, open('/var/tmp/hardwall_potential.json', 'w'), indent=1)"
    build/spodvr /var/tmp/hardwall_potential.json --threads 8 --output-dir /var/tmp/hardwall_potential
    uv run --locked python -I -B src/gasphase/figures/make_figures.py /var/tmp/hardwall_potential

Usage: python3 -I -B src/gasphase/figures/make_figures.py RUN_DIR   (needs numpy and matplotlib: run it with
`uv run --locked python -I -B ...`)
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent

# Colours and rc settings copied from tools/plots/style.py (INK ... ORANGE, DPI) and the wave-packet designer's
# apply_style (tools/wave_packet_designer/plots.py); copied, not imported, so that this script stands alone.
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e3df"
SURFACE = "#ffffff"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
DPI = 110
WIDTH_IN = 8.0  # 880 px at DPI 110: GitHub's README column

E_SHOWN_EV = 1.6  # top of the left panel's energy axis
R_SHOWN_BOHR = (0.6, 3.6)  # the left panel's r range: the well; the grid goes on to r_last
LOG_FLOOR = 2e-9  # bottom of the right panel's log axis
MARKER_PT = 5.5  # 8.4 px at DPI 110: the smallest marker whose open or filled face reads
EDGE_PT = 0.95  # 1.45 px: the edge of the open markers
STATE_COLOR = {0: BLUE, 1: ORANGE}  # the same state in both panels


def apply_style() -> None:
    """The figure style of the repository's validation figures (white, recessive axes and grid)."""
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "font.size": 9.0,
            "axes.titlesize": 10.0,
            "axes.edgecolor": MUTED,
            "axes.labelcolor": INK,
            "xtick.color": INK_2,
            "ytick.color": INK_2,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "lines.linewidth": 1.6,
            "legend.frameon": False,
            "legend.fontsize": 8.0,
        }
    )


def read_run(run_dir: Path) -> dict:
    """The r grid, V_asy (eV), the states (v, j, E in eV, u) and the run's metadata."""
    info = json.loads((run_dir / "run_info.json").read_text())
    potential = np.loadtxt(run_dir / "asymptotic_potential.dat")
    states_file = run_dir / "gas_phase_states.dat"
    header = next(line for line in states_file.read_text().splitlines() if line.startswith("# r grid:"))
    rows = np.loadtxt(states_file)
    to_ev = info["unit_factors"]["hartree_to_ev"]
    return {
        "info": info,
        "r": np.array([float(x) for x in header.split(":", 1)[1].split()]),
        "v_asy_ev": potential[:, 3],
        "v": rows[:, 0].astype(int),
        "j": rows[:, 1].astype(int),
        "energy_ev": rows[:, 2] * to_ev,
        "u": rows[:, 3:],
    }


def state(run: dict, v: int, j: int) -> int:
    """The row of state (v, j)."""
    return int(np.flatnonzero((run["v"] == v) & (run["j"] == j))[0])


def available_energy_ev(run: dict) -> float:
    """E_avail = EMaxWP + (E(VIn, JIn) - E(0, 0)) + EParallel, in eV."""
    info = run["info"]
    e00 = run["energy_ev"][state(run, 0, 0)]
    e_initial = info["initial_state"]["energy_ev"]
    return info["wave_packet"]["e_max_ev"] + (e_initial - e00) + info["wave_packet"]["e_parallel_ev"]


def turning_points(r: np.ndarray, v: np.ndarray, energy: float) -> tuple[float, float]:
    """The classical turning points of a level: where V(r) = E on either side of the minimum (linear
    interpolation between grid points)."""
    i_min = int(np.argmin(v))
    inner = next(i for i in range(i_min, 0, -1) if v[i - 1] >= energy)
    outer = next(i for i in range(i_min, len(v) - 1) if v[i + 1] >= energy)
    r_in = np.interp(energy, [v[inner], v[inner - 1]], [r[inner], r[inner - 1]])
    r_out = np.interp(energy, [v[outer], v[outer + 1]], [r[outer], r[outer + 1]])
    return float(r_in), float(r_out)


def draw_levels(ax: plt.Axes, run: dict) -> None:
    """Left panel: V_asy(r), the J = 0 levels between their turning points (open solid, closed dotted) and the
    line E(0, 0) + E_avail."""
    r, v_asy = run["r"], run["v_asy_ev"]
    e00 = run["energy_ev"][state(run, 0, 0)]
    e_top = e00 + available_energy_ev(run)
    ax.plot(r, v_asy, color=INK_2, linewidth=1.4)
    ax.annotate("V_asy(r)", (2.05, 1.45), color=INK)
    for v in range(3):
        energy = run["energy_ev"][state(run, v, 0)]
        is_open = energy <= e_top
        r_in, r_out = turning_points(r, v_asy, energy)
        color = STATE_COLOR.get(v, MUTED) if is_open else MUTED
        ax.plot([r_in, r_out], [energy, energy], color=color, linestyle="-" if is_open else ":", linewidth=2.0)
        above = 1000 * (energy - e00)
        label = f"v = {v}: {above:.0f} meV, " + ("open" if is_open else "closed")
        ax.annotate(label, (r_out, energy), xytext=(4, 0), textcoords="offset points", va="center", color=INK)
    ax.axhline(e_top, color=INK, linestyle="-.", linewidth=1.2)
    ax.annotate(
        f"E(0, 0) + E_avail\n(E_avail = {1000 * (e_top - e00):.0f} meV)",
        (2.15, e_top),
        xytext=(0, 3),
        textcoords="offset points",
        ha="left",
        va="bottom",
        color=INK,
    )
    r_start = r[0]
    ax.axvline(r_start, color=MUTED, linestyle="--", linewidth=1.0)
    ax.annotate(f"RStart\n{r_start:.2f}", (r_start, 0.05), xytext=(4, 0), textcoords="offset points", color=INK_2)
    ax.annotate(
        f"r grid ends at {r[-1]:.2f} bohr →\n(V_asy there: {v_asy[-1]:.2f} eV)",
        (R_SHOWN_BOHR[1], 0.05),
        xytext=(-4, 0),
        textcoords="offset points",
        ha="right",
        color=INK_2,
    )
    ax.set_xlim(*R_SHOWN_BOHR)
    ax.set_ylim(0.0, E_SHOWN_EV)
    ax.set_xlabel("r (bohr)")
    ax.set_ylabel("energy (eV, from the potential's zero)")
    ax.set_title("What is open: V_asy(r) and the J = 0 levels")


def draw_tails(ax: plt.Axes, run: dict) -> None:
    """Right panel: |u_00| and |u_10| on a log axis, filled markers where u > 0, open where u < 0; the r absorber,
    the cut inner tail, the alternating plateau, the wrapped tail and the sign rule."""
    r = run["r"]
    for v, style in ((0, "-"), (1, "--")):
        color = STATE_COLOR[v]
        u = run["u"][state(run, v, 0)]
        magnitude = np.where(np.abs(u) > LOG_FLOOR, np.abs(u), np.nan)
        ax.plot(r, magnitude, color=color, linestyle=style, linewidth=1.0, zorder=2, label=f"|u_{v}0(r)|")
        positive, negative = u > 0, u < 0
        ax.plot(r[positive], np.abs(u[positive]), "o", color=color, markersize=MARKER_PT, zorder=3)
        ax.plot(
            r[negative],
            np.abs(u[negative]),
            "o",
            color=color,
            markersize=MARKER_PT,
            markerfacecolor=SURFACE,
            markeredgewidth=EDGE_PT,
            zorder=3,
        )
    absorber = run["info"]["absorbers"]["list"]
    r_opt = next(a["start"] for a in absorber if a["name"] == "ROpt")
    ax.axvspan(r_opt, r[-1], color=GRID, alpha=0.6, linewidth=0)
    ax.annotate("r absorber", (r_opt, 1.5), xytext=(4, 0), textcoords="offset points", color=INK_2, va="center")
    ax.annotate("v = 0", (1.55, 1.9), color=INK, ha="right")
    ax.annotate("v = 1", (2.05, 0.6), color=INK)

    u00 = run["u"][state(run, 0, 0)]
    peak = np.max(np.abs(u00))
    ax.annotate(
        f"inner tail cut at RStart:\n{u00[0]:.1e} (peak {peak:.2f})",
        (r[0], u00[0]),
        xytext=(3.0, 3e-2),
        textcoords="data",
        color=INK,
        arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.8},
    )
    middle = (r > 4.6) & (r < 5.1)
    level = float(np.median(np.abs(u00[middle])))
    late = float(np.abs(u00[np.argmin(np.abs(r - 8.45))]))
    ax.annotate(
        f"plateau: sign alternates\nat every point;\n|u| ≈ {level:.1e} near 5 bohr,\nrising to ~{late:.0e} at 8.5",
        (5.0, level),
        xytext=(4.4, 3e-6),
        textcoords="data",
        ha="left",
        va="bottom",
        color=INK,
        arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.8},
    )
    ax.annotate(
        f"wrapped tail; sign rule:\nu(r_last) > 0 ({u00[-1]:.1e})",
        (r[-1], u00[-1]),
        xytext=(5.0, 4e-3),
        textcoords="data",
        ha="left",
        color=INK,
        arrowprops={"arrowstyle": "-", "color": MUTED, "linewidth": 0.8},
    )
    ax.set_yscale("log")
    ax.set_ylim(LOG_FLOOR, 5.0)
    ax.set_xlim(r[0] - 0.2, r[-1] + 0.2)
    ax.set_xlabel("r (bohr)")
    ax.set_ylabel("|u(r)| (bohr$^{-1/2}$)")
    ax.set_title("The tails and the sign rule: J = 0")
    handles = [
        Line2D([], [], color=BLUE, linestyle="-", label="v = 0"),
        Line2D([], [], color=ORANGE, linestyle="--", label="v = 1"),
        Line2D([], [], color=INK_2, marker="o", linestyle="", markersize=MARKER_PT, label="u > 0"),
        Line2D(
            [],
            [],
            color=INK_2,
            marker="o",
            linestyle="",
            markersize=MARKER_PT,
            markerfacecolor=SURFACE,
            markeredgewidth=EDGE_PT,
            label="u < 0",
        ),
    ]
    ax.legend(handles=handles, loc="lower left", ncols=2, handlelength=1.8, columnspacing=1.0)


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("usage: python3 -I -B src/gasphase/figures/make_figures.py RUN_DIR")
    run = read_run(Path(sys.argv[1]))
    apply_style()
    figure, (left, right) = plt.subplots(1, 2, figsize=(WIDTH_IN, 4.2), layout="constrained")
    draw_levels(left, run)
    draw_tails(right, run)
    figure.savefig(HERE / "gas_phase.png", dpi=DPI)


if __name__ == "__main__":
    main()

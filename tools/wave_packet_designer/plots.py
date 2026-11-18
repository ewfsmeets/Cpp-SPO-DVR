"""The figure of an Evaluation: one function per panel, each drawing on a given matplotlib Axes.

Presentation layer (README, "Architecture"): the panels draw what evaluate() computed and compute nothing.
Used by the window (gui_mpl.py) and by `wave_packet_designer.py check --png`; explanations.py says what each
panel shows. Style as the repository's validation figures (docs/validation_thesis_fig5_3_dz/scripts/style.py):
the first three slots of the dataviz reference palette, which pass its colour-vision checks (validated
2026-11-18). Colour follows the grid (blue the scattering grid, aqua the specular grid, orange r); line styles
and labels separate everything else, so identity never rests on colour alone. One quantity per axes; values
below a panel's floor are not drawn (masked), so nothing lies along the bottom axis that looks like data.
"""

import textwrap
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import NullFormatter

from evaluation import SYMMETRY_GIVEN, SYMMETRY_NOT_GIVEN, Evaluation
from report import format_bytes, status_line, verdict_lines

INK, INK_2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df", "#ffffff"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GRID_COLOR = {"Z scattering": BLUE, "Z specular": AQUA, "r": ORANGE}
LOW_Z = ("ZOptLow", "ZSpOptLow")  # drawn hatched, to tell them from the main absorber on the same grid
STATE_STYLES = ["-", "--", "-.", ":"]  # analysed v; each curve is also labelled with its v
README_LEVELS = {"ZOpt": 4e-5, "ZSpOpt": 4e-9}  # R + T of the H2 examples, inputs/README.md section 7
LANE_HEIGHT = 0.07  # axes fraction of each grid lane at the top of the Z panel
VERDICT_WIDTH, VERDICT_LINES = 54, 26  # characters per line and lines of the verdict panel
Panel = tuple[Axes, Callable[[Axes, Evaluation], None]]


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
            "lines.linewidth": 2.0,
            "legend.frameon": False,
            "legend.fontsize": 8.0,
        }
    )


def report_figure(e: Evaluation) -> Figure:
    """A new figure with every panel of the evaluation (`check --png`)."""
    apply_style()
    figure = plt.figure(figsize=(18.0, 8.5), layout="constrained")
    draw_panels(add_panels(figure, figure.add_gridspec(2, 4)), e)
    figure.suptitle(status_line(e), color=INK, fontsize=10.5)
    return figure


def add_panels(figure: Figure, grid: GridSpec) -> list[Panel]:
    """The panel axes on a 2 x 4 grid of the figure, each with the function that draws it."""
    return [
        (figure.add_subplot(grid[0, :2]), draw_z_layout),
        (figure.add_subplot(grid[0, 2]), draw_momentum),
        (figure.add_subplot(grid[0, 3]), draw_time),
        (figure.add_subplot(grid[1, 0]), draw_absorber_response),
        (figure.add_subplot(grid[1, 1]), draw_r_layout),
        (figure.add_subplot(grid[1, 2]), draw_memory),
        (figure.add_subplot(grid[1, 3]), draw_verdict),
    ]


def draw_panels(panels: list[Panel], e: Evaluation) -> None:
    """Clears each panel and draws the evaluation on it (the window redraws so after every change)."""
    for ax, draw in panels:
        ax.clear()
        draw(ax, e)


def masked(values: np.ndarray, floor: float) -> np.ndarray:
    """The values with those below `floor` set to NaN: not drawn."""
    return np.where(values >= floor, values, np.nan)


def mark(ax: Axes, x: float, label: str, height: float, style: str = "-") -> None:
    """A vertical line across the axes with its label at `height` (axes fraction)."""
    ax.axvline(x, color=INK_2, ls=style, lw=1.1)
    ax.text(x, height, f" {label}", transform=ax.get_xaxis_transform(), color=INK_2, fontsize=8, va="top")


def draw_z_layout(ax: Axes, e: Evaluation) -> None:
    """The Z design: one lane per grid at the top (the scattering grid with its interaction region and
    absorber, the specular grid with the packet's allowed region and its absorber), and below them |psi(Z)|^2
    of the initial packet, which lives on the specular grid only."""
    ax.set_title("Z: the two grids, and the initial wave packet on the specular grid")
    ax.set_xlabel("Z (bohr)")
    ax.set_ylabel("|psi(Z)|^2 (1/bohr)")
    ax.set_yscale("log")
    g, d = e.grids, e.densities
    if g is None:
        return
    floor = 1e-16
    if d is not None:
        ax.plot(d.z, masked(d.z_density, floor), color=INK, lw=1.5)
    peak = float(np.nanmax(d.z_density)) if d is not None else 1.0
    ax.set_ylim(floor, peak * 1e5)  # room above the packet for the lanes and the labels
    ax.set_yticks([t for t in ax.get_yticks() if floor <= t <= 10.0 * peak])  # no ticks beside the lanes
    draw_lanes(ax, e)
    p, s = e.parameters, e.packet_summary
    mark(ax, g.z_analysis, "Zana", 0.76)
    mark(ax, p.z0, "Z0", 0.68, "--")
    if s is not None:
        mark(ax, p.z0 + 5.2 * s.sigma_z, "outer tail", 0.60, ":")
    ax.set_xlim(g.z_specular.start, g.z_specular.last)


def draw_lanes(ax: Axes, e: Evaluation) -> None:
    """The two grids as lanes at the top of the Z panel, each with its regions labelled inside them."""
    g = e.grids
    lane = ax.get_xaxis_transform()  # x in bohr, y in axes fraction
    lanes = (
        ("scattering grid", g.z_scattering, "Z scattering", 1.0),
        ("specular grid", g.z_specular, "Z specular", 0.92),
    )
    for name, grid, kind, top in lanes:
        bottom, middle = top - LANE_HEIGHT, top - LANE_HEIGHT / 2
        ax.fill_between([grid.start, grid.last], bottom, top, transform=lane, color=GRID, lw=0)
        ax.text(grid.start, middle, f" {name}", transform=lane, va="center", fontsize=8, color=INK)
        if kind == "Z scattering":
            ax.fill_between([grid.start, e.z_asymptotic], bottom, top, transform=lane, color=MUTED, alpha=0.3, lw=0)
            ax.text(
                e.z_asymptotic,
                middle,
                "interaction region ",
                transform=lane,
                va="center",
                ha="right",
                fontsize=8,
                color=INK,
            )
        for a in e.absorbers:
            if a.grid == kind:
                ax.fill_between(
                    [a.start, a.end],
                    bottom,
                    top,
                    transform=lane,
                    color=GRID_COLOR[kind],
                    alpha=0.55,
                    lw=0,
                    hatch="//" if a.name in LOW_Z else None,
                )
                ax.text((a.start + a.end) / 2, middle, a.name, transform=lane, va="center", ha="center", fontsize=8)
    if e.packet is not None:
        z = g.z_specular
        low, high = z.point(e.packet.info.i_low), z.point(e.packet.info.i_high)
        ax.annotate(
            "", (low, 0.835), (high, 0.835), xycoords=lane, arrowprops={"arrowstyle": "<->", "color": INK_2, "lw": 1.0}
        )
        ax.text(high, 0.835, " allowed region (TolSp)", transform=lane, va="center", fontsize=8, color=INK_2)


def draw_momentum(ax: Axes, e: Evaluation) -> None:
    """|psi(k)|^2 of the initial packet against the grid's momentum limits +-k_grid = +-pi/dZ, with the
    packet's mean momentum -k0 and the momenta of EMinWP and EMaxWP."""
    ax.set_title("Z momentum: does the grid hold the packet?")
    ax.set_xlabel("k (1/bohr)")
    ax.set_ylabel("|psi(k)|^2 (bohr)")
    ax.set_yscale("log")
    m, s = e.densities, e.packet_summary
    if m is None or s is None:
        return
    floor = 1e-20
    ax.plot(m.k, masked(m.k_density, floor), color=INK, lw=1.5)
    ax.set_ylim(floor, 1e3 * float(m.k_density.max()))
    ax.set_xlim(-1.08 * s.grid_momentum_limit, 1.08 * s.grid_momentum_limit)
    k0, width = e.derived.kz, 2.0 * s.sigma_k
    mark(ax, -s.grid_momentum_limit, "-k_grid", 0.97, ":")
    ax.axvline(s.grid_momentum_limit, color=INK_2, ls=":", lw=1.1)
    ax.text(
        s.grid_momentum_limit,
        0.97,
        "+k_grid ",
        transform=ax.get_xaxis_transform(),
        ha="right",
        va="top",
        fontsize=8,
        color=INK_2,
    )
    axes_x = ax.get_xaxis_transform()
    for k, name, height in ((-k0 - width, "EMaxWP", 0.05), (-k0, "-k0", 0.11), (-k0 + width, "EMinWP", 0.05)):
        ax.plot([k, k], [0.0, 0.04], transform=axes_x, color=INK, lw=1.2)
        ax.text(
            k,
            height,
            name,
            transform=axes_x,
            ha="right" if name == "EMaxWP" else "left" if name == "EMinWP" else "center",
            fontsize=7.5,
            color=INK_2,
        )
    ax.axvspan(0.0, s.grid_momentum_limit, color=ORANGE, alpha=0.08, lw=0)
    ax.text(
        0.5 * s.grid_momentum_limit,
        0.86,
        f"outward (k > 0):\n{m.outward_fraction:.1e} of the norm",
        transform=axes_x,
        ha="center",
        va="top",
        fontsize=8,
        color=INK,
    )


def draw_absorber_response(ax: Axes, e: Evaluation) -> None:
    """Reflection R (solid) and transmission T (dotted) of each absorber over the energies it must absorb,
    with the levels of the README's H2 examples."""
    ax.set_title("Absorbers: what is not absorbed (lower is better)")
    ax.set_xlabel("kinetic energy of the incoming wave (eV)")
    ax.set_ylabel("fraction reflected (R) or transmitted (T)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_minor_formatter(NullFormatter())
    floor = 1e-14
    for r in e.absorber_responses:
        color = GRID_COLOR[next(a.grid for a in e.absorbers if a.name == r.name)]
        marker = "s" if r.name in LOW_Z else "o"
        ax.plot(
            r.energies,
            masked(np.array(r.reflection), floor),
            color=color,
            ls="-",
            marker=marker,
            ms=3.5,
            label=f"{r.name} R",
        )
        ax.plot(
            r.energies,
            masked(np.array(r.transmission), floor),
            color=color,
            ls=":",
            marker=marker,
            ms=3.5,
            mfc=SURFACE,
            label=f"{r.name} T",
        )
        if r.name in README_LEVELS:
            ax.axhline(README_LEVELS[r.name], color=color, lw=0.9, ls="--", alpha=0.7)
            ax.text(
                0.99,
                README_LEVELS[r.name],
                f"README {r.name} R+T ",
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="bottom",
                fontsize=7.5,
                color=INK_2,
            )
    if e.absorber_responses:
        low, high = ax.get_ylim()
        ax.set_ylim(max(low, floor), high * 1e3)  # room for the legend
        ax.legend(loc="upper center", ncols=3, fontsize=7.5)


def draw_r_layout(ax: Axes, e: Evaluation) -> None:
    """For each analysed vibrational level v, the norm beyond r of its most extended analysed bound state: what
    an r absorber starting at r would remove. Against it the limit (1e-6), ROptStart and the r grid."""
    title = "r: would the r absorber remove bound molecules?"
    ax.set_title(title + (" (indicative)" if e.gas_phase_indicative else ""))
    ax.set_xlabel("r (bohr)")
    ax.set_ylabel("norm of the bound state beyond r")
    ax.set_yscale("log")
    floor = 1e-10
    ax.set_ylim(floor, 3.0)
    if e.grids is None or e.gas_phase is None:
        return
    r = e.gas_phase.r
    for v, tail in enumerate(e.bound_state_tails):
        ax.plot(r, masked(tail, floor), color=INK, lw=1.4, ls=STATE_STYLES[v % len(STATE_STYLES)])
        label_at = int(np.argmax(tail < 10.0 ** -(2 + v)))  # each label at its own level, so they do not overlap
        ax.text(r[label_at], tail[label_at], f" v={v}", fontsize=7.5, color=INK, va="bottom")
    if e.gas_reference is not None:  # the converged states: where they part from the program's, dR is too coarse
        for v, tail in enumerate(e.reference_tails):
            ax.plot(e.gas_reference.r, masked(tail, floor), color=MUTED, lw=1.0, ls=STATE_STYLES[v % len(STATE_STYLES)])
        ax.plot([], [], color=INK, lw=1.4, label="on the input's r grid (the program's)")
        ax.plot([], [], color=MUTED, lw=1.0, label="converged (finer r grid)")
        ax.legend(loc="upper right", fontsize=7.5)
    ax.axhline(1e-6, color=ORANGE, lw=1.0, ls="--")
    ax.text(r[0], 1e-6, " limit 1e-6", va="bottom", fontsize=8, color=INK_2)
    r_absorber = next((a for a in e.absorbers if a.name == "ROpt"), None)
    if r_absorber is not None:
        ax.axvspan(r_absorber.start, r_absorber.end, color=ORANGE, alpha=0.15, lw=0)
        mark(ax, r_absorber.start, "ROptStart", 0.97)
    g = e.grids.r
    mark(ax, g.start, "RStart", 0.12, ":")
    ax.axvline(g.last, color=INK_2, ls=":", lw=1.1)
    ax.text(
        g.last, 0.12, "r grid end ", transform=ax.get_xaxis_transform(), ha="right", va="top", fontsize=8, color=INK_2
    )


def draw_memory(ax: Axes, e: Evaluation) -> None:
    """The peak memory of the run on the CPU and on the GPU, against the memory available (if known)."""
    ax.set_title("Memory: does the run fit?")
    ax.set_xlabel("GiB")
    r = e.resources
    if r is None:
        return
    bars = []
    if r.memory is not None:
        bars.append(("CPU", r.memory.peak_bytes, e.settings.ram_bytes, "RAM"))
    if r.gpu_memory is not None:
        bars.append(("GPU", r.gpu_memory.peak_bytes, e.settings.gpu_bytes, "free GPU"))
    sizes = [size for _, size, _, _ in bars] + [cap for _, _, cap, _ in bars if cap is not None]
    for i, (_, size, capacity, capacity_name) in enumerate(bars):
        ax.barh(i, size / 2**30, color=INK_2, height=0.45)
        ax.text(size / 2**30, i, f" {format_bytes(size)}", va="center", fontsize=8.5, color=INK)
        if capacity is not None:
            ax.plot([capacity / 2**30] * 2, [i - 0.35, i + 0.35], color=ORANGE, lw=2.0)
            ax.text(
                capacity / 2**30,
                i - 0.38,
                f"{capacity_name} {format_bytes(capacity)}",
                ha="center",
                va="bottom",
                fontsize=7.5,
                color=INK_2,
            )
    ax.set_yticks(range(len(bars)), [name for name, _, _, _ in bars])
    ax.set_xlim(0.0, 1.3 * max(sizes + [2**30]) / 2**30)
    ax.set_ylim(len(bars) - 0.4, -0.8)
    ax.grid(axis="y", visible=False)
    n_xy = e.parameters.nx * e.parameters.ny
    symmetry = {
        SYMMETRY_GIVEN: f"potential at {r.symmetry.n_representatives} of {n_xy} (x, y) points",
        SYMMETRY_NOT_GIVEN: "UseSymmetry: reduction not given (an upper bound)",
    }.get(r.symmetry_source, "no surface symmetry")
    ax.text(
        0.99, 0.02, f"{symmetry}; {r.threads} threads", transform=ax.transAxes, ha="right", fontsize=7.5, color=INK_2
    )


def draw_time(ax: Axes, e: Evaluation) -> None:
    """Z against time for the packet's fast, central and slow parts: in from Z0, reflected at Z_turn and, if
    elastically, out at the same speed; and an outgoing wave of normal energy E_n,min from the slow part's
    reflection. Against them Zana and tEnd (inputs/README.md section 8)."""
    ax.set_title("Time: back out across Zana by tEnd?")
    ax.set_xlabel("t (a.u.)")
    ax.set_ylabel("Z (bohr)")
    b = e.time_budget
    if b is None:
        return
    for path, style in zip(b.paths, ("-", "--", ":", "-"), strict=True):
        slow_out = path is b.paths[-1]
        ax.plot(
            path.times,
            path.z,
            color=ORANGE if slow_out else INK,
            ls=style,
            lw=1.6 if slow_out else 1.2,
            label=path.label,
        )
    ax.axhline(e.grids.z_analysis, color=INK_2, lw=1.0)
    ax.text(0.0, e.grids.z_analysis, " Zana", va="bottom", fontsize=8, color=INK_2)
    ax.axhline(b.z_turn, color=INK_2, lw=1.0, ls=":")
    ax.text(0.0, b.z_turn, " Z_turn", va="bottom", fontsize=8, color=INK_2)
    mark(ax, b.t_end, "tEnd", 0.5)
    ax.axvline(b.t_end, color=INK, lw=2.0)
    for t, name, side in (
        (b.t_elastic, "elastic", "left"),
        (b.t_needed, f"E_n {b.lowest_outgoing_energy:g} eV", "right"),
    ):
        ax.plot([t], [e.grids.z_analysis], marker="o", color=INK, ms=4)
        ax.text(t, e.grids.z_analysis, f" {name} \n {t:.0f} ", fontsize=7.5, color=INK, va="top", ha=side)
    ax.set_xlim(0.0, 1.1 * max(b.t_end, b.t_needed))
    ax.set_ylim(0.0, 1.7 * e.parameters.z0)  # room above the paths for the legend
    ax.legend(loc="upper right", fontsize=7, ncols=2)


def draw_verdict(ax: Axes, e: Evaluation) -> None:
    """The verdict as text: what spodvr would do, the recommendations missed, the key numbers (report.py)."""
    ax.axis("off")
    ax.set_title("Verdict")
    lines = []
    for paragraph in verdict_lines(e):
        lines += textwrap.wrap(paragraph, VERDICT_WIDTH, subsequent_indent="  ") or [""]
    if len(lines) > VERDICT_LINES:
        lines = lines[: VERDICT_LINES - 1] + ["  ... (the full list: `check`)"]
    ax.text(0.0, 1.0, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=7.0, color=INK)

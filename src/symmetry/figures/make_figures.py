"""The orbit map of src/symmetry/README.md: the C6v orbits of an 8 × 8 grid (nφ = 12) on the hexagonal unit cell.

The 60° cell in Cartesian coordinates (units of the lattice constant) with its lattice vectors, the grid points
coloured by the size of their orbit and labelled with the orbit's letter (the representative in capitals, ringed),
the top, bridge and hollow sites, and the C6v mirror lines through the lattice points of the cell.

A short Python transcription of the orbit walk of src/symmetry/grid_symmetry.cpp (GridSymmetry::build_orbits), with
the operations of point_group.cpp and the lattice matrices of surface_lattice.cpp, for illustration only, not the
implementation. Before drawing, it checks that its orbits reproduce the README's 8 × 8 text map letter for letter and
the counts 10 of 64 and 19 of 144 (12 × 12) quoted there, and stops otherwise. Writes orbits_c6v_8x8.png next to
this script.

Usage: uv run --locked python -I -B src/symmetry/figures/make_figures.py   (run time ≈ 1 s; needs numpy, matplotlib)
"""

import math
import re
import sys
from fractions import Fraction
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
README = HERE.parent / "README.md"
OUTPUT = HERE / "orbits_c6v_8x8.png"

# Colours and settings copied from tools/plots/style.py (INK, INK_2, MUTED, SURFACE, SERIES_COLORS, DPI).
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8985"
SURFACE = "#ffffff"
DPI = 110
# One colour per orbit size, from SERIES_COLORS, in the order that passes the dataviz palette validator (light mode,
# no warnings: lightness, chroma, CVD separation of neighbours, contrast against white).
SIZE_COLORS = {1: "#cc79a7", 3: "#d55e00", 6: "#009e73", 12: "#0072b2"}
SIZE_MEANING = {1: "top site", 3: "bridge sites", 6: "on a mirror line", 12: "general points"}

SKEW = math.pi / 3  # the 60° cell; lattice constant 1
# The hollow sites in fractional coordinates; with the CRP C3v stacking the second layer lies under (1/3, 1/3).
HOLLOWS = {"hcp hollow": ("^", (Fraction(1, 3), Fraction(1, 3))), "fcc hollow": ("v", (Fraction(2, 3), Fraction(2, 3)))}
C6V_GENERATORS = [(False, 2), (True, 0)]  # rotation 60°, mirror at 0° (point_group.cpp); (mirror, k), angle k·30°


# --- the group and its action (surface_operation.cpp, point_group.cpp, surface_lattice.cpp) -----------------------


def compose(a, b):
    """a ∘ b: first b, then a (surface_operation.cpp compose)."""
    sign_a = -1 if a[0] else 1
    return (a[0] != b[0], (a[1] + sign_a * b[1]) % 12)


def group_closure(generators):
    """The group in the canonical order: identity, rotations by increasing k, mirrors by increasing k."""
    group = [(False, 0)]
    i = 0
    while i < len(group):
        for g in generators:
            product = compose(g, group[i])
            if product not in group:
                group.append(product)
        i += 1
    return sorted(group, key=lambda g: (g[0], g[1]))


def inverse(g):
    return g if g[0] else (False, (-g[1]) % 12)


def to_cartesian(u, w):
    """Fractional (u, w) along a₁ = (1, 0) and a₂ = (cos γ, sin γ) → Cartesian."""
    return u + w * math.cos(SKEW), w * math.sin(SKEW)


def to_fractional(x, y):
    w = y / math.sin(SKEW)
    return x - w * math.cos(SKEW), w


def apply_cartesian(g, x, y):
    beta = 2 * math.pi * g[1] / 12
    c, s = math.cos(beta), math.sin(beta)
    y = -y if g[0] else y
    return c * x - s * y, s * x + c * y


def lattice_matrix(g):
    """A_g: the images of a₁ and a₂ in the basis (a₁, a₂), integer for a lattice symmetry."""
    u1, w1 = to_fractional(*apply_cartesian(g, *to_cartesian(1, 0)))
    u2, w2 = to_fractional(*apply_cartesian(g, *to_cartesian(0, 1)))
    entries = [u1, u2, w1, w2]  # a11, a12, a21, a22
    rounded = [round(e) for e in entries]
    assert all(abs(e - r) < 1e-10 for e, r in zip(entries, rounded, strict=True)), f"{g} is not a lattice symmetry"
    return rounded


# --- the orbit walk (grid_symmetry.cpp) -------------------------------------------------------------------------


def usable_operations(nx, ny, nphi):
    """The operations that map the (x, y) and the φ grid onto themselves, with their lattice matrices."""
    usable = []
    for g in group_closure(C6V_GENERATORS):
        a = lattice_matrix(g)
        if (g[1] * nphi) % 12 or (a[1] * nx) % ny or (a[2] * ny) % nx:
            continue
        usable.append((g, a))
    return usable


def orbits(nx, ny, nphi):
    """{(ix, iy): orbit number} and the representatives, walking the points in memory order ixy = iy + ny·ix."""
    usable = usable_operations(nx, ny, nphi)
    orbit_of = {}
    representatives = []
    for ix in range(nx):
        for iy in range(ny):
            if (ix, iy) in orbit_of:
                continue
            representatives.append((ix, iy))
            for _, a in usable:
                image = ((a[0] * ix + (a[1] * nx // ny) * iy) % nx, ((a[2] * ny // nx) * ix + a[3] * iy) % ny)
                orbit_of.setdefault(image, len(representatives) - 1)
    return orbit_of, representatives


def letter_map(nx, ny, nphi):
    """{(ix, iy): letter}, the representative's in capitals, as the README prints it."""
    orbit_of, representatives = orbits(nx, ny, nphi)
    letters = {}
    for point, k in orbit_of.items():
        letter = chr(ord("a") + k)
        letters[point] = letter.upper() if point in representatives else letter
    return letters


# --- the checks against the README -------------------------------------------------------------------------------


def readme_map():
    """{(ix, iy): letter} of the README's text map (the ```text block with the rows iy=…)."""
    blocks = re.findall(r"```text\n(.*?)```", README.read_text(encoding="utf-8"), re.S)
    block = next((b for b in blocks if "iy=0 " in b), None)
    if block is None:
        sys.exit(f"no text map with a row 'iy=0' in {README}")
    letters = {}
    for line in block.splitlines():
        row = re.match(r"iy=(\d+)\s+(.*)$", line)
        if row:
            for ix, letter in enumerate(row.group(2).split()):
                letters[(ix, int(row.group(1)))] = letter
    return letters


def check_against_readme():
    computed, quoted = letter_map(8, 8, 12), readme_map()
    differences = sorted(p for p in set(computed) | set(quoted) if computed.get(p) != quoted.get(p))
    if differences:
        detail = ", ".join(f"(ix, iy) = {p}: README {quoted.get(p)}, orbit walk {computed.get(p)}" for p in differences)
        sys.exit(f"the README's 8 × 8 map differs from the orbit walk at {len(differences)} points: {detail}")
    text = README.read_text(encoding="utf-8")
    for n, expected in ((8, 10), (12, 19)):
        count = len(orbits(n, n, 12)[1])
        if count != expected or f"{expected} of {n * n}" not in text:
            sys.exit(f"{n} × {n}: the orbit walk gives {count} representatives; the README must quote {expected}")


# --- the figure ----------------------------------------------------------------------------------------------------


def clip_to_cell(point, direction):
    """The part of the line point + t·direction inside the cell, in Cartesian coordinates, or None.

    Clipped in fractional coordinates, where the cell is the unit square (the map is linear, so lines stay lines)."""
    p, d = np.array(to_fractional(*point)), np.array(to_fractional(*direction))
    t_low, t_high = -np.inf, np.inf
    for k in range(2):
        if abs(d[k]) < 1e-12:
            if not -1e-12 <= p[k] <= 1 + 1e-12:
                return None
            continue
        t1, t2 = (0 - p[k]) / d[k], (1 - p[k]) / d[k]
        t_low, t_high = max(t_low, min(t1, t2)), min(t_high, max(t1, t2))
    if t_high - t_low < 1e-9:
        return None
    ends = [to_cartesian(*(p + t * d)) for t in (t_low, t_high)]
    return [e[0] for e in ends], [e[1] for e in ends]


def mirror_segments():
    """The C6v mirror lines (at 0°, 30°, …, 150°) through every lattice point that can reach the cell, clipped."""
    segments, seen = [], set()
    for i in range(-1, 3):
        for j in range(-1, 3):
            origin = to_cartesian(i, j)
            for k in range(6):
                angle = k * math.pi / 6
                segment = clip_to_cell(origin, (math.cos(angle), math.sin(angle)))
                if segment is None:
                    continue
                key = tuple(round(v, 6) for end in sorted(zip(*segment, strict=True)) for v in end)
                if key not in seen:
                    seen.add(key)
                    segments.append(segment)
    return segments


def draw():
    letters = letter_map(8, 8, 12)
    orbit_of, _ = orbits(8, 8, 12)
    sizes = {k: list(orbit_of.values()).count(k) for k in set(orbit_of.values())}

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "legend.fontsize": 8, "text.color": INK})
    fig, ax = plt.subplots(figsize=(8.0, 5.0), dpi=DPI)
    fig.patch.set_facecolor(SURFACE)

    for xs, ys in mirror_segments():
        ax.plot(xs, ys, color=MUTED, linewidth=0.6, linestyle=(0, (4, 3)), zorder=1)
    outline = [to_cartesian(u, w) for u, w in ((0, 0), (1, 0), (1, 1), (0, 1), (0, 0))]
    ax.plot(*zip(*outline, strict=True), color=INK, linewidth=0.8, zorder=2)

    for (ix, iy), letter in letters.items():
        x, y = to_cartesian(ix / 8, iy / 8)
        color = SIZE_COLORS[sizes[orbit_of[(ix, iy)]]]
        representative = letter.isupper()
        ax.plot(
            x,
            y,
            marker="o",
            markersize=17 if representative else 13,
            markerfacecolor=color,
            markeredgecolor=INK if representative else color,
            markeredgewidth=2.0 if representative else 0.8,
            alpha=1.0,
            zorder=3,
        )
        ax.plot(x, y, marker="o", markersize=10 if representative else 8, color=SURFACE, zorder=4)
        ax.text(
            x,
            y,
            letter,
            ha="center",
            va="center_baseline",
            fontsize=8,
            fontweight="bold" if representative else "normal",
            color=INK,
            zorder=5,
        )

    # The sites: top and bridge lie on grid points, the hollows at a third of the cell between them.
    for marker, (u, w) in HOLLOWS.values():
        ax.plot(*to_cartesian(u, w), marker=marker, markersize=7, color=INK, zorder=4)
    labelled_sites = {
        "top": ((0, 0), (-0.05, -0.08), "right"),
        "bridge": ((0.5, 0), (0.5, -0.08), "center"),
    }
    for label, ((u, w), (tx, ty), align) in labelled_sites.items():
        x, y = to_cartesian(u, w)
        ax.annotate(
            label,
            (x, y),
            xytext=(tx, ty),
            ha=align,
            va="top",
            fontsize=7,
            color=INK_2,
            arrowprops={"arrowstyle": "-", "color": INK_2, "linewidth": 0.6, "shrinkB": 9},
        )

    # The lattice vectors, from the origin, just outside the cell's edges.
    for (u, w), name, offset in (((1, 0), "a₁", (0, -0.17)), ((0, 1), "a₂", (-0.17, 0.0))):
        ex, ey = to_cartesian(u, w)
        ox, oy = offset
        ax.annotate(
            "",
            (ex + ox, ey + oy),
            xytext=(ox, oy),
            arrowprops={"arrowstyle": "-|>", "color": INK, "linewidth": 1.0, "shrinkA": 0, "shrinkB": 0},
        )
        ax.text(
            (ex / 2 + ox) - 0.04 * (u == 0),
            ey / 2 + oy - 0.04 * (u == 1),
            name,
            fontsize=9,
            color=INK,
            ha="right" if u == 0 else "center",
            va="top" if u == 1 else "center",
        )

    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=10,
            markerfacecolor=SURFACE,
            markeredgecolor=color,
            markeredgewidth=2.5,
            label=f"orbit of {size}: {SIZE_MEANING[size]}",
        )
        for size, color in SIZE_COLORS.items()
    ]
    handles += [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markersize=11,
            markerfacecolor=SURFACE,
            markeredgecolor=INK,
            markeredgewidth=2.0,
            label="representative (capital letter)",
        ),
        *(
            Line2D(
                [],
                [],
                marker=marker,
                linestyle="none",
                markersize=7,
                color=INK,
                label=f"{name} ({u.numerator}/{u.denominator}, {w.numerator}/{w.denominator})",
            )
            for name, (marker, (u, w)) in HOLLOWS.items()
        ),
        Line2D([], [], color=MUTED, linewidth=0.6, linestyle=(0, (4, 3)), label="mirror lines"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, borderaxespad=0.2, handletextpad=0.4)

    ax.set_title("C6v orbits on an 8 × 8 grid (nφ = 12)", color=INK)
    ax.set_aspect("equal")
    ax.set_xlim(-0.33, 1.53)
    ax.set_ylim(-0.27, 0.92)
    ax.axis("off")
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.01)
    fig.savefig(OUTPUT, dpi=DPI, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    check_against_readme()
    draw()
    print(f"wrote {OUTPUT}")

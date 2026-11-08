"""Shared figure style: white background, one fixed colour + line style + marker per program and
per code pair, legible at 880 px (GitHub content width).

Colours: the first three slots of the dataviz reference palette (blue, orange, aqua; documented
there to pass the all-pairs colour-vision checks), plus violet for the third code pair and green
for the strict-build GPU run (added 2026-11-08). Identity
never rests on colour alone: every series also has its own line style and marker, a legend, and
direct labels where they fit. (The palette validator needs node, which this machine lacks; the
slot choice follows the documented validation.)"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e3df"
SURFACE = "#ffffff"

PROGRAM = {
    "fortran": {"label": "Fortran 5.19 (gfortran -O3, FFTW)", "short": "Fortran", "color": "#2a78d6", "ls": "-", "marker": "o"},
    "cpu": {"label": "C++ CPU (build-fast: fast-math, native, oneMKL)", "short": "C++ CPU", "color": "#eb6834", "ls": "--", "marker": "s"},
    "gpu": {"label": "C++ GPU (build-cuda-fast: FMA kernels, cuFFT)", "short": "C++ GPU", "color": "#1baf7a", "ls": ":", "marker": "^"},
    "gpu_strict": {"label": "C++ GPU strict (build-cuda: no fast math, no FMA; added 2026-11-08)", "short": "C++ GPU strict", "color": "#008300", "ls": (0, (5, 1.5, 1, 1.5, 1, 1.5)), "marker": "h"},
    "thesis": {"label": "Thesis Fig. 5.3, DZ PES (dotted curve, digitised)", "short": "Thesis", "color": "#52514e", "ls": "-", "marker": "D"},
}
PAIR = {
    "cpu-fortran": {"label": "C++ CPU − Fortran", "color": "#eb6834", "ls": "--", "marker": "s"},
    "gpu-fortran": {"label": "C++ GPU − Fortran", "color": "#1baf7a", "ls": ":", "marker": "^"},
    "gpu-cpu": {"label": "C++ GPU − C++ CPU", "color": "#6250d6", "ls": "-.", "marker": "v"},
    "gpu_strict-fortran": {"label": "C++ GPU strict − Fortran", "color": "#008300", "ls": "-", "marker": "h"},
}


def setup():
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.size": 10.5, "axes.titlesize": 11.5, "axes.labelsize": 10.5, "legend.fontsize": 9,
        "xtick.labelsize": 9.5, "ytick.labelsize": 9.5,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK_2, "ytick.color": INK_2,
        "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.6,
        "legend.frameon": False, "savefig.dpi": 110, "font.family": "DejaVu Sans",
        # reproducible PNGs: no timestamps / software tags in the metadata
        "savefig.bbox": "tight",
    })


def save(fig, path):
    fig.savefig(path, metadata={"Software": None, "Creation Time": None}, pil_kwargs={"optimize": True})
    plt.close(fig)

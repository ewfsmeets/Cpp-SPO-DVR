#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib>=3.9", "PyQt6"]
# ///
"""Wave-packet designer, interactive window: the figure of `check --png`, with sliders for the recipe's knobs,
the verdict, and explanations.

Usage (from any directory; uv installs numpy, matplotlib and PyQt6 once, from the header above):
  uv run tools/wave_packet_designer/gui_mpl.py [input.json] [--zasymp Z --point-group G]
      [--threads N] [--representatives N [--operations N]] [--ram GiB] [--gpu GiB]
Without an input it opens inputs/thesis_fig5_3_dz/input.json. --ram and --gpu replace this machine's RAM and
free GPU memory (nvidia-smi) in the memory checks.

The window opens on the input as read. Moving a slider changes the input in memory; when dZ, ZSpOptStart or
an absorber length moves, the Z grids follow as the recipe places them, and tEnd moves in whole averaging
intervals. Hovering over a panel or a control shows its explanation (explanations.py); "Explain" lays the
explanations over the panels. "Suggest" applies the recipe of inputs/README.md (section 2), "Reset" returns
to the input, "Save" writes the changed input to the file named in the box (never over an existing file).
The window only wires: every number comes from evaluate() (README, "Architecture").
"""

import argparse
import math
import sys
import textwrap
from dataclasses import dataclass, fields, replace
from pathlib import Path

sys.dont_write_bytecode = True  # no __pycache__ in the repository, also when not run with -B
sys.path.insert(0, str(Path(__file__).resolve().parent))  # the tool's modules, found next to this file

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.backends import BackendFilter, backend_registry  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.widgets import Button, Slider, TextBox  # noqa: E402

from evaluation import SETTINGS_STAGE, evaluate, fitted_parameters, fitted_time, suggest  # noqa: E402
from explanations import CONTROLS, PANELS  # noqa: E402
from input_file import InputFile, InputFileError, load_input_file, save_changed_copy  # noqa: E402
from plots import INK, INK_2, add_panels, apply_style, draw_panels  # noqa: E402
from program.config import InputParameters  # noqa: E402
from program.json_input import key_of_field, parse_input  # noqa: E402
from program.point_group import GROUP_ORDER  # noqa: E402
from report import status_line, verdict_lines  # noqa: E402
from settings import Settings, this_machine  # noqa: E402

DEFAULT_INPUT = Path(__file__).resolve().parent.parent.parent / "inputs/thesis_fig5_3_dz/input.json"
UPDATE_DELAY = 250  # ms after the last slider move: one evaluation (about 0.5 s) per move, not per pixel
# Two sliders set lengths, not input keys: the absorbers then end on the last point of the fitted grids.
SCATTERING_LENGTH, SPECULAR_LENGTH = "scattering absorber length", "specular absorber length"
# The sliders whose change re-fits the Z grids, and the one that moves tEnd in whole averaging intervals.
GRID_SLIDERS = ("dz", "zsp_opt_start", SCATTERING_LENGTH, SPECULAR_LENGTH)
TIME_SLIDER = "t_end"
NOT_INPUT_KEYS = (SCATTERING_LENGTH, SPECULAR_LENGTH, TIME_SLIDER)  # sliders that are fitted, not copied
SLIDER_ROWS, SLIDER_COLUMNS = 6, 2
TEXT_WIDTH = 66  # characters per line in the verdict and help boxes
VERDICT_LINES = 20  # at most this many lines of the verdict are shown
HELP_PROMPT = "Hover over a panel, a slider or a button for its explanation."


@dataclass(frozen=True)
class SliderSpec:
    """One slider: an InputParameters field (or one of the two lengths) and its range."""

    label: str
    field: str
    low: float
    high: float
    step: float
    initial: float


def slider(label: str, field: str, low: float, high: float, step: float, initial: float) -> SliderSpec:
    """A slider whose range holds the input's value on its grid of steps: matplotlib's Slider moves the initial
    value onto low + n * step, so low is shifted to make the input's value such a point."""
    on_grid_low = initial - step * math.floor((initial - low) / step) if initial >= low else initial
    return SliderSpec(label, field, on_grid_low, max(high, initial), step, initial)


def slider_specs(p: InputParameters) -> list[SliderSpec]:
    """The sliders for an input: the recipe's knobs, around the input's values. A target energy has a slider
    only where the input tunes that absorber by energy."""
    specs = [
        slider("EMinWP (eV)", "e_min_wave_packet", 0.02, 3.0, 0.01, p.e_min_wave_packet),
        slider("EMaxWP (eV)", "e_max_wave_packet", 0.05, 6.0, 0.01, p.e_max_wave_packet),
        slider("dZ (bohr)", "dz", 0.02, 0.25, 0.005, p.dz),
        slider("Z0 (bohr)", "z0", p.z0 - 8.0, p.z0 + 8.0, 0.05, p.z0),
        slider(
            "ZSpOptStart (bohr)", "zsp_opt_start", p.zsp_opt_start - 8.0, p.zsp_opt_start + 8.0, 0.02, p.zsp_opt_start
        ),
        slider("ZOpt length (bohr)", SCATTERING_LENGTH, 1.0, 12.0, 0.02, p.z_opt_end - p.z_opt_start),
        slider("ZSpOpt length (bohr)", SPECULAR_LENGTH, 1.0, 20.0, 0.02, p.zsp_opt_end - p.zsp_opt_start),
        slider("ROptStart (bohr)", "r_opt_start", p.r_start, p.r_opt_end, 0.05, p.r_opt_start),
        slider("tEnd (a.u.)", TIME_SLIDER, p.t_start + p.dt_averages, 4.0 * p.t_end, p.dt_averages, p.t_end),
    ]
    energies = (
        ("ZOptEnergy (eV)", "z_opt_energy", 1.5),
        ("ZSpOptEnergy (eV)", "zsp_opt_energy", 3.0),
        ("ROptEnergy (eV)", "r_opt_energy", 1.5),
    )
    for label, field, high in energies:
        value = getattr(p, field)
        if value is not None:
            specs.append(slider(label, field, 0.02, high, 0.01, value))
    return specs


class DesignerWindow:
    """The panels, the sliders, the buttons, the verdict and the help box on one figure; every change
    re-evaluates and redraws."""

    def __init__(self, figure: Figure, input_file: InputFile, parameters: InputParameters, settings: Settings):
        self.figure, self.input_file, self.settings, self.original = figure, input_file, settings, parameters
        grid = figure.add_gridspec(2, 3, left=0.05, right=0.985, top=0.93, bottom=0.47, hspace=0.42, wspace=0.27)
        self.panels = add_panels(figure, grid)
        self.specs = slider_specs(parameters)
        if len(self.specs) > SLIDER_ROWS * SLIDER_COLUMNS:
            raise ValueError(f"{len(self.specs)} sliders do not fit the {SLIDER_ROWS} x {SLIDER_COLUMNS} places")
        self.sliders = {spec.field: self.add_slider(i, spec) for i, spec in enumerate(self.specs)}
        self.explaining = False
        self.add_buttons()
        self.message = figure.text(0.05, 0.05, "", color=INK_2, fontsize=8.5)
        self.verdict = figure.text(0.70, 0.415, "", color=INK, fontsize=7.8, va="top")
        self.help = figure.text(
            0.70,
            0.015,
            wrapped([self.help_for(None)], 14),
            color=INK,
            fontsize=7.8,
            va="bottom",
            bbox={"facecolor": "#f6f5f2", "edgecolor": "#e4e3df", "pad": 4},
        )
        # A slider move restarts this timer; when it runs out, update() evaluates once.
        self.timer = figure.canvas.new_timer(interval=UPDATE_DELAY)
        self.timer.single_shot = True
        self.timer.add_callback(self.update)
        figure.canvas.mpl_connect("motion_notify_event", self.on_hover)
        self.update()

    def add_slider(self, index: int, spec: SliderSpec) -> Slider:
        """A slider in one of two columns below the panels; its black tick marks the input's value."""
        column, row = divmod(index, SLIDER_ROWS)
        ax = self.figure.add_axes((0.13 + 0.33 * column, 0.385 - 0.042 * row, 0.19, 0.024))
        widget = Slider(
            ax, spec.label, spec.low, spec.high, valinit=spec.initial, valstep=spec.step, color=INK_2, initcolor=INK
        )
        widget.on_changed(self.schedule_update)
        return widget

    def add_buttons(self) -> None:
        """Suggest, Reset, Explain and Save, and the box with the file name for Save."""
        self.buttons: dict[str, Button] = {}  # kept: a widget without a reference stops responding
        actions = (
            ("Suggest", self.suggest),
            ("Reset", self.reset),
            ("Explain", self.toggle_explanations),
            ("Save", self.save),
        )
        for i, (name, action) in enumerate(actions):
            button = Button(self.figure.add_axes((0.05 + 0.08 * i, 0.105, 0.07, 0.035)), name)
            button.on_clicked(action)
            self.buttons[name] = button
        self.path_box = TextBox(self.figure.add_axes((0.40, 0.105, 0.26, 0.035)), "file ", initial="designed.json")

    def schedule_update(self, _value: object = None) -> None:
        """Restarts the update timer (on every slider move), and clears the last message."""
        self.message.set_text("")
        self.timer.stop()
        self.timer.start()

    def slider_values(self) -> dict[str, float]:
        """Each slider's value, rounded to 10 decimals: a slider adds its steps in binary (0.39999999999999997)."""
        return {spec.field: round(float(self.sliders[spec.field].val), 10) for spec in self.specs}

    def current_parameters(self) -> InputParameters:
        """The input with the sliders' values. The Z grids are re-fitted only when a grid slider differs from
        the input, and tEnd only when its slider does: an untouched window is the input as read."""
        values = self.slider_values()
        initial = {spec.field: round(spec.initial, 10) for spec in self.specs}
        p = replace(self.original, **{f: v for f, v in values.items() if f not in NOT_INPUT_KEYS})
        if any(values[field] != initial[field] for field in GRID_SLIDERS):
            p = fitted_parameters(p, self.settings, values[SCATTERING_LENGTH], values[SPECULAR_LENGTH])
        if values[TIME_SLIDER] != initial[TIME_SLIDER]:
            p = fitted_time(p, values[TIME_SLIDER])
        return p

    def update(self, _value: object = None) -> None:
        """Evaluates the current input and redraws the panels, the title, the verdict and the explanations."""
        self.parameters = self.current_parameters()
        self.evaluation = evaluate(self.parameters, self.settings)
        draw_panels(self.panels, self.evaluation)
        if self.explaining:
            self.draw_explanations()
        self.figure.suptitle(status_line(self.evaluation), color=INK, fontsize=10.5)
        self.verdict.set_text(wrapped(verdict_lines(self.evaluation), VERDICT_LINES))
        self.figure.canvas.draw_idle()

    def set_sliders(self, p: InputParameters) -> None:
        """Moves the sliders to the values of `p` (widening a range that does not hold one), then updates."""
        for spec in slider_specs(p):
            if spec.field in self.sliders:
                widget = self.sliders[spec.field]
                if not widget.valmin <= spec.initial <= widget.valmax:
                    widget.valmin, widget.valmax = min(widget.valmin, spec.initial), max(widget.valmax, spec.initial)
                    widget.ax.set_xlim(widget.valmin, widget.valmax)
                widget.eventson = False
                widget.set_val(spec.initial)
                widget.eventson = True
        self.update()

    def suggest(self, _event: object = None) -> None:
        """Applies the recipe of inputs/README.md to the current input."""
        s = suggest(self.parameters, self.settings)
        if s.blockers:
            self.message.set_text("Suggest needs first: " + "; ".join(f"{f.key}: {f.message}" for f in s.blockers))
            self.figure.canvas.draw_idle()
            return
        self.set_sliders(replace(self.parameters, **s.changes))
        self.message.set_text("Suggested by the recipe of inputs/README.md, section 2.")

    def reset(self, _event: object = None) -> None:
        """Back to the input as read."""
        self.set_sliders(self.original)
        self.message.set_text("Back to the input as read.")

    def save(self, _event: object = None) -> None:
        """Writes the input with every changed key to the file named in the box (never over an existing one)."""
        self.timer.stop()
        self.update()  # a slider moved less than UPDATE_DELAY ago: save what it shows
        names = [f.name for f in fields(InputParameters)]
        changed = [name for name in names if getattr(self.parameters, name) != getattr(self.original, name)]
        changes = {key_of_field(name): getattr(self.parameters, name) for name in changed}
        try:
            target = save_changed_copy(self.input_file, changes, Path(self.path_box.text))
            self.message.set_text(f"Written to {target} ({len(changes)} keys changed).")
        except InputFileError as error:
            self.message.set_text(f"Not written: {error}")
        self.figure.canvas.draw_idle()

    def toggle_explanations(self, _event: object = None) -> None:
        """Shows or hides the explanations on the panels."""
        self.explaining = not self.explaining
        self.update()

    def draw_explanations(self) -> None:
        """Each panel's explanation over the panel, in a box (the panels are redrawn without it)."""
        for ax, draw in self.panels:
            width = 150 if ax.get_position().width > 0.5 else 66
            ax.text(
                0.01,
                0.98,
                textwrap.fill(PANELS[draw.__name__], width),
                transform=ax.transAxes,
                va="top",
                fontsize=7.5,
                color=INK,
                zorder=10,
                bbox={"facecolor": "white", "alpha": 0.93, "edgecolor": "#e4e3df", "pad": 4},
            )

    def help_for(self, ax: Axes | None) -> str:
        """The explanation of the panel, slider or button that `ax` belongs to; elsewhere the prompt and what
        the verdict is."""
        for panel_ax, draw in self.panels:
            if ax is panel_ax:
                return PANELS[draw.__name__]
        for field, widget in self.sliders.items():
            if ax is widget.ax:
                return CONTROLS[field] + "\n" + CONTROLS["sliders"]
        for name, button in self.buttons.items():
            if ax is button.ax:
                return CONTROLS[name]
        return HELP_PROMPT + "\n" + CONTROLS["verdict"]

    def on_hover(self, event: object) -> None:
        """Shows the explanation of what the mouse is over in the help box."""
        text = wrapped([self.help_for(getattr(event, "inaxes", None))], 14)
        if text != self.help.get_text():
            self.help.set_text(text)
            self.figure.canvas.draw_idle()


def wrapped(paragraphs: list[str], max_lines: int) -> str:
    """The paragraphs wrapped to TEXT_WIDTH, at most `max_lines` lines (then an ellipsis)."""
    lines = []
    for paragraph in paragraphs:
        for part in paragraph.split("\n"):
            lines += textwrap.wrap(part, TEXT_WIDTH, subsequent_indent="  ") or [""]
    return "\n".join(lines if len(lines) <= max_lines else lines[: max_lines - 1] + ["  ..."])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", nargs="?", default=str(DEFAULT_INPUT), help="input.json (default: the DZ input)")
    parser.add_argument("--zasymp", type=float, help="Z_asymptotic (bohr) of a CRP data set")
    parser.add_argument("--point-group", choices=sorted(GROUP_ORDER), help="point group of a CRP data set")
    parser.add_argument("--threads", type=positive, default=Settings.threads, help="OpenMP threads (default 8)")
    parser.add_argument("--representatives", type=positive, help='with UseSymmetry: "(x, y) points computed"')
    parser.add_argument("--operations", type=positive, help='the names under "Operations used", identity included')
    parser.add_argument("--ram", type=float, help="memory for the run (GiB); default: this machine's RAM")
    parser.add_argument("--gpu", type=float, help="GPU memory for the run (GiB); default: this GPU's free memory")
    args = parser.parse_args(argv)
    if plt.get_backend().lower() in backend_registry.list_builtin(BackendFilter.NON_INTERACTIVE):
        print(
            f"error: matplotlib has no window backend here ({plt.get_backend()}); run this file with "
            "`uv run tools/wave_packet_designer/gui_mpl.py`, which installs PyQt6, on a machine with a display",
            file=sys.stderr,
        )
        return 2
    try:
        input_file = load_input_file(Path(args.input))
    except InputFileError as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    parsed = parse_input(input_file.values)
    if parsed.parameters is None:
        print(f"error: spodvr cannot read {args.input}: run `check` on it", file=sys.stderr)
        return 1
    ram, gpu = this_machine()
    ram = int(args.ram * 2**30) if args.ram is not None else ram
    gpu = int(args.gpu * 2**30) if args.gpu is not None else gpu
    settings = Settings(args.zasymp, args.point_group, args.threads, args.representatives, args.operations, ram, gpu)
    first = evaluate(parsed.parameters, settings)
    if first.stopped_at == SETTINGS_STAGE:  # a CRP data set without its values, or impossible symmetry counts
        for finding in first.findings:
            print(f"error: {finding.message}", file=sys.stderr)
        return 2
    apply_style()
    figure = plt.figure(figsize=(16.0, 11.0))
    figure.canvas.manager.set_window_title(f"Wave-packet designer: {args.input}")
    window = DesignerWindow(figure, input_file, parsed.parameters, settings)  # noqa: F841 (keeps the widgets alive)
    plt.show()
    return 0


def positive(text: str) -> int:
    """An integer >= 1 (argparse type)."""
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError(f"{text} is not a positive integer")
    return value


if __name__ == "__main__":
    sys.exit(main())

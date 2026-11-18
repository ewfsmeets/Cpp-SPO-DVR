"""Case plots: the figure (plots.py), the verdict text (report.py), `check --png`, and the window's wiring
(gui_mpl.py), all without a window: test_designer.py selects matplotlib's Agg backend."""

import os
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt

from evaluation import evaluate, suggest
from explanations import CONTROLS, PANELS
from gui_mpl import DesignerWindow
from input_file import load_input_file
from plots import add_panels, apply_style, report_figure
from program.json_input import parse_input
from report import status_line, verdict_lines
from settings import Settings
from tests.helpers import DZ_INPUT, ER_INPUT, HARDWALL_INPUT, TOOL_DIR, check, check_equal, load_parameters


def texts(ax: plt.Axes) -> str:
    """All text drawn on the axes, joined."""
    return " ".join(t.get_text() for t in ax.texts)


def test_panels_draw_their_content(work: Path) -> None:
    """For the DZ input (with 32 GiB RAM and 9 GiB free GPU memory given): the Z panel has both grid lanes, the
    packet and its marks; the momentum panel the packet and +-k_grid; the time panel Zana, Z_turn, tEnd and the
    elastic and E_n,min recrossings; the absorber panel R and T of its three absorbers; the r panel one curve
    per analysed v, the limit and the two legend entries; the memory panel a bar and a capacity mark for the
    CPU and the GPU; the verdict panel the verdict."""
    e = evaluate(load_parameters(DZ_INPUT), Settings(representatives=19, ram_bytes=32 * 2**30, gpu_bytes=9 * 2**30))
    figure = report_figure(e)
    z, k, time, absorbers, r, memory, verdict = figure.axes
    for expected in ("scattering grid", "specular grid", "Zana", "Z0", "outer tail", "allowed region"):
        check(expected in texts(z), f"Z panel: {expected} missing")
    check(len(z.lines) >= 1 and len(k.lines) >= 1, "packet curves")
    check("k_grid" in texts(k) and "outward" in texts(k), texts(k))
    labels = sorted(line.get_label() for line in absorbers.lines if not line.get_label().startswith("_"))
    check_equal(labels, sorted(f"{a} {x}" for a in ("ROpt", "ZOpt", "ZSpOpt") for x in "RT"), "R and T curves")
    check_equal(sum(1 for t in r.texts if t.get_text().startswith(" v=")), len(e.bound_state_tails), "v labels")
    check("limit 1e-6" in texts(r) and "ROptStart" in texts(r), texts(r))
    legend = [t.get_text() for t in r.get_legend().get_texts()]
    check(len(legend) == 2 and "converged" in legend[1], f"r legend {legend}")
    check(len(memory.patches) == 2 and "RAM" in texts(memory) and "free GPU" in texts(memory), texts(memory))
    for expected in ("Zana", "Z_turn", "tEnd", "elastic", "E_n 0.03 eV"):
        check(expected in texts(time), f"time panel: {expected} missing")
    check("spodvr accepts the input." in texts(verdict), "verdict panel")
    plt.close(figure)


def test_figure_of_each_kind_of_input(work: Path) -> None:
    """The figure draws without error for the repository inputs, a potential-only run (no packet) and an input
    the program rejects (Z0 16.78, the TolSp check)."""
    for name, p in {
        "er": load_parameters(ER_INPUT),
        "hardwall": load_parameters(HARDWALL_INPUT),
        "potential_only": load_parameters(DZ_INPUT, only_potential=True),
        "tolsp": load_parameters(DZ_INPUT, z0=16.78),
    }.items():
        figure = report_figure(evaluate(p, Settings()))
        figure.savefig(work / f"{name}.png", dpi=50)
        plt.close(figure)


def test_title_and_verdict_give_the_reason(work: Path) -> None:
    """When spodvr stops, the title gives its first error (ZSpOptStart 14.08 lies below Z0: "Z0 > ZSpOptStart.")
    and the verdict says the design is not checked; an accepted input lists what it misses with value and limit."""
    stopped = evaluate(load_parameters(DZ_INPUT, zsp_opt_start=14.08), Settings())
    check("spodvr stops: Z0: Z0 > ZSpOptStart." in status_line(stopped), status_line(stopped))
    check(any("not checked" in line for line in verdict_lines(stopped)), "design not checked")
    lines = verdict_lines(evaluate(load_parameters(DZ_INPUT), Settings()))
    check(any("lowest outgoing normal energy" in line and "needs <= 0.03" in line for line in lines), str(lines))


def test_check_draws_the_figure(work: Path) -> None:
    """`check --png` writes the figure and keeps the report and the exit status of `check`."""
    target = work / "check.png"
    environment = os.environ | {"MPLCONFIGDIR": str(work / "matplotlib")}
    command = [sys.executable, "-I", "-B", str(TOOL_DIR / "wave_packet_designer.py"), "check", str(DZ_INPUT)]
    result = subprocess.run([*command, "--png", str(target)], capture_output=True, text=True, env=environment)
    check(result.returncode == 0 and target.exists(), f"{result.returncode}: {result.stdout}{result.stderr}")
    check("spodvr accepts the input" in result.stdout, "the report")


def new_window(path: Path = DZ_INPUT) -> DesignerWindow:
    """The window on an input, drawn on a figure of the Agg backend."""
    input_file = load_input_file(path)
    apply_style()
    parameters = parse_input(input_file.values).parameters
    return DesignerWindow(plt.figure(figsize=(16.0, 11.0)), input_file, parameters, Settings())


def test_window_opens_and_resets_to_the_input(work: Path) -> None:
    """The window shows the input as read, also one not laid out by the recipe (hardwall_small), and Reset
    returns to it after a change."""
    for path in (DZ_INPUT, HARDWALL_INPUT):
        window = new_window(path)
        check_equal(window.parameters, window.original, f"{path.parent.name} on opening")
        window.sliders["dz"].set_val(0.1)
        window.update()  # the timer that does this after a slider move does not run on the Agg backend
        check(window.parameters != window.original, "changed by the slider")
        window.reset()
        check_equal(window.parameters, window.original, f"{path.parent.name} after Reset")


def test_window_follows_its_sliders(work: Path) -> None:
    """Moving dZ to 0.1 re-fits the Z grids: Zana stays on Z_asymptotic = 15 (point 150), the scattering
    absorber of 4.08 bohr then ends on the last point of nZ 192 (4.1 bohr on that grid); moving tEnd to 25 100
    gives 25 500 (whole averaging intervals of 500) with nStep 10 200; the program accepts the input. dR 0.1
    with nR 80 moves ROptEnd to the last r point."""
    window = new_window()
    window.sliders["dz"].set_val(0.1)
    window.sliders["t_end"].set_val(25100.0)
    window.update()
    p = window.parameters
    check_equal((p.dz, p.z_analysis, p.nz, round(p.z_opt_end - p.z_opt_start, 10)), (0.1, 15.0, 192, 4.1), "Z")
    check_equal((p.t_end, p.n_steps), (25500.0, 10200), "tEnd and nStep")
    check(window.evaluation.accepted, f"{window.evaluation.findings}")
    window.sliders["dr"].set_val(0.1)
    window.sliders["nr"].set_val(80)
    window.update()
    check_equal((window.parameters.nr, window.parameters.r_opt_end), (80, 8.4), "nR 80, dR 0.1: ROptEnd 0.5 + 79 dR")


def test_suggest_is_not_clamped_by_the_sliders(work: Path) -> None:
    """Suggest applies the recipe's values even beyond a slider's range (EMaxWP 0.45: Z0 about 40 bohr, far
    beyond the input's Z0 + 8): the window's input equals the suggestion."""
    window = new_window()
    window.sliders["e_max_wave_packet"].set_val(0.45)
    window.update()
    expected = replace(window.parameters, **suggest(window.parameters, Settings()).changes)
    window.suggest()
    check_equal(window.parameters, expected, "the window after Suggest")
    check(window.evaluation.accepted, f"{window.evaluation.findings}")


def test_window_saves_the_changed_keys(work: Path) -> None:
    """Save writes the changed keys (the copy reads back as the window's input) and refuses an existing file."""
    window = new_window()
    window.sliders["e_max_wave_packet"].set_val(1.2)
    window.suggest()
    target = work / "designed.json"
    window.path_box.set_val(str(target))
    window.save()
    check_equal(parse_input(load_input_file(target).values).parameters, window.parameters, "the saved copy")
    window.save()
    check("Not written" in window.message.get_text(), window.message.get_text())


def test_e_n_min_slider_sets_the_time_check(work: Path) -> None:
    """The DZ input's tEnd 17 500 converges outgoing energies down to 0.06 eV: with E_n,min 0.03 eV (the default)
    the time check fails, with 0.07 eV it passes."""
    window = new_window()
    rule = "lowest outgoing normal energy converged by tEnd"
    check(rule in [c.name for c in window.evaluation.design if c.passed is False], "fails at 0.03 eV")
    window.sliders["lowest_outgoing_energy"].set_val(0.07)
    window.update()
    check(rule not in [c.name for c in window.evaluation.design if c.passed is False], "passes at 0.07 eV")


def test_everything_has_an_explanation(work: Path) -> None:
    """Every panel, slider and button has an explanation, hovering over one finds it, and Explain lays the
    panel texts over the panels."""
    figure = plt.figure()
    drawn = [draw.__name__ for _, draw in add_panels(figure, figure.add_gridspec(2, 4))]
    check_equal(sorted(drawn), sorted(PANELS), "panels with an explanation")
    window = new_window()
    for field, widget in window.sliders.items():
        check(field in CONTROLS and window.help_for(widget.ax).startswith(CONTROLS[field]), f"slider {field}")
    for name, button in window.buttons.items():
        check_equal(window.help_for(button.ax), CONTROLS[name], f"button {name}")
    check_equal(window.help_for(window.panels[0][0]), PANELS["draw_z_layout"], "the Z panel")
    window.toggle_explanations()
    check(any(t.get_text().startswith("Z: the two") for t in window.panels[0][0].texts), "Explain shows the text")


TESTS = [
    test_panels_draw_their_content,
    test_figure_of_each_kind_of_input,
    test_title_and_verdict_give_the_reason,
    test_check_draws_the_figure,
    test_window_opens_and_resets_to_the_input,
    test_window_follows_its_sliders,
    test_suggest_is_not_clamped_by_the_sliders,
    test_window_saves_the_changed_keys,
    test_e_n_min_slider_sets_the_time_check,
    test_everything_has_an_explanation,
]

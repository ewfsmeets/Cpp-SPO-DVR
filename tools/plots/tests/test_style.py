"""The look: the scale rule, colours only in style.py, and the designer's colour constants."""

import ast
import re
from pathlib import Path

import numpy as np

import style

TOOL_DIR = Path(__file__).resolve().parent.parent
DESIGNER_PLOTS = TOOL_DIR.parent / "wave_packet_designer" / "plots.py"


def test_wide_range_gets_a_log_scale_with_a_floor() -> None:
    """Values spanning more than 1e3 get a log axis whose lower limit is max(smallest positive, 1e-12 max);
    the values <= 0 a log axis cannot show are counted."""
    scale = style.scale_and_limits([np.array([1e-20, 1e-5, 1.0, 0.0, -1.0])])
    assert scale.kind == "log" and scale.lower == 1e-12 and scale.not_shown == 2


def test_narrow_and_signed_ranges_stay_linear() -> None:
    """A range of exactly LOG_RANGE, and a signed quantity, are drawn on a linear axis; just above, log."""
    assert style.scale_and_limits([np.array([1.0, style.LOG_RANGE])]).kind == "linear"
    assert style.scale_and_limits([np.array([1.0, 1.001 * style.LOG_RANGE])]).kind == "log"
    assert style.scale_and_limits([np.array([-1e-6, 1.0])], signed=True).kind == "linear"


def test_masked_hides_values_below_the_floor() -> None:
    """Values below the floor become NaN (not drawn)."""
    assert np.isnan(style.masked(np.array([1e-9, 1.0]), 1e-8)[0])


def test_no_colour_outside_style() -> None:
    """No hex colour or colour-map name in any module of the tool but style.py (and this test)."""
    pattern = re.compile(r"#[0-9a-fA-F]{6}\b|\b(viridis|cividis|PuOr|jet|plasma|magma|inferno|Blues|Reds|tab10)\b")
    found = []
    for path in sorted(TOOL_DIR.rglob("*.py")):
        if path.name in ("style.py", "test_style.py") or "__pycache__" in path.parts:
            continue
        found += [f"{path.relative_to(TOOL_DIR)}: {m.group(0)}" for m in pattern.finditer(path.read_text())]
    assert not found, "\n".join(found)


def designer_constants() -> dict[str, object]:
    """The constants the designer's plots.py assigns at module level as `A, B = "#...", "#..."` or `A = ...`,
    read with ast and literal_eval (the two tools never share a process)."""
    constants = {}
    for node in ast.parse(DESIGNER_PLOTS.read_text()).body:
        if not (isinstance(node, ast.Assign) and len(node.targets) == 1):
            continue
        target = node.targets[0]
        names = [t.id for t in target.elts] if isinstance(target, ast.Tuple) else [getattr(target, "id", "")]
        try:
            values = ast.literal_eval(node.value)
        except ValueError:  # not a literal (a dict of names, a call): not a colour constant
            continue
        constants.update(zip(names, values if len(names) > 1 else [values], strict=True))
    return constants


def test_colours_equal_the_designers() -> None:
    """INK … AQUA have the designer's values (the same look and meanings in both tools)."""
    theirs = designer_constants()
    for name in ("INK", "INK_2", "MUTED", "GRID", "SURFACE", "BLUE", "ORANGE", "AQUA"):
        assert getattr(style, name) == theirs[name], name

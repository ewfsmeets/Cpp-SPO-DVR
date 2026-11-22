"""Tests (ctest WavePacketDesigner.InputFile): reading, and saving changed values without touching the rest of the text
(input_file.py)."""

import json
from dataclasses import replace
from pathlib import Path

from input_file import InputFileError, load_input_file, parse_json_object, save_changed_copy
from program.json_input import key_of_field, parse_input
from tests.helpers import DZ_INPUT, check, check_equal


def test_saved_copy_changes_only_the_values(work: Path) -> None:
    """Changing dZ and nZ rewrites only the lines that hold them; nZ stays a JSON integer."""
    source = load_input_file(DZ_INPUT)
    target = save_changed_copy(source, {"dZ": 0.1, "nZ": 192}, work / "changed.json")
    old_lines, new_lines = source.text.splitlines(), target.read_text().splitlines()
    check_equal(len(new_lines), len(old_lines), "number of lines")
    changed = [i for i, (a, b) in enumerate(zip(old_lines, new_lines, strict=True)) if a != b]
    check_equal(len(changed), 1, "changed lines (dZ and nZ share one line)")
    values = json.loads(target.read_text())
    check(values["dZ"] == 0.1 and values["nZ"] == 192 and isinstance(values["nZ"], int), "new values")


def test_new_key_is_added(work: Path) -> None:
    """A key that is not in the file goes before the closing brace."""
    target = save_changed_copy(load_input_file(DZ_INPUT), {"ZAnaSp": 21.0}, work / "added.json")
    check_equal(json.loads(target.read_text())["ZAnaSp"], 21.0, "added key")


def test_changed_parameters_are_saved_under_their_keys(work: Path) -> None:
    """The front ends' way of saving: changed InputParameters fields go into a copy under their JSON keys
    (key_of_field), and reading the copy gives the changed parameters."""
    source = load_input_file(DZ_INPUT)
    changed = replace(parse_input(source.values).parameters, z0=18.5, e_max_wave_packet=1.05, nz_specular=270)
    changes = {key_of_field(name): getattr(changed, name) for name in ("z0", "e_max_wave_packet", "nz_specular")}
    check_equal(sorted(changes), ["EMaxWP", "Z0", "nZSpecular"], "keys")
    target = save_changed_copy(source, changes, work / "front_end.json")
    check_equal(parse_input(load_input_file(target).values).parameters, changed, "parameters read back")


def test_existing_file_is_not_overwritten(work: Path) -> None:
    """Saving refuses to replace an existing file unless allowed."""
    source = load_input_file(DZ_INPUT)
    save_changed_copy(source, {"Z0": 18.5}, work / "copy.json")
    try:
        save_changed_copy(source, {"Z0": 19.0}, work / "copy.json")
    except InputFileError:
        save_changed_copy(source, {"Z0": 19.0}, work / "copy.json", overwrite=True)
        check_equal(json.loads((work / "copy.json").read_text())["Z0"], 19.0, "overwritten when allowed")
        return
    check(False, "an existing file was overwritten")


def test_comments_and_nan_are_refused() -> None:
    """Comments (which spodvr accepts) and NaN (which it refuses) give a clear error, not a wrong value."""
    for text, expected in (('{"a": 1 // note\n}', "comments"), ('{"a": NaN}', "NaN")):
        try:
            parse_json_object(text, "test")
        except InputFileError as error:
            check(expected in str(error), f"message for {expected}: {error}")
            continue
        check(False, f"{expected} was accepted")


def test_repeated_key_is_refused() -> None:
    """A key given twice in one object is an error, also inside a list entry (spodvr and the Fortran refuse it too);
    the same key in two different objects is not."""
    repeated = (
        '{"nZ": 176, "dZ": 0.08, "nZ": 192}',
        '{"a": 1, "b": {}, "a": 2}',
        '{"PlotSites": [{"name": "a", "x": 0, "x": 1}]}',
        '{"a": [[{"b": 1, "b": 2}]]}',
    )
    for text in repeated:
        try:
            parse_json_object(text, "test.json")
        except InputFileError as error:
            check("given twice" in str(error) and "test.json" in str(error), f"message for {text}: {error}")
            continue
        check(False, f"a repeated key was accepted: {text}")
    check_equal(parse_json_object('{"a": {"c": 1}, "c": 2}', "t"), {"a": {"c": 1}, "c": 2}, "same key, two objects")
    check_equal(parse_json_object('{"a": [{"b": 1}, {"b": 2}]}', "t"), {"a": [{"b": 1}, {"b": 2}]}, "two entries")

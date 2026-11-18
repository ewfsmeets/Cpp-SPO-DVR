"""Reads an input.json and saves changed values back into a copy of its text.

Part of the data layer (README, "Architecture"). Saving edits the text in place, so the layout of the file
(several keys per line, the order, the blank lines) and the integer keys stay as they were; only the changed
values differ.
"""

import json
import re
from dataclasses import dataclass
from pathlib import Path


class InputFileError(Exception):
    """An input file the designer cannot read or write."""


@dataclass(frozen=True)
class InputFile:
    """The text of an input file and its JSON object."""

    path: Path
    text: str
    values: dict


def load_input_file(path: Path) -> InputFile:
    """Reads and parses an input file; refuses comments and NaN/Infinity, which json.loads would treat differently."""
    try:
        text = path.read_text()
    except OSError as error:
        raise InputFileError(f"cannot read {path}: {error.strerror}") from error
    return InputFile(path, text, parse_json_object(text, str(path)))


def parse_json_object(text: str, name: str) -> dict:
    """The JSON object of an input file's text."""

    def refuse_constant(constant: str) -> None:
        raise InputFileError(f"{name}: {constant} is not valid JSON (spodvr refuses it too)")

    try:
        values = json.loads(text, parse_constant=refuse_constant)
    except json.JSONDecodeError as error:
        hint = (
            " (comments, which spodvr accepts, are not supported by the designer)"
            if "//" in text or "/*" in text
            else ""
        )
        raise InputFileError(f"{name}: invalid JSON at line {error.lineno}: {error.msg}{hint}") from error
    if not isinstance(values, dict):
        raise InputFileError(f"{name}: the input must be a JSON object")
    return values


# One JSON value: a string, a number, true, false or null.
JSON_VALUE = r'"(?:[^"\\]|\\.)*"|-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?|true|false|null'


def with_changed_values(text: str, changes: dict[str, object]) -> str:
    """The text with each key's value replaced, or the key added before the closing brace if it is absent."""
    for key, value in changes.items():
        encoded = json.dumps(value)
        pattern = re.compile(r'("' + re.escape(key) + r'"\s*:\s*)(' + JSON_VALUE + r")")
        matches = list(pattern.finditer(text))
        if len(matches) > 1:
            raise InputFileError(f"key {key} occurs {len(matches)} times")
        if matches:
            value_start, value_end = matches[0].span(2)  # group 2 is the old value
            text = text[:value_start] + encoded + text[value_end:]
        else:
            text = add_key(text, key, encoded)
    return text


def add_key(text: str, key: str, encoded_value: str) -> str:
    """Adds `"key": value` as a new line before the object's closing brace."""
    closing = text.rstrip().rfind("}")
    body = text[:closing].rstrip()
    separator = "," if not body.endswith("{") else ""
    return f'{body}{separator}\n  "{key}": {encoded_value}\n' + text[closing:]


def save_changed_copy(source: InputFile, changes: dict[str, object], target: Path, overwrite: bool = False) -> Path:
    """Writes the source text with `changes` applied to `target`; refuses to replace an existing file unless asked.

    `changes` maps JSON keys to new values (program/json_input.py key_of_field gives the key of a parameter
    field). The front ends' "Save as" goes through here.
    """
    if target.exists() and not overwrite:
        raise InputFileError(f"{target} exists; choose another name or allow overwriting")
    text = with_changed_values(source.text, changes)
    parse_json_object(text, str(target))  # never write a file that is not valid JSON
    target.write_text(text)
    return target

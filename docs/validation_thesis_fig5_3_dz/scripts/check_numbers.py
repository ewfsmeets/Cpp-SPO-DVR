#!/usr/bin/env python3
"""Number check: every number in the prose of README.md, stages.md and appendix.md (outside the
generated tables, code spans and code blocks) must match a value in scripts/results.json at the
precision written, or be a listed constant (criteria, item numbers, input labels).
Prints the numbers that match nothing; exit 1 if there are any."""
from __future__ import annotations

import json
import math
import re
import sys

from common import DOC, RESULTS

FILES = ["README.md", "stages.md", "appendix.md"]
# constants that are definitions, not results: criteria and tolerances, digits, item numbers,
# physical labels (v=0, j=0, Fig. 5.3), units and formats, citation volumes and pages
CONSTANTS = {1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-6, 1e-8, 1e-4, 1e-19, 1e-2, 0.01, 1.01, 2.0, 5.3, 1e-15, 1e-16,
             64.0, 880.0, 2026.0, 1.0, 0.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
             17.0, 18.0, 19.0, 20.0, 0.0, 100.0, 1995.0, 6280.0, 102.0, 360.0, 390.0}
SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def numbers_in_results():
    vals = []

    def walk(x):
        if isinstance(x, bool):
            return
        if isinstance(x, (int, float)):
            vals.append(float(x))
        elif isinstance(x, dict):
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
        elif isinstance(x, str):
            for m in re.findall(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", x):
                try:
                    vals.append(float(m))
                except ValueError:
                    pass
    walk(json.loads(RESULTS.read_text()))
    return vals


def prose(text):
    text = re.sub(r"<!-- BEGIN table:(\S+) -->.*?<!-- END table:\1 -->", "", text, flags=re.S)
    text = re.sub(r"```.*?```", "", text, flags=re.S)
    text = re.sub(r"`[^`]*`", "", text)
    text = re.sub(r"\]\([^)]*\)", "]", text)          # link targets
    text = re.sub(r'<a id="[^"]*"></a>', "", text)
    text = re.sub(r"\(\d{4}\)", "", text)
    text = re.sub(r"^\|.*\|$", "", text, flags=re.M)  # hand-written tables are checked too, see below
    return text


TOKEN = re.compile(r"(?<![\w.])([-+−]?\d+(?:[.,]\d+)?(?:[eE][-+−]?\d+)?)(?:\s*×\s*10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+))?(\s*%)?(\s+million)?")


def half_unit(written: str) -> float:
    """Half a unit of the last written digit: 1.4e-11 -> 5e-13, 0.0072 -> 5e-5, 636373 -> 0.5."""
    t = written.lower().replace("−", "-").lstrip("+-")
    mant, _, exp = t.partition("e")
    decimals = len(mant.split(".")[1]) if "." in mant else 0
    return 0.5 * 10.0 ** ((int(exp) if exp else 0) - decimals)


def matches(x, written, vals, scale=1.0):
    """x matches if some value (or its percentage / per-hour / ratio form) rounds to the written
    text: |v - x| <= 0.51 x half a unit of the last written digit."""
    if x in CONSTANTS or abs(x) in CONSTANTS:
        return True
    tol = 1.02 * half_unit(written) * scale
    for v in vals:
        for cand in (v, abs(v), 100 * abs(v), v / 3600, 1000 * v):
            if abs(cand - x) <= tol or abs(cand - abs(x)) <= tol:
                return True
    return False


def main():
    vals = numbers_in_results()
    unmatched = []
    for name in FILES:
        text = (DOC / name).read_text()
        body = prose(text)
        # also the hand-written tables (not generated): rows outside table markers
        tables = "\n".join(l for l in re.sub(r"<!-- BEGIN table:(\S+) -->.*?<!-- END table:\1 -->", "", text, flags=re.S)
                           .splitlines() if l.startswith("|"))
        for part in (body, re.sub(r"`[^`]*`", "", tables)):
            for line in part.splitlines():
                for m in TOKEN.finditer(line):
                    raw = m.group(1).replace("−", "-").replace(",", "")
                    try:
                        x = float(raw)
                    except ValueError:
                        continue
                    scale = 1.0
                    if m.group(2):
                        scale = 10 ** int(m.group(2).translate(SUP))
                        x = x * scale
                    if m.group(3):
                        x = x / 100
                        scale = 0.01
                    if m.group(4):
                        x = x * 1e6
                        scale = 1e6
                    if not matches(x, raw, vals, scale):
                        unmatched.append((name, m.group(0).strip(), line.strip()[:110]))
    seen = set()
    for name, tok, ctx in unmatched:
        if (name, tok) in seen:
            continue
        seen.add((name, tok))
        print(f"{name}: {tok!r:14s} … {ctx}")
    print(f"check_numbers: {len(seen)} number(s) without a match in results.json")
    sys.exit(1 if seen else 0)


if __name__ == "__main__":
    main()

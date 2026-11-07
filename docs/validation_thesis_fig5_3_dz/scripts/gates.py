#!/usr/bin/env python3
"""Phase 0: summary of the regression-gate runs (ctest -R Gate) of the two builds that ran,
plus the strict build's run of the one gate the fast builds fail. Reads the saved ctest output
in runs/thesis_fig5_3_dz_validation/gates/; writes data/items/gates.json and a copy of the
summary lines to data/gates_ctest_summary.txt."""
from __future__ import annotations

import re

from common import DATA, WORK, save_item

FILES = {"build-fast": "gates_build-fast.txt", "build-cuda-fast": "gates_build-cuda-fast.txt",
         "build (strict), hard-wall vasy gate only": "gates_build_hardwall_vasy.txt"}


def parse(text: str) -> dict:
    tests = []
    for m in re.finditer(r"Test +#(\d+): (\S+) \.+\s*(\*{0,3}\w+)\s+([\d.]+) sec", text):
        tests.append({"id": int(m.group(1)), "name": m.group(2), "result": m.group(3).strip("*"),
                      "seconds": float(m.group(4))})
    failures = []
    for m in re.finditer(r"actual: ([\d.e+-]+) vs ([\d.e+-]+)\n(.*)", text):
        failures.append({"actual": float(m.group(1)), "limit": float(m.group(2)), "detail": m.group(3).strip()})
    return {"tests": tests, "failures": failures,
            "passed": sum(t["result"] == "Passed" for t in tests),
            "failed": sum(t["result"] == "Failed" for t in tests),
            "skipped": sum(t["result"] == "Skipped" for t in tests)}


def main():
    out = {}
    lines = []
    for label, name in FILES.items():
        text = (WORK / "gates" / name).read_text()
        r = parse(text)
        out[label] = r
        lines.append(f"## {label}: {r['passed']} passed, {r['failed']} failed, {r['skipped']} skipped")
        for t in r["tests"]:
            lines.append(f"{t['result']:8s} {t['name']}")
        for f in r["failures"]:
            lines.append(f"  failure: {f['actual']:.4e} > {f['limit']:.0e} ({f['detail']})")
        lines.append("")
    save_item("gates", out)
    (DATA / "gates_ctest_summary.txt").write_text("\n".join(lines).rstrip() + "\n")
    for label, r in out.items():
        print(f"{label}: passed {r['passed']}, failed {r['failed']}, skipped {r['skipped']}",
              [f"{f['actual']:.3e}" for f in r["failures"]])


if __name__ == "__main__":
    main()

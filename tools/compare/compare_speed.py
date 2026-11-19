#!/usr/bin/env python3
"""Compares the speed and host peak memory of two spodvr builds from the logs of alternating runs
(tools/benchmark/compare_speed.sh): each log holds spodvr's output ("Propagation time: X s") and
/usr/bin/time -v's ("Maximum resident set size (kbytes): K", the host's resident memory; GPU memory is not
measured). The new build passes when its median propagation time is at most max(1 %, spread) slower than the
base build's, the spread being the larger of the two builds' (max - min) / median, and its largest host peak
memory at most 16 MiB above the base's. A faster or smaller new build passes.

Usage: python3 -I -B tools/compare/compare_speed.py --base LOG... --new LOG...
Output: the medians, spreads and peak memory of both builds, then "OK: ..." or one line per failed bound
("SLOWER: ...", "MORE MEMORY: ...").
Exit status: 0 within the bounds; 1 outside them, or a log that cannot be read ("error: <what>; <fix>" on
stderr); 2 bad arguments.
Environment: none read. Run time: well under a second.
Docs: docs/CONTRIBUTING.md, section "Additions to the time loop".
"""

import argparse
import re
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path

MIN_ALLOWED_SLOWDOWN = 0.01  # the time bound is at least 1 % of the base median, however quiet the runs
MEMORY_SLACK_KIB = 16 * 1024  # 16 MiB
PROPAGATION_TIME_LINE = re.compile(r"Propagation time:\s*([0-9.]+)\s*s")
PEAK_MEMORY_LINE = re.compile(r"Maximum resident set size \(kbytes\):\s*([0-9]+)")


class Refusal(Exception):
    """A log that cannot be read, printed as "error: <what>; <fix>" with exit status 1."""


@dataclass
class Runs:
    """The propagation times (s) and host peak memory (KiB) of the runs of one build."""

    times_s: list[float]
    peaks_kib: list[int]

    @property
    def median_s(self) -> float:
        """The median propagation time."""
        return statistics.median(self.times_s)

    @property
    def spread(self) -> float:
        """(max - min) / median of the propagation times."""
        return (max(self.times_s) - min(self.times_s)) / self.median_s


def read_log(path: str) -> tuple[float, int]:
    """(propagation time in s, host peak memory in KiB) of one run's log."""
    try:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
    except OSError as error:
        raise Refusal(f"{path}: cannot be read ({error.strerror}); give the logs compare_speed.sh wrote") from error
    time, peak = PROPAGATION_TIME_LINE.search(text), PEAK_MEMORY_LINE.search(text)
    if not time:
        raise Refusal(f"{path}: no 'Propagation time' line; the run did not finish, see the log")
    if not peak:
        raise Refusal(f"{path}: no peak memory line; run spodvr under /usr/bin/time -v")
    if float(time.group(1)) <= 0.0:
        raise Refusal(f"{path}: a propagation time of 0 s; run more steps")
    return float(time.group(1)), int(peak.group(1))


def read_runs(paths: list[str]) -> Runs:
    """The Runs of one build's logs."""
    values = [read_log(p) for p in paths]
    return Runs([t for t, _ in values], [m for _, m in values])


def verdict(base: Runs, new: Runs) -> tuple[bool, list[str]]:
    """(passed, lines): the new build at most max(1 %, spread) slower than the base's median time and at most
    16 MiB above its peak memory; one line per failed bound, or one "OK" line."""
    allowed = max(MIN_ALLOWED_SLOWDOWN, base.spread, new.spread)
    change = (new.median_s - base.median_s) / base.median_s
    extra_kib = max(new.peaks_kib) - max(base.peaks_kib)
    failed = []
    if change > allowed:
        failed.append(f"SLOWER: median {change:+.2%}, more than the allowed {allowed:.2%}")
    if extra_kib > MEMORY_SLACK_KIB:
        failed.append(
            f"MORE MEMORY: peak {extra_kib / 1024:+.1f} MiB, more than the allowed {MEMORY_SLACK_KIB // 1024}"
        )
    if failed:
        return False, failed
    return True, [f"OK: median {change:+.2%} (allowed +{allowed:.2%}), peak memory {extra_kib / 1024:+.1f} MiB"]


def main(argv: list[str]) -> int:
    """Prints both builds' numbers and the verdict; returns the exit status."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", nargs="+", required=True, help="the base build's run logs")
    parser.add_argument("--new", nargs="+", required=True, help="the new build's run logs")
    args = parser.parse_args(argv)
    try:
        base, new = read_runs(args.base), read_runs(args.new)
    except Refusal as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    for name, runs in (("base", base), ("new", new)):
        print(
            f"{name}: median {runs.median_s:.2f} s over {len(runs.times_s)} runs, spread {runs.spread:.2%}, "
            f"host peak memory {max(runs.peaks_kib) / 1024:.1f} MiB"
        )
    passed, lines = verdict(base, new)
    print("\n".join(lines))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

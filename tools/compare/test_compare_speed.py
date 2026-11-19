"""Tests of compare_speed.py (the speed and memory comparison of rule "time-loop additions are free when off")
on synthetic run logs: the bound max(1 %, spread) one-sided, the memory slack, and the refusals."""

import pytest

import compare_speed


def log(tmp_path, name, seconds, peak_kib):
    """A run log with spodvr's propagation-time line and /usr/bin/time -v's peak-memory line."""
    path = tmp_path / name
    path.write_text(
        f"  Propagation time: {seconds} s (0:30, 50 steps)\n\tMaximum resident set size (kbytes): {peak_kib}\n"
    )
    return path


def logs(tmp_path, prefix, seconds, peak_kib=1_000_000):
    """One log per propagation time."""
    return [log(tmp_path, f"{prefix}{k}.log", s, peak_kib) for k, s in enumerate(seconds)]


def compare(tmp_path, base_s, new_s, base_kib=1_000_000, new_kib=1_000_000):
    """The verdict of new against base (each a list of propagation times)."""
    base = compare_speed.read_runs(logs(tmp_path, "b", base_s, base_kib))
    new = compare_speed.read_runs(logs(tmp_path, "n", new_s, new_kib))
    return compare_speed.verdict(base, new)


def test_equal_builds_pass(tmp_path):
    """The same times (spread 0.67 %): within the bound; main exits 0."""
    base, new = logs(tmp_path, "b", [30.0, 30.1, 29.9]), logs(tmp_path, "n", [30.0, 30.2, 29.8])
    assert compare_speed.main(["--base", *map(str, base), "--new", *map(str, new)]) == 0


def test_slower_than_the_bound_fails(tmp_path):
    """Negative test: quiet runs (spread 0.1 %), so the bound is the 1 % floor: 1.2 % slower fails with SLOWER
    (main exits 1), 0.8 % slower passes."""
    quiet = [30.0, 30.03, 30.01]
    passed, lines = compare(tmp_path, quiet, [30.36, 30.39, 30.37])
    assert not passed and lines == [lines[0]] and lines[0].startswith("SLOWER")
    assert compare(tmp_path, quiet, [30.24, 30.27, 30.25])[0]
    base, new = logs(tmp_path, "b2", quiet), logs(tmp_path, "n2", [30.36, 30.39, 30.37])
    assert compare_speed.main(["--base", *map(str, base), "--new", *map(str, new)]) == 1


def test_noisy_runs_widen_the_bound(tmp_path):
    """2 % slower passes when either build's spread is 5 % (the spread is the noise), fails when neither is."""
    assert compare(tmp_path, [30.0, 31.5, 30.2], [30.6, 30.8, 30.7])[0]
    assert compare(tmp_path, [30.0, 30.1, 30.05], [30.6, 32.1, 30.7])[0]
    assert not compare(tmp_path, [30.0, 30.1, 30.05], [30.6, 30.8, 30.7])[0]


def test_faster_build_passes(tmp_path):
    """A new build 10 % faster passes: the bound is on the extra cost only."""
    assert compare(tmp_path, [30.0, 30.03, 30.01], [27.0, 27.02, 27.01])[0]


def test_more_memory_fails(tmp_path):
    """17 MiB more peak memory fails, 15 MiB passes; slower and larger together report both."""
    plus_17, plus_15 = 1_000_000 + 17 * 1024, 1_000_000 + 15 * 1024
    passed, lines = compare(tmp_path, [30.0, 30.0], [30.0, 30.0], new_kib=plus_17)
    assert not passed and lines[0].startswith("MORE MEMORY")
    assert compare(tmp_path, [30.0, 30.0], [30.0, 30.0], new_kib=plus_15)[0]
    passed, lines = compare(tmp_path, [30.0, 30.0], [33.0, 33.0], new_kib=plus_17)
    assert [line.split(":")[0] for line in lines] == ["SLOWER", "MORE MEMORY"]


def test_unreadable_logs_are_refused(tmp_path, capsys):
    """A log without the propagation time, one without the memory line and a missing log: refused, main exits 1
    with "error:" on stderr naming the file."""
    no_time = tmp_path / "no_time.log"
    no_time.write_text("\tMaximum resident set size (kbytes): 5\n")
    no_memory = tmp_path / "no_memory.log"
    no_memory.write_text("  Propagation time: 30.0 s\n")
    good = logs(tmp_path, "g", [30.0])
    for bad in (no_time, no_memory, tmp_path / "missing.log"):
        with pytest.raises(compare_speed.Refusal):
            compare_speed.read_log(bad)
        assert compare_speed.main(["--base", str(bad), "--new", *map(str, good)]) == 1
        assert capsys.readouterr().err.startswith(f"error: {bad}")

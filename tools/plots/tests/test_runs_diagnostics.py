"""The verdict on the fixture runs: each raises the findings its state calls for. The fixtures are small and
not fully converged, and their bases are small (dailight_oblique: j' = JMax = 2 and the outermost diffraction
channels carry flux; both inherit gate-case r grids whose gas-phase states reach the grid ends); the stopped
run has only the files a run that stopped early has."""

from pathlib import Path

import pytest

from readers.segments import joined_info, ordered_segments
from verdict import run_verdict

EXPECTED = {
    "hardwall_tiny": [
        "unanalysed flux, upper bound",
        "P_scattered still moving",
        "unitarity error",
        "gas-phase states reach the ends of the r grid",
        "slowest converged energy",
        "norm left",
    ],
    "dailight_oblique": [
        "P_scattered still moving",
        "unitarity error",
        "flux at the edge of the basis: j' = 2",
        "flux at the edge of the basis: outermost diffraction channels",
        "gas-phase states reach the ends of the r grid",
        "slowest converged energy",
        "norm left",
    ],
    "potential_only": ["not checked: the propagation"],
}


def verdict_of(runs: list[Path]) -> list:
    """The verdict of a run given by its segment directories (any order)."""
    infos, _ = ordered_segments(runs)
    return run_verdict([i.directory for i in infos], joined_info(infos))


@pytest.mark.parametrize("run", sorted(EXPECTED))
def test_findings_of_each_fixture(fixture_runs: Path, run: str) -> None:
    """The names of the findings, worst first."""
    assert [f.name for f in verdict_of([fixture_runs / run])] == EXPECTED[run]


def test_stopped_run_is_incomplete(stopped_run: Path) -> None:
    """A run that stopped: one error naming the missing files and how far the BK file got (its cut last
    record left out), and no other check."""
    (finding,) = verdict_of([stopped_run])
    assert (finding.severity, finding.name) == ("error", "run incomplete")
    assert finding.value == "averages.dat, scattering_summary.tsv, s_matrix.tsv missing; BK records up to t = 500 a.u."


def test_split_run_has_the_unsplit_verdict(fixture_runs: Path) -> None:
    """The two segments of hardwall_tiny_split, joined, give the findings and values of hardwall_tiny."""
    split = fixture_runs / "hardwall_tiny_split"
    joined = verdict_of([split / "second", split / "first"])
    whole = verdict_of([fixture_runs / "hardwall_tiny"])
    assert [(f.name, f.value) for f in joined] == [(f.name, f.value) for f in whole]


def test_missing_checkpoint_is_not_checked(fixture_runs: Path, tmp_path: Path) -> None:
    """A run without checkpoint.bin (CheckpointFormat "fortran"): the unanalysed flux is "not checked"."""
    run = tmp_path / "no_checkpoint"
    run.mkdir()
    for path in (fixture_runs / "hardwall_tiny").iterdir():
        if path.is_file() and path.name != "checkpoint.bin":
            (run / path.name).symlink_to(path.resolve())
    names = [f.name for f in verdict_of([run])]
    assert "not checked: unanalysed flux" in names and "unanalysed flux, upper bound" not in names

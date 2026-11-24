"""The verdict on the fixture runs: each raises the findings its state calls for. The fixtures are small and
not fully converged, and their bases are small (dailight_oblique: j' = JMax = 2 and the outermost diffraction
channels carry flux; both inherit gate-case r grids whose gas-phase states reach the grid ends); the stopped
run has only the files a run that stopped early has."""

import copy
import shutil
from dataclasses import replace
from pathlib import Path

import pytest

from readers.binary import read_bk_coefficients
from readers.run_info import RunInfo, read_run_info
from readers.segments import joined_info, ordered_segments
from readers.tables import read_averages
from verdict import absorber_reflection_finding, frames_at_edges, run_verdict

EXPECTED = {
    "hardwall_tiny": [
        "unanalysed flux, upper bound",
        "P_scattered still moving",
        "unitarity error",
        # 2.36e-3 in P at 0.22 eV: the 3.52 bohr ZSpOpt has R 4.86e-7 at 0.2 eV (section 7's target: R + T <= 4e-9)
        "absorber reflection",
        "gas-phase states reach the ends of the r grid",
        "slowest converged energy",
        "norm left",
        "norm by region at the end",
    ],
    "dailight_oblique": [
        "P_scattered still moving",
        "unitarity error",
        "absorber reflection",  # 2.56e-3 in P at 0.30 eV (P 0.65): the same 3.52 bohr ZSpOpt
        "flux at the edge of the basis: j' = 2",
        "flux at the edge of the basis: outermost diffraction channels",
        "gas-phase states reach the ends of the r grid",
        # frames/ with FrameExtras: 13 % of the norm in the outer orders (n or m = 2) of the 4 x 4 (x, y) grid
        "packet at the edge of the basis: outer diffraction orders (momentum grid, above the first frame's)",
        "norm budget does not close",  # norm_budget/: the angular transform at the surface (JMax 2)
        "slowest converged energy",
        "packet at the edge of the basis: top J",  # frames/, information (JMax 2: J = 2 holds up to half the norm)
        "norm left",
        "norm by region at the end",
    ],
    "potential_only": ["not checked: the propagation"],
    "pes_cut": ["not checked: the propagation"],  # a write2DPES run (run mode 0) stops after the cut
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


def test_damaged_densities_are_a_finding(fixture_runs: Path, tmp_path: Path) -> None:
    """A complete run whose densities/manifest.json is damaged: the verdict says the regions were not checked
    instead of failing; the other findings stay."""
    run = tmp_path / "hardwall_tiny"
    shutil.copytree(fixture_runs / "hardwall_tiny", run)
    (run / "densities" / "manifest.json").write_text("{", encoding="utf-8")
    names = [f.name for f in verdict_of([run])]
    assert "not checked: norm by region" in names and "norm by region at the end" not in names
    assert "unitarity error" in names


def test_norm_by_region_lists_the_rest_norm(fixture_runs: Path) -> None:
    """The specular wave inside its analysis line, listed last, equals the rest norm of averages.dat at the same
    time (the program's sum over the same cells)."""
    finding = next(f for f in verdict_of([fixture_runs / "dailight_oblique"]) if f.name == "norm by region at the end")
    listed = float(finding.value.rsplit(": ", 1)[1])
    rest_norm = read_averages(fixture_runs / "dailight_oblique" / "averages.dat").rest_norm[-1]
    assert abs(listed - rest_norm) <= 0.05 * rest_norm  # printed with 2 digits


def test_unreadable_frames_are_not_checked(fixture_runs: Path, tmp_path: Path) -> None:
    """A run whose frames/manifest.json is not JSON: the basis edges of frames/ are "not checked" (with the reader's
    reason), not an exception; the intact run gives findings or none, never that."""
    fixture = fixture_runs / "dailight_oblique"
    run = tmp_path / "broken_frames"
    run.mkdir()
    for path in fixture.iterdir():
        if path.name != "frames":
            (run / path.name).symlink_to(path.resolve())
    (run / "frames").mkdir()
    (run / "frames" / "manifest.json").write_text("{", encoding="utf-8")
    (finding,) = frames_at_edges([run])
    assert finding.name == "not checked: the basis edges of frames/" and "not JSON" in finding.value
    assert all(f is None or not f.name.startswith("not checked") for f in frames_at_edges([fixture]))


def test_absorber_reflection_checks_and_an_absorber_that_is_off(fixture_runs: Path) -> None:
    """hardwall_tiny with its run_info changed: an absorber power other than 2 and a missing BK file are "not
    checked"; ZSpOpt off (strength 0) returns a = 1 and warns; a low-Z absorber in the list is named."""
    run = fixture_runs / "hardwall_tiny"
    info, bk = read_run_info(run), read_bk_coefficients(run / "bk_coefficients.bin")

    def changed(change: dict) -> RunInfo:
        """The run's RunInfo with its absorbers changed (a deep copy)."""
        data = copy.deepcopy(info.data)
        data["absorbers"].update(change)
        return replace(info, data=data)

    power = absorber_reflection_finding(run, changed({"power": 3}), bk)
    assert (power.name, power.value) == ("not checked: absorber reflection", "power 3, not 2")
    assert absorber_reflection_finding(run, info, None).value == "no bk_coefficients.bin or s_matrix.tsv"
    absorbers = copy.deepcopy(info.data["absorbers"]["list"])
    next(a for a in absorbers if a["name"] == "ZSpOpt")["strength"] = 0.0
    off = absorber_reflection_finding(run, changed({"list": absorbers}), bk)
    assert off.severity == "warning" and float(off.value.split()[2]) > 1.0  # "up to <error> in P": 2 a P with a = 1
    low_z = absorbers + [{"name": "ZOptLow", "start": 0.0, "end": 1.0, "strength": 0.01}]
    assert "low-Z absorbers on: not estimated" in absorber_reflection_finding(run, changed({"list": low_z}), bk).value

"""The command line on the fixture runs: list, plot, report and verdict, their exit statuses and files."""

import json
import subprocess
import sys
from pathlib import Path

from explanations import FIGURES, TITLES

TOOL = Path(__file__).resolve().parent.parent / "spodvr_plot.py"


def run(*arguments: object) -> subprocess.CompletedProcess:
    """spodvr_plot.py with the arguments, as a program."""
    command = [sys.executable, "-I", "-B", str(TOOL), *map(str, arguments)]
    return subprocess.run(command, capture_output=True, text=True, check=False)


def test_list_names_every_figure() -> None:
    """`list` prints one line per figure."""
    result = run("list")
    assert result.returncode == 0
    for name in ("verdict", "results", "averages", "gas_phase", "bk"):
        assert any(line.startswith(name) for line in result.stdout.splitlines()), name


def listing(directory: Path) -> dict[str, int]:
    """Every file below a directory with its size (to see that a command wrote nothing there)."""
    return {str(p.relative_to(directory)): p.stat().st_size for p in directory.rglob("*") if p.is_file()}


def test_report_writes_only_into_out(fixture_runs: Path, tmp_path: Path) -> None:
    """`report` writes exactly report.md, its manifest and the files the manifest lists into --out (a movie as
    GIF and PNG, the GIF shown; the frame movies of the two sites under one heading with the explanation once),
    nothing into the run directory; a second report replaces its own files. On the smaller hard-wall run, since the
    report is made twice."""
    run_dir, out = fixture_runs / "hardwall_tiny", tmp_path / "report"
    before = listing(run_dir)
    assert run("report", run_dir, "--out", out).returncode == 0
    manifest = json.loads((out / "written_by_spodvr_plot.json").read_text(encoding="utf-8"))["files"]
    assert {p.name for p in out.iterdir()} == {*manifest, "written_by_spodvr_plot.json"}
    assert {"report.md", "results.png", "results.pdf", "bk.png"} <= set(manifest) and "verdict.png" not in manifest
    assert {"density_movie.gif", "density_movie.png"} <= set(manifest) and "density_movie.pdf" not in manifest
    report = (out / "report.md").read_text(encoding="utf-8")
    assert "## Verdict" in report and "(density_movie.gif)" in report
    frame_section = report.split(f"## {TITLES['frame_movie']}\n")[1].split("\n## ")[0]
    assert report.count(f"## {TITLES['frame_movie']}\n") == 1 and frame_section.count(FIGURES["frame_movie"]) == 1
    assert "(frame_movie_site1.gif)" in frame_section and "(frame_movie_site2.gif)" in frame_section
    assert run("report", run_dir, "--out", out).returncode == 0
    assert listing(run_dir) == before


def test_files_the_tool_did_not_write_stop_it_before_writing(fixture_runs: Path, tmp_path: Path) -> None:
    """A results.png or report.md the tool did not write: `report` and `plot` exit 2, the file stays
    untouched, and nothing is written."""
    for command, foreign in (("report", "results.png"), ("report", "report.md"), ("plot", "results.png")):
        out = tmp_path / f"{command}_{foreign}"
        out.mkdir()
        (out / foreign).write_text("mine", encoding="utf-8")
        arguments = ("results",) if command == "plot" else ()
        result = run(command, *arguments, fixture_runs / "dailight_oblique", "--out", out)
        assert result.returncode == 2 and "not written by this tool" in result.stderr
        assert listing(out) == {foreign: 4}


def test_report_of_a_run_without_results(fixture_runs: Path, stopped_run: Path, tmp_path: Path) -> None:
    """A potential-only run and a run that stopped: exit 0, the verdict says why, and the figures the run's
    files do not allow are listed as not drawn."""
    for run_dir, verdict in ((fixture_runs / "potential_only", "not checked"), (stopped_run, "run incomplete")):
        out = tmp_path / run_dir.name
        assert run("report", run_dir, "--out", out).returncode == 0
        report = (out / "report.md").read_text(encoding="utf-8")
        assert verdict in report and "## Not drawn" in report


def test_report_of_a_run_without_bk_records(fresh_run: Path, tmp_path: Path) -> None:
    """A run before its first analysis step: the report is written, says the BK figure could not be drawn and
    exits 1; the verdict exits 0 and says there is no BK record."""
    result = run("report", fresh_run, "--out", tmp_path / "report")
    assert result.returncode == 1 and "no BK record yet" in result.stderr
    assert "no BK record yet" in (tmp_path / "report" / "report.md").read_text(encoding="utf-8")
    verdict = run("verdict", fresh_run)
    assert verdict.returncode == 0 and "no BK record" in verdict.stdout


def test_out_that_is_a_file_exits_2(fixture_runs: Path, tmp_path: Path) -> None:
    """--out naming an existing file: exit 2 with "error: …", no traceback."""
    (tmp_path / "a_file").write_text("", encoding="utf-8")
    result = run("plot", "gas_phase", fixture_runs / "potential_only", "--out", tmp_path / "a_file")
    assert result.returncode == 2 and result.stderr.startswith("error: ")


def test_verdict_and_plot(fixture_runs: Path, stopped_run: Path, tmp_path: Path) -> None:
    """`verdict` prints the findings and exits 0 even with an error finding; `plot` writes one figure, and
    refuses (exit 2) a figure whose files the run lacks."""
    result = run("verdict", stopped_run)
    assert result.returncode == 0 and result.stdout.startswith("ERROR")
    assert run("plot", "gas_phase", fixture_runs / "potential_only", "--out", tmp_path).returncode == 0
    assert run("plot", "gas_phase", fixture_runs / "potential_only", "--out", tmp_path).returncode == 0  # its own
    assert (tmp_path / "gas_phase.png").is_file()
    assert run("plot", "results", fixture_runs / "potential_only", "--out", tmp_path).returncode == 2


def test_unreadable_run_exits_2(tmp_path: Path) -> None:
    """A directory without run_info.json: exit 2 with "error: …"."""
    result = run("verdict", tmp_path)
    assert result.returncode == 2 and result.stderr.startswith("error: ")


def test_options_reach_the_figures(fixture_runs: Path, tmp_path: Path) -> None:
    """--energy picks the analysis energy of the one-energy panels (its title says which); --channels sets the
    number of delay curves; a --channels below 1 is refused (exit 2)."""
    run_dir = fixture_runs / "dailight_oblique"
    assert run("plot", "state_to_state", run_dir, "--out", tmp_path / "a", "--energy", "0.62").returncode == 0
    assert run("plot", "time_delay", run_dir, "--out", tmp_path / "b", "--channels", "3").returncode == 0
    assert run("plot", "time_delay", run_dir, "--out", tmp_path / "c", "--channels", "0").returncode == 2
    assert (tmp_path / "a" / "state_to_state.png").is_file() and (tmp_path / "b" / "time_delay.png").is_file()

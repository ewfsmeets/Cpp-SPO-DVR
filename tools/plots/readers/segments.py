"""A run split into restart segments, read as one run: the segments in time order, and their records joined.

Joining is concatenation: every joined quantity is a value at its own time, never a running sum. A segment
must start where the previous one ended (the program's tolerant comparison, 1e-12); a record at that shared
time is dropped from the later segment after checking that it equals the earlier one to REPEAT_TOLERANCE (a
restart on the same device and build repeats it bitwise; on another device the reductions round differently).
(A restart's bk_coefficients.bin already holds every earlier record: read the last segment's, never join them.)
"""

from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import numpy as np

from copies.real_compare import rne
from readers.arrays import Series, read_series
from readers.errors import ReadError
from readers.run_info import RunInfo, read_run_info
from readers.tables import Averages, first_non_finite, read_averages

# A repeated record's tolerance, relative to the record's largest |value|: the CPU's and the GPU's densities of
# the same wave function differ by <= 1e-14 per cell (CudaBackend.DensityRecordsOfTheSameWaveMatchTheCpu). For an
# averages.dat row the largest value is <Z> or <r> (bohr), so its norm is checked to about 1e-11 absolute.
REPEAT_TOLERANCE = 1e-12


def ordered_segments(directories: list[Path]) -> tuple[list[RunInfo], list[str]]:
    """The segments' run_info in time order, and warnings (a segment whose RestartDir does not name the
    previous segment's directory: directories get moved and RestartDir may be relative, so only the last
    path component is compared, and a mismatch is not an error).

    Raises ReadError when a segment does not start at the previous one's tEnd.
    """
    infos = sorted((read_run_info(d) for d in directories), key=lambda info: info.t_start_au)
    warnings = []
    for earlier, later in pairwise(infos):
        if rne(later.t_start_au, earlier.t_end_au):
            raise ReadError(
                f"{later.directory}: starts at {later.t_start_au}, but {earlier.directory} ends at "
                f"{earlier.t_end_au}; give the segments of one run"
            )
        restart = later.data.get("restart_dir")
        if restart is None or Path(restart).name != earlier.directory.resolve().name:
            warnings.append(f"{later.directory}: RestartDir {restart} is not {earlier.directory}")
    return infos, warnings


def joined_info(infos: list[RunInfo]) -> RunInfo:
    """The run_info of the whole run: the last segment's, with the first segment's start time."""
    return replace(infos[-1], t_start_au=infos[0].t_start_au)


def dropped_repeats(times: list[np.ndarray]) -> list[int]:
    """Per segment, how many of its first records the join drops: 1 when its first time equals the last time kept
    so far (it repeats that record), else 0."""
    dropped, last_time = [], None
    for t in times:
        repeats = len(t) > 0 and last_time is not None and not rne(t[0], last_time)
        dropped.append(1 if repeats else 0)
        if len(t) > dropped[-1]:
            last_time = t[-1]
    return dropped


def repeats_record(record: np.ndarray, kept: np.ndarray) -> bool:
    """Whether `record` equals `kept` to REPEAT_TOLERANCE of kept's largest |value| (NaN in either: False)."""
    difference = np.max(np.abs(record - kept), initial=0.0)
    return bool(difference <= REPEAT_TOLERANCE * np.max(np.abs(kept), initial=0.0))


def check_repeats(times: list[np.ndarray], values: list[np.ndarray], dropped: list[int], what: str) -> None:
    """ReadError naming `what` when a dropped record does not repeat the last record kept before it (to
    REPEAT_TOLERANCE), or when it holds a NaN or inf (then the propagation diverged before the split)."""
    last_value = None
    for t, v, skip in zip(times, values, dropped, strict=True):
        if skip and not np.all(np.isfinite(v[0])):
            raise ReadError(f"{what}: the record at t = {t[0]} holds a NaN or inf; the run diverged before the restart")
        if skip and not repeats_record(np.asarray(v[0]), np.asarray(last_value)):
            raise ReadError(f"{what}: the record at t = {t[0]} differs between two segments; they are not one run")
        if len(v) > skip:
            last_value = v[-1]


def join_records(times: list[np.ndarray], values: list[np.ndarray], what: str) -> tuple[np.ndarray, np.ndarray]:
    """The records of consecutive segments (times and values along axis 0) as one series: a segment's first record
    that repeats the last one kept is dropped (dropped_repeats) after check_repeats. Returns: (times, values)."""
    dropped = dropped_repeats(times)
    check_repeats(times, values, dropped, what)
    return kept(times, dropped), kept(values, dropped)


def kept(arrays: list[np.ndarray], dropped: list[int]) -> np.ndarray:
    """The segments' arrays without their dropped first records, concatenated along axis 0."""
    return np.concatenate([a[skip:] for a, skip in zip(arrays, dropped, strict=True)])


def read_joined_averages(directories: list[Path]) -> tuple[Averages, list[float]]:
    """averages.dat of the segments of one run, joined, and the start times of the segments (a.u.)."""
    infos, _ = ordered_segments(directories)
    parts = [read_averages(info.directory / "averages.dat") for info in infos]
    columns = [np.column_stack([a.norm, a.z_mean_bohr, a.r_mean_bohr, a.rest_norm]) for a in parts]
    times, rows = join_records([a.times_au for a in parts], columns, "averages.dat")
    norm, z_mean_bohr, r_mean_bohr, rest_norm = rows.T
    joined = Averages(times, norm, z_mean_bohr, r_mean_bohr, rest_norm, first_non_finite(rows))
    return joined, [info.t_start_au for info in infos]


def joined_first_non_finite(first_bad: list[int | None], lengths: list[int], dropped: list[int]) -> int | None:
    """The joined series' first non-finite record from each segment's (check_repeats has refused a dropped record
    that is not finite, so a segment's first non-finite record is never a dropped one)."""
    offset = 0
    for bad, length, skip in zip(first_bad, lengths, dropped, strict=True):
        if bad is not None:
            assert bad >= skip, "a non-finite dropped record should have been refused"
            return offset + bad - skip
        offset += length - skip
    return None


def read_joined_series(directories: list[Path], name: str) -> tuple[Series, list[float]]:
    """The series directory `name` (e.g. "densities") of the segments of one run, joined record by record
    (a record repeated at a segment boundary must equal the kept one, see check_repeats, and appears once), and
    the segments' start times. One segment is returned as read (memory-mapped, not copied). Raises ReadError
    when the segments' manifests differ."""
    infos, _ = ordered_segments(directories)
    parts = [read_series(info.directory / name) for info in infos]
    starts = [info.t_start_au for info in infos]
    if len(parts) == 1:
        return parts[0], starts
    for part in parts[1:]:
        if part.manifest != parts[0].manifest:
            raise ReadError(
                f"{part.directory}/manifest.json: differs from {parts[0].directory}/manifest.json (datasets, grids, "
                "record interval or descriptions); give the segments of one run"
            )
    times = [p.times_au for p in parts]
    dropped = dropped_repeats(times)
    what = f"{name}/ of {', '.join(str(info.directory) for info in infos)}"
    arrays = {}
    for dataset in parts[0].arrays:
        values = [p.arrays[dataset] for p in parts]
        check_repeats(times, values, dropped, what)
        arrays[dataset] = kept(values, dropped)
    first_bad = joined_first_non_finite([p.first_non_finite for p in parts], [len(t) for t in times], dropped)
    return replace(parts[0], times_au=kept(times, dropped), arrays=arrays, first_non_finite=first_bad), starts

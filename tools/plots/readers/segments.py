"""A run split into restart segments, read as one run: the segments in time order, and their records joined.

Joining is concatenation: every joined quantity is a value at its own time, never a running sum. A segment
must start where the previous one ended (the program's tolerant comparison, 1e-12); a record at that shared
time is dropped from the later segment after checking that it equals the earlier one bitwise. (A restart's
bk_coefficients.bin already holds every earlier record: read the last segment's, never join them.)
"""

from dataclasses import replace
from itertools import pairwise
from pathlib import Path

import numpy as np

from copies.real_compare import rne
from readers.errors import ReadError
from readers.run_info import RunInfo, read_run_info
from readers.tables import Averages, first_non_finite, read_averages


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


def join_records(times: list[np.ndarray], values: list[np.ndarray], what: str) -> tuple[np.ndarray, np.ndarray]:
    """The records of consecutive segments (times and values along axis 0) as one series.

    A segment whose first time equals the previous segment's last time repeats that record: it is dropped
    after checking that it is bitwise equal (else ReadError). Returns: (times, values).
    """
    all_times = [times[0]]
    all_values = [values[0]]
    for t, v in zip(times[1:], values[1:], strict=True):
        repeats_last_record = len(t) > 0 and len(all_times[-1]) > 0 and not rne(t[0], all_times[-1][-1])
        if repeats_last_record:
            if v[0].tobytes() != all_values[-1][-1].tobytes():  # bitwise, also NaN and the sign of zero
                raise ReadError(f"{what}: the record at t = {t[0]} differs between two segments; they are not one run")
            t, v = t[1:], v[1:]
        all_times.append(t)
        all_values.append(v)
    return np.concatenate(all_times), np.concatenate(all_values)


def read_joined_averages(directories: list[Path]) -> tuple[Averages, list[float]]:
    """averages.dat of the segments of one run, joined, and the start times of the segments (a.u.)."""
    infos, _ = ordered_segments(directories)
    parts = [read_averages(info.directory / "averages.dat") for info in infos]
    columns = [np.column_stack([a.norm, a.z_mean_bohr, a.r_mean_bohr, a.rest_norm]) for a in parts]
    times, rows = join_records([a.times_au for a in parts], columns, "averages.dat")
    norm, z_mean_bohr, r_mean_bohr, rest_norm = rows.T
    joined = Averages(times, norm, z_mean_bohr, r_mean_bohr, rest_norm, first_non_finite(rows))
    return joined, [info.t_start_au for info in infos]

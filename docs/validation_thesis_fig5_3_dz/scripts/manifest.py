#!/usr/bin/env python3
"""sha256 manifest of the raw files the report is based on (heavy: reads ~20 GB): every file of the
three run directories and of the strict-build GPU run, their raw logs (uncompressed;
data/logs/*.log.gz hold them), the onlyAnalysis outputs, the executables and the thesis PDF. Writes data/sha256_manifest.txt."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path

from common import DATA, RAW_LOGS, REPO, RUNS, RUN_DIRS

STRICT = RUNS / "thesis_fig5_3_dz_cuda"          # strict-build GPU run (strict_gpu.py)
EXTRA = [REPO / "build-fast" / "spodvr", REPO / "build-cuda-fast" / "spodvr", REPO / "build-cuda" / "spodvr",
         Path.home() / "work/SPO-DVR_v5.19/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x",
         Path.home() / "work/mf_somers_phd_thesis.pdf",
         REPO / "inputs/thesis_fig5_3_dz/input.json", REPO / "inputs/thesis_fig5_3_dz/DVR.inp",
         # the reference executable's objects linked by the angular-matrix driver (2026-11-08)
         *[Path.home() / f"work/SPO-DVR_v5.19/procfiles.dir/OBJS.dir/{n}.o" for n in
           ("transform", "angulargrid", "legendre", "ffts", "fftw_2_scs", "constants", "errortrap", "MPI_wrapper", "x86_hack")]]


def sha256(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def show(p: Path) -> str:
    s = str(p)
    home = str(Path.home())
    if s.startswith(str(REPO) + "/"):
        return s[len(str(REPO)) + 1:]
    return s.replace(home, "~")


def main():
    files = []
    for d in RUN_DIRS.values():
        files += sorted(p for p in d.iterdir() if p.is_file())
    files += sorted(p for p in STRICT.iterdir() if p.is_file())
    files += list(RAW_LOGS.values()) + [RUNS / "thesis_fig5_3_dz_cuda.log"]
    for k in ("cpu", "gpu", "fortran"):
        d = RUNS / f"thesis_fig5_3_dz_oa_{k}" / "fine_analysis"
        files += sorted(p for p in d.rglob("*") if p.is_file())
    files += EXTRA
    lines = ["# sha256 of the raw files behind this report (scripts/manifest.py). Sizes in bytes.",
             "# data/logs/*.log.gz are `gzip -9 -n` of the raw logs listed here (sha256 of the uncompressed log)."]
    for p in files:
        lines.append(f"{sha256(p)}  {os.path.getsize(p):>12d}  {show(p)}")
    (DATA / "sha256_manifest.txt").write_text("\n".join(lines) + "\n")
    print(len(files), "files hashed")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Phase 0: provenance of the three runs (builds, binaries, commits, input files, hardware,
timings from /usr/bin/time and the programs' own timers). Writes data/items/phase0.json.

Needs the repository, the build directories and runs/ (extract mode)."""
from __future__ import annotations

import datetime as dt
import json
import os
import re
import subprocess
from pathlib import Path

from common import RAW_LOGS, REPO, RUN_DIRS, INPUT_DIR, clean_log_lines, save_item

BUILDS = {"cpu": REPO / "build-fast", "gpu": REPO / "build-cuda-fast"}
FORTRAN_EXE = Path.home() / "work/SPO-DVR_v5.19/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x"
CACHE_KEYS = ["CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CUDA_COMPILER", "SPODVR_FAST_MATH", "SPODVR_ARCH",
              "SPODVR_USE_CUDA", "SPODVR_USE_MKL", "SPODVR_FFT_LIBRARY", "SPODVR_PGO", "MKL_INCLUDE"]


def sh(*cmd, cwd=REPO) -> str:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, check=False).stdout.strip()


def mtime(path: Path) -> str:
    return dt.datetime.fromtimestamp(os.path.getmtime(path)).astimezone().isoformat(timespec="seconds")


def cmake_cache(build: Path) -> dict:
    out = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        m = re.match(r"^([A-Z_]+):[A-Z]+=(.*)$", line)
        if m and m.group(1) in CACHE_KEYS:
            out[m.group(1)] = m.group(2)
    return out


def script_header(raw: Path) -> dict:
    with open(raw, "rb") as f:
        head = f.read(600).decode(errors="replace")
    with open(raw, "rb") as f:
        f.seek(max(0, os.path.getsize(raw) - 400))
        tail = f.read().decode(errors="replace")
    start = re.search(r"Script started on (\S+ \S+)", head).group(1)
    command = re.search(r'COMMAND="([^"]*)"', head).group(1)
    done = re.search(r"Script done on (\S+ \S+) \[COMMAND_EXIT_CODE=\"(\d+)\"\]", tail)
    return {"start": start, "command": " ".join(command.split()), "end": done.group(1), "exit_code": int(done.group(2))}


def time_v(lines: list[str]) -> dict:
    keys = {"User time (seconds)": "user_s", "System time (seconds)": "system_s",
            "Percent of CPU this job got": "cpu_percent", "Elapsed (wall clock) time (h:mm:ss or m:ss)": "wall",
            "Maximum resident set size (kbytes)": "max_rss_kb", "Exit status": "exit_status"}
    out = {}
    for line in lines:
        s = line.strip()
        for k, name in keys.items():
            if s.startswith(k + ":"):
                out[name] = s[len(k) + 1:].strip()
    w = [float(x) for x in out["wall"].split(":")]
    out["wall_s"] = w[0] * 3600 + w[1] * 60 + w[2] if len(w) == 3 else w[0] * 60 + w[1]
    out["user_s"] = float(out["user_s"])
    out["system_s"] = float(out["system_s"])
    out["cpu_percent"] = float(out["cpu_percent"].rstrip("%"))
    out["max_rss_kb"] = int(out["max_rss_kb"])
    out["max_rss_gib"] = out["max_rss_kb"] / 1024 ** 2
    return out


def cpp_timing(lines: list[str]) -> dict:
    for line in lines:
        m = re.search(r"Propagation time: ([\d.]+) s \(([\d:]+), (\d+) steps, ([\d.]+) s/step\)", line)
        if m:
            return {"propagation_s": float(m.group(1)), "steps": int(m.group(3)), "s_per_step": float(m.group(4))}
    return {}


def fortran_timing(lines: list[str]) -> dict:
    out = {}
    active = False
    for line in lines:
        if line.strip() == "CPU TIMES":
            active = True
            continue
        if active:
            m = re.match(r"^\s*(\S.*?)\s{2,}([\d.]+)(?:\s*\(\w\))?\s*$", line)
            if m:
                out.setdefault(m.group(1), float(m.group(2)))
            if line.strip().startswith("Command being timed"):
                break
    return out


def main():
    res = {"generated": dt.datetime.now().astimezone().isoformat(timespec="seconds")}
    # git: commit of the builds = last commit touching src/ before the binaries were linked
    res["git_head"] = sh("git", "rev-parse", "--short", "HEAD")
    res["git_head_date"] = sh("git", "log", "-1", "--format=%ad", "--date=iso-strict")
    res["git_last_src_commit"] = sh("git", "log", "-1", "--format=%h", "--", "src/")
    res["git_last_src_commit_date"] = sh("git", "log", "-1", "--format=%ad", "--date=iso-strict", "--", "src/")
    base = res["git_last_src_commit"]
    res["src_diff_since_build_commit"] = sh("git", "diff", "--stat", f"{base}..HEAD", "--", "src/", "tests/",
                                            "CMakeLists.txt", "cmake/")
    res["worktree_status_src"] = sh("git", "status", "--porcelain", "--", "src/", "tests/", "CMakeLists.txt", "cmake/")
    res["input_commit"] = sh("git", "log", "-1", "--format=%h %ad", "--date=iso-strict", "--", "inputs/thesis_fig5_3_dz/input.json")
    res["input_dirty"] = sh("git", "status", "--porcelain", "--", "inputs/thesis_fig5_3_dz/")
    builds = {}
    for k, b in BUILDS.items():
        newer = sh("find", "src", "tests", "CMakeLists.txt", "cmake", "-newer", str(b / "spodvr"), "-type", "f")
        builds[k] = {"dir": b.name, "cache": cmake_cache(b), "spodvr_mtime": mtime(b / "spodvr"),
                     "libspodvr_core_mtime": mtime(b / "libspodvr_core.a"),
                     "sources_newer_than_binary": newer.splitlines()}
    builds["fortran"] = {"exe": str(FORTRAN_EXE).replace(str(Path.home()), "~"), "mtime": mtime(FORTRAN_EXE),
                         "ldd": [l.split()[0] for l in sh("ldd", str(FORTRAN_EXE)).splitlines()],
                         "flags": "-O3 -fopenmp -fdefault-real-8 (GFORTRAN.mk), FFTW_ESTIMATE plans"}
    res["builds"] = builds
    runs = {}
    for k in ("cpu", "gpu", "fortran"):
        lines = clean_log_lines(k)
        r = script_header(RAW_LOGS[k])
        r["time"] = time_v(lines)
        r["binary_older_than_run"] = None
        runs[k] = r
    runs["cpu"]["program_timing"] = cpp_timing(clean_log_lines("cpu"))
    runs["gpu"]["program_timing"] = cpp_timing(clean_log_lines("gpu"))
    runs["fortran"]["program_timing"] = fortran_timing(clean_log_lines("fortran"))
    fl = clean_log_lines("fortran")
    runs["fortran"]["header"] = next(l for l in fl if "SPO-DVR program version" in l).strip("* ").strip()
    for k, b in (("cpu", "spodvr_mtime"), ("gpu", "spodvr_mtime")):
        runs[k]["binary_older_than_run"] = dt.datetime.fromisoformat(builds[k][b]) < dt.datetime.fromisoformat(runs[k]["start"])
    runs["fortran"]["binary_older_than_run"] = dt.datetime.fromisoformat(builds["fortran"]["mtime"]) < dt.datetime.fromisoformat(runs["fortran"]["start"])
    for k, t in (("cpu", builds["cpu"]["spodvr_mtime"]), ("gpu", builds["gpu"]["spodvr_mtime"]), ("fortran", builds["fortran"]["mtime"])):
        runs[k]["binary_age_h"] = (dt.datetime.fromisoformat(runs[k]["start"]) - dt.datetime.fromisoformat(t)).total_seconds() / 3600
    # libraries reported by the C++ program
    for k in ("cpu", "gpu"):
        lines = clean_log_lines(k)
        runs[k]["libraries"] = next(l.strip() for l in lines if l.strip().startswith("BLAS/LAPACK:"))
        runs[k]["banner"] = next(l.strip() for l in lines if l.strip().startswith("SPO-DVR (C++) version"))
        gpu = [l.strip() for l in lines if l.strip().startswith("GPU NVIDIA")]
        if gpu:
            runs[k]["device"] = gpu[0]
    res["runs"] = runs
    # run directories: final files present
    res["files"] = {k: sorted(p.name for p in d.iterdir()) for k, d in RUN_DIRS.items()}
    # input files: DVR.inp of the Fortran run vs the committed one
    res["dvr_inp_identical_to_committed"] = (RUN_DIRS["fortran"] / "DVR.inp").read_bytes() == (INPUT_DIR / "DVR.inp").read_bytes()
    # hardware
    lscpu = sh("lscpu")
    def lsc(key):
        m = re.search(rf"^{re.escape(key)}:\s*(.*)$", lscpu, re.M)
        return m.group(1).strip() if m else ""
    res["hardware"] = {
        "cpu_model": lsc("Model name"), "logical_cpus": lsc("CPU(s)"), "threads_per_core": lsc("Thread(s) per core"),
        "cores_per_socket": lsc("Core(s) per socket"), "numa_nodes": lsc("NUMA node(s)"), "l3": lsc("L3 cache"),
        "mem_total_gib": round(int(re.search(r"MemTotal:\s+(\d+)", Path("/proc/meminfo").read_text()).group(1)) / 1024 ** 2, 1),
        "gpu": sh("nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"),
        "kernel": sh("uname", "-r"),
        "os": re.search(r'PRETTY_NAME="([^"]+)"', Path("/etc/os-release").read_text()).group(1),
        "gxx": sh("g++", "--version").splitlines()[0], "gfortran": sh("gfortran", "--version").splitlines()[0],
        "nvcc": sh("/usr/local/cuda/bin/nvcc", "--version").splitlines()[-2],
        "fftw_package": sh("dpkg-query", "-W", "-f=${Version}", "libfftw3-double3"),
    }
    save_item("phase0", res)
    print(json.dumps({k: res[k] for k in ("git_head", "git_last_src_commit", "src_diff_since_build_commit",
                                         "worktree_status_src", "dvr_inp_identical_to_committed")}, indent=1))
    for k, r in runs.items():
        print(k, r["start"], "->", r["end"], "exit", r["exit_code"], "binary older:", r["binary_older_than_run"],
              "wall", r["time"]["wall"], "rss GiB %.2f" % r["time"]["max_rss_gib"])
    for k in ("cpu", "gpu"):
        print(k, builds[k]["cache"].get("SPODVR_FAST_MATH"), builds[k]["cache"].get("SPODVR_ARCH"),
              "newer sources:", builds[k]["sources_newer_than_binary"])


if __name__ == "__main__":
    main()

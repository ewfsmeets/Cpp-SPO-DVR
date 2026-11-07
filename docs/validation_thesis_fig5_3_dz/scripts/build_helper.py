#!/usr/bin/env python3
"""Compile and link a validation helper against the libraries of the build that ran.

Usage: build_helper.py <build dir name> <helper.cpp>   (e.g. build-fast potential_check.cpp)

Compile flags: those of src/app/run.cpp in <build>/compile_commands.json (same defines,
includes, -O3 -ffast-math -march=native, LTO). Link: the `spodvr` link line of
<build>/build.ninja with main.cpp.o replaced by the helper object. The libraries hold slim LTO
objects, so the machine code is generated at this link: same source and flags as the run, not
the same machine code (hence the bitwise self-checks). Output (binary and object) goes to
runs/thesis_fig5_3_dz_validation/helpers/<build>/, never into the build directory."""
from __future__ import annotations

import json
import re
import shlex
import subprocess
import sys

from common import REPO, SCRIPTS, WORK


def compile_flags(build) -> list[str]:
    cc = json.loads((build / "compile_commands.json").read_text())
    entry = next(e for e in cc if e["file"].endswith("src/app/run.cpp"))
    args = shlex.split(entry["command"])
    out = []
    skip = False
    for a in args[1:]:
        if skip:
            skip = False
            continue
        if a in ("-o", "-c"):
            skip = True
            continue
        out.append(a)
    # the core library's include of oneMKL's FFTW headers (fft/fftw_fft.hpp is included by helpers)
    core = next(e for e in cc if e["file"].endswith("src/fft/fftw_fft.cpp"))
    for a in shlex.split(core["command"]):
        if a.startswith("-I") and a not in out:
            out.append(a)
    return args[0], out


def link_line(build) -> tuple[str, list[str]]:
    text = (build / "build.ninja").read_text()
    block = re.search(r"^build spodvr: .*?\n((?:  .*\n)+)", text, re.M).group(1)
    flags = re.search(r"^\s+FLAGS = (.*)$", block, re.M).group(1)
    libs = re.search(r"^\s+LINK_LIBRARIES = (.*)$", block, re.M).group(1)
    path = re.search(r"^\s+LINK_PATH = (.*)$", block, re.M)
    return shlex.split(flags), (shlex.split(path.group(1)) if path else []) + shlex.split(libs)


def main():
    build = REPO / sys.argv[1]
    source = SCRIPTS / sys.argv[2]
    name = source.stem
    out_dir = WORK / "helpers" / build.name
    out_dir.mkdir(parents=True, exist_ok=True)
    compiler, cflags = compile_flags(build)
    obj = out_dir / f"{name}.o"
    exe = out_dir / name
    cmd = [compiler, *cflags, "-o", str(obj), "-c", str(source)]
    subprocess.run(cmd, check=True)
    flags, libs = link_line(build)
    libs = [str(build / l) if l.endswith(".a") and not l.startswith("/") else l for l in libs]
    cmd = [compiler, *flags, str(obj), "-o", str(exe), *libs]
    subprocess.run(cmd, check=True, cwd=out_dir)
    (out_dir / f"{name}.build.txt").write_text(
        "compile: " + " ".join([compiler, *cflags]) + "\nlink: " + " ".join([compiler, *flags, "<obj>", *libs]) + "\n")
    print(f"built {exe.relative_to(REPO)}")


if __name__ == "__main__":
    main()

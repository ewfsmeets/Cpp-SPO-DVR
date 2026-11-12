#!/usr/bin/env python3
"""Compile and link a C++ helper of the flux replay against the libraries of a build.

Compile flags: those of src/app/run.cpp in BUILD_DIR/compile_commands.json (defines, include
directories incl. nlohmann_json, optimisation and floating-point flags). Link: the `spodvr` link
line of BUILD_DIR/build.ninja with main.cpp.o replaced by the helper's object. The flux analysis
itself (src/analysis/flux_analysis.cpp) is strict floating point in every build; use the default
build (`build`) for the bitwise comparisons.

Usage: python3 -I -B build_helper.py BUILD_DIR HELPER.cpp OUT_DIR
  e.g. build_helper.py build docs/validation_flux_replay/replay_flux.cpp /tmp/flux_replay/helpers
Output: OUT_DIR/<helper name> (and its object), never inside BUILD_DIR; the last stdout line names it.
Environment: none read. Run time: about 30 s (the link runs the link-time optimisation).
Docs: docs/validation_flux_replay/README.md.
"""
import json
import pathlib
import re
import shlex
import subprocess
import sys


def compile_command(build: pathlib.Path) -> tuple[str, list[str]]:
    """The compiler and the flags of src/app/run.cpp, without -o/-c and their arguments."""
    entries = json.loads((build / "compile_commands.json").read_text())
    entry = next((e for e in entries if e["file"].endswith("src/app/run.cpp")), None)
    if entry is None:
        raise LookupError(f"{build}/compile_commands.json has no src/app/run.cpp")
    args = shlex.split(entry["command"])
    flags = []
    skip = False
    for arg in args[1:]:
        if skip:
            skip = False
        elif arg in ("-o", "-c"):
            skip = True
        else:
            flags.append(arg)
    return args[0], flags


def link_command(build: pathlib.Path) -> tuple[list[str], list[str]]:
    """The link flags and libraries of the spodvr executable (library paths made absolute)."""
    text = (build / "build.ninja").read_text()
    block = re.search(r"^build spodvr: .*?\n((?:  .*\n)+)", text, re.M)
    if block is None:
        raise LookupError(f"{build}/build.ninja has no spodvr link step")
    variables = dict(re.findall(r"^\s+(\w+) = (.*)$", block.group(1), re.M))
    flags = shlex.split(variables.get("FLAGS", "")) + [
        f for f in shlex.split(variables.get("LINK_FLAGS", "")) if not f.startswith("-Wl,--dependency-file")
    ]
    libraries = shlex.split(variables.get("LINK_PATH", "")) + shlex.split(variables.get("LINK_LIBRARIES", ""))
    libraries = [str(build / lib) if lib.endswith(".a") and not lib.startswith("/") else lib for lib in libraries]
    return flags, libraries


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: build_helper.py BUILD_DIR HELPER.cpp OUT_DIR", file=sys.stderr)
        return 2
    build = pathlib.Path(sys.argv[1]).resolve()
    source = pathlib.Path(sys.argv[2]).resolve()
    out_dir = pathlib.Path(sys.argv[3]).resolve()
    if out_dir == build or build in out_dir.parents:
        print(f"error: OUT_DIR {out_dir} is inside the build directory; name a scratch directory", file=sys.stderr)
        return 1
    try:
        compiler, compile_flags = compile_command(build)
        link_flags, libraries = link_command(build)
    except (OSError, LookupError) as error:
        print(f"error: {error}; name a configured and built build directory", file=sys.stderr)
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)
    obj = out_dir / f"{source.stem}.o"
    exe = out_dir / source.stem
    try:
        subprocess.run([compiler, *compile_flags, "-o", str(obj), "-c", str(source)], check=True)
        subprocess.run([compiler, *link_flags, str(obj), "-o", str(exe), *libraries], check=True, cwd=out_dir)
    except subprocess.CalledProcessError as error:
        print(f"error: {error.cmd[0]} failed (exit {error.returncode}); see its messages above", file=sys.stderr)
        return 1
    print(f"built {exe}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

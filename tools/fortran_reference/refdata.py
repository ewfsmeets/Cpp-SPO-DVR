#!/usr/bin/env python3
"""Provenance of the Fortran reference data: the MANIFEST.tsv of each reference case.

Every file a Fortran run generates for the regression tests has one line in the MANIFEST.tsv of
its case directory (tests/reference/<case>/, tests/crplib/reference/<case>/): its path, sha256 and
how it was made. The columns are described in tools/fortran_reference/README.md; ctest
ReferenceManifest.FilesMatch checks the hashes.

  manifest record CASE_DIR FILE... --script NAME (--install DIR --binary FILE | --pre-manifest)
                  [--input PATH]... [--runtime-drift TEXT]
      Writes the lines of FILE... (paths relative to CASE_DIR, just generated) into
      CASE_DIR/MANIFEST.tsv, replacing their old lines; creates the MANIFEST if needed. The run
      scripts call it through record_manifest (reference_config.sh). An --input directory (a CRP
      data set) is hashed as a tree (tree_sha256).
  manifest delete CASE_DIR...
      Deletes the files each CASE_DIR/MANIFEST.tsv lists, the directories this leaves empty, and
      the MANIFEST.tsv (regenerate_all.sh, before it regenerates). Every MANIFEST is read and
      checked before anything is deleted.
Both refuse the committed inputs of a case (is_input) and a directory that holds no Fortran
output (NOT_A_CASE), so that a wrong call can never make a regeneration delete an input.

Usage: python3 -I -B tools/fortran_reference/refdata.py manifest {record,delete} ... (above; -h)
Exit status: 0 done, 1 refused ("error: <what>; <fix>"), 2 bad arguments.
Environment: none read. Run time: under a second. Docs: tools/fortran_reference/README.md.
"""

import argparse
import datetime
import hashlib
import os
import subprocess
import sys
from pathlib import Path, PurePosixPath

REPO_DIR = Path(__file__).resolve().parents[2]
MANIFEST_NAME = "MANIFEST.tsv"
MANIFEST_COLUMNS = ("path", "sha256", "script", "commit", "fortran_ref", "input_sha", "runtime_drift",
                    "cpu", "date")
MANIFEST_WRITER_LINES = (
    "# Provenance of the Fortran reference files of this case, written by",
    "# tools/fortran_reference/refdata.py manifest. Columns: tools/fortran_reference/README.md.",
)
BUILDINFO_COLUMNS = ("key", "value")
EMPTY = "-"
SHORT_HASH = 12  # hex digits of the secondary hashes (programs, inputs, CPU flags)
PRE_MANIFEST_REF = "5.19 (pre-MANIFEST)"
# Everything but the reference trees: their tracked files change in every regeneration, so they
# do not make the commit "dirty".
OUTSIDE_REFERENCE_TREES = ":(exclude,glob)tests/**/reference/**"
# The committed inputs of a case (see is_input) and the directories of the reference trees that
# hold no Fortran output (fortran_unformatted/: hand-made gfortran test records).
CASE_INPUT_FILES = ("DVR.inp", "input.json")
CASE_INPUT_SUFFIX = ".in"  # the CRP driver's point lists
CASE_INPUT_DIR = "dataset"  # a CRP data set
NOT_A_CASE = ("fortran_unformatted",)


class Refusal(Exception):
    """A refusal, printed as "error: <what>; <fix>" with exit status 1."""


def sha256_of(path):
    """The sha256 (hex) of the file at path."""
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_sha256(directory):
    """The sha256 of the sorted lines "<sha256>  <path>" (sha256sum format, paths relative to
    directory) of every file below it: a data set hashes the same wherever it is staged. Shell:
    (cd DIR && find . -type f -printf '%P\\n' | LC_ALL=C sort | xargs sha256sum) | sha256sum"""
    files = sorted(str(f.relative_to(directory)) for f in Path(directory).rglob("*") if f.is_file())
    listing = "".join(f"{sha256_of(Path(directory) / f)}  {f}\n" for f in files)
    return hashlib.sha256(listing.encode()).hexdigest()


def read_tsv(path, columns):
    """The rows of a TSV file ('#' writer lines, then the header `columns`) as dicts."""
    rows, header_seen = [], False
    for number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if not header_seen:
            if tuple(fields) != columns:
                raise Refusal(f"{path}:{number}: the header is not {'<TAB>'.join(columns)}; "
                              "restore the file from git")
            header_seen = True
        elif len(fields) != len(columns):
            raise Refusal(f"{path}:{number}: {len(fields)} fields, expected {len(columns)}; "
                          "restore the file from git")
        else:
            rows.append(dict(zip(columns, fields)))
    if not header_seen:
        raise Refusal(f"{path} has no header row; restore the file from git")
    return rows


def write_manifest(path, rows):
    """Writes MANIFEST rows sorted by path, via a temporary file and a rename."""
    lines = list(MANIFEST_WRITER_LINES) + ["\t".join(MANIFEST_COLUMNS)]
    for row in sorted(rows, key=lambda r: r["path"]):
        lines.append("\t".join(row[column] for column in MANIFEST_COLUMNS))
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def checked_field(name, value):
    """value as a MANIFEST field: never empty, no tab or line break (TSV), and no ';', '[' or ']'
    (ReferenceManifest.cmake reads the lines as a CMake list: ';' splits a line, an unbalanced
    bracket joins it with the next)."""
    if value == "" or any(c in value for c in "\t\n\r;[]"):
        raise Refusal(f"{name} {value!r} is empty or holds a tab, line break, ';', '[' or ']'; change it")
    return value


def is_input(path):
    """True for a committed input of a case (never generated, never deleted): the case's DVR.inp
    and input.json, the CRP point lists (*.in) and data set (dataset/). Everything else in a case
    directory is generated, including fine_<variant>/DVR.inp."""
    parts = PurePosixPath(path).parts
    at_top = len(parts) == 1
    return (parts[0] == CASE_INPUT_DIR or (at_top and parts[0] in CASE_INPUT_FILES)
            or (at_top and parts[0].endswith(CASE_INPUT_SUFFIX)))


def check_case_dir(case_dir):
    """Refuses a directory of the reference trees that is not a case of Fortran outputs."""
    if case_dir.name in NOT_A_CASE:
        raise Refusal(f"{case_dir} holds no Fortran output; name a reference case directory")


def case_relative(case_dir, name):
    """The case-relative path `name` as a string, refused if it leaves the case directory or names
    the MANIFEST or a committed input."""
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or path.name == MANIFEST_NAME or str(path) == ".":
        raise Refusal(f"{name} is not a file path inside {case_dir}; name it relative to the case directory")
    if is_input(path):
        raise Refusal(f"{name} is a committed input of {case_dir}, not a generated file; "
                      "record only the files a run wrote")
    return str(path)


def git(*arguments):
    """The stripped stdout of a git command in this repository."""
    result = subprocess.run(["git", "-C", str(REPO_DIR), *arguments], capture_output=True, text=True)
    if result.returncode != 0:
        raise Refusal(f"git {' '.join(arguments)} failed in {REPO_DIR}: {result.stderr.strip()}; "
                      "run from a git clone or worktree")
    return result.stdout.strip()


def repository_commit():
    """`git describe --always` of this repository, plus "-dirty" when a tracked file outside the
    reference trees differs from it (regenerating changes those trees, so they cannot count)."""
    changed = git("--no-optional-locks", "status", "--porcelain", "--untracked-files=no", "--", ".",
                  OUTSIDE_REFERENCE_TREES)
    return git("describe", "--always") + ("-dirty" if changed else "")


def fortran_ref(install, binary):
    """'<install id> <program>:<sha256, 12 digits>': the install (BUILDINFO.tsv `id`) and the
    program that ran (its executable, or the CRP driver linked against its objects)."""
    buildinfo = Path(install) / "BUILDINFO.tsv"
    if not buildinfo.is_file():
        raise Refusal(f"{buildinfo} not found; point --install at a reference install")
    keys = {row["key"]: row["value"] for row in read_tsv(buildinfo, BUILDINFO_COLUMNS)}
    if "id" not in keys:
        raise Refusal(f"{buildinfo} has no id row; restore the install's BUILDINFO.tsv")
    return f"{keys['id']} {Path(binary).name}:{sha256_of(binary)[:SHORT_HASH]}"


def input_hashes(inputs):
    """'<name>:<sha256, 12 digits>' of each input as run, comma-separated ('-': none). The name is
    the path as given if relative, else its last component; a directory hashes as a tree."""
    entries = []
    for given in inputs:
        path = Path(given)
        name = path.name if path.is_absolute() else given
        digest = tree_sha256(path) if path.is_dir() else sha256_of(path)
        entries.append(f"{name}:{digest[:SHORT_HASH]}")
    return ",".join(entries) or EMPTY


def cpu_description():
    """'<model name> flags:<sha256 of the flags, 12 digits>' of the first CPU (/proc/cpuinfo):
    libm and FFTW pick code paths by CPU feature, so a result can depend on the flags."""
    info = {}
    try:
        for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
            key, _, value = line.partition(":")
            info.setdefault(key.strip(), value.strip())
    except OSError:
        return EMPTY
    if "model name" not in info or "flags" not in info:
        return EMPTY
    flags = hashlib.sha256(info["flags"].encode()).hexdigest()[:SHORT_HASH]
    return f"{info['model name']} flags:{flags}"


def provenance(args):
    """The columns after path and sha256, the same for every file of one record call."""
    if args.pre_manifest:
        values = dict(commit=EMPTY, fortran_ref=PRE_MANIFEST_REF, input_sha=EMPTY,
                      runtime_drift=EMPTY, cpu=EMPTY, date=EMPTY)
    else:
        values = dict(commit=repository_commit(), fortran_ref=fortran_ref(args.install, args.binary),
                      input_sha=input_hashes(args.input), runtime_drift=args.runtime_drift,
                      cpu=cpu_description(),
                      date=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d"))
    values["script"] = args.script
    return {name: checked_field(name, value) for name, value in values.items()}


def record(args):
    """manifest record: replaces the lines of the named files in the case's MANIFEST.tsv."""
    case_dir = Path(args.case_dir).resolve()
    if not case_dir.is_dir():
        raise Refusal(f"case directory {args.case_dir} not found; name an existing one")
    check_case_dir(case_dir)
    common = provenance(args)
    new_rows = {}
    for name in args.files:
        path = case_relative(case_dir, name)
        if not (case_dir / path).is_file():
            raise Refusal(f"{case_dir / path} is not a file; record only files the run wrote")
        new_rows[path] = dict(path=checked_field("path", path), sha256=sha256_of(case_dir / path), **common)
    manifest = case_dir / MANIFEST_NAME
    old_rows = read_tsv(manifest, MANIFEST_COLUMNS) if manifest.exists() else []
    kept = [row for row in old_rows if row["path"] not in new_rows]
    write_manifest(manifest, kept + list(new_rows.values()))
    print(f"* {manifest}: {len(new_rows)} line(s) written, {len(kept)} kept")


def delete(args):
    """manifest delete: the listed files, the directories left empty, then the MANIFESTs."""
    plans = []
    for name in args.case_dirs:
        case_dir = Path(name).resolve()
        check_case_dir(case_dir)
        manifest = case_dir / MANIFEST_NAME
        if not manifest.is_file():
            raise Refusal(f"{manifest} not found; name case directories that have a {MANIFEST_NAME}")
        paths = [case_relative(case_dir, row["path"]) for row in read_tsv(manifest, MANIFEST_COLUMNS)]
        plans.append((case_dir, manifest, paths))
    for case_dir, manifest, paths in plans:
        missing = 0
        for path in paths:
            target = case_dir / path
            if not target.exists():
                missing += 1
                continue
            target.unlink()
            for parent in target.parents:
                if parent == case_dir or any(parent.iterdir()):
                    break
                parent.rmdir()
        manifest.unlink()
        note = f" ({missing} listed file(s) were already missing)" if missing else ""
        print(f"* {case_dir}: {len(paths) - missing} file(s) and {MANIFEST_NAME} deleted{note}")


def parse_arguments(argv):
    """The command line; argparse exits with status 2 and a usage line on bad arguments."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    manifest = commands.add_parser("manifest", help="write or delete MANIFEST.tsv lines")
    actions = manifest.add_subparsers(dest="action", required=True)

    add = actions.add_parser("record", help="replace the lines of just generated files")
    add.add_argument("case_dir", metavar="CASE_DIR", help="the case directory (holds MANIFEST.tsv)")
    add.add_argument("files", metavar="FILE", nargs="+", help="generated file, relative to CASE_DIR")
    add.add_argument("--script", metavar="NAME", required=True, help="the script that generated the files")
    add.add_argument("--install", metavar="DIR", help="the reference install that ran (its BUILDINFO.tsv gives the id)")
    add.add_argument("--binary", metavar="FILE", help="the program that ran: the install's executable or the CRP driver")
    add.add_argument("--input", metavar="PATH", action="append", default=[], help="an input as run (repeatable): a file, or a directory (a data set, hashed as a tree); "
                          "named by the path as given if relative, else by its last component")
    add.add_argument("--runtime-drift", metavar="TEXT", default=EMPTY,
                     help="run-time packages that differ from BUILDINFO.tsv (default: -)")
    add.add_argument("--pre-manifest", action="store_true",
                     help=f'files made before MANIFESTs existed: fortran_ref "{PRE_MANIFEST_REF}", "-" '
                          "for the columns not recorded then")

    remove = actions.add_parser("delete", help="delete the listed files and the MANIFESTs")
    remove.add_argument("case_dirs", metavar="CASE_DIR", nargs="+", help="a case directory with a MANIFEST.tsv")

    args = parser.parse_args(argv)
    if args.action == "record":
        run_options = args.install or args.binary or args.input or args.runtime_drift != EMPTY
        if args.pre_manifest and run_options:
            add.error("--pre-manifest takes none of --install, --binary, --input, --runtime-drift")
        if not args.pre_manifest and not (args.install and args.binary):
            add.error("--install and --binary are required (or --pre-manifest)")
    return args


def main(argv):
    """Runs one subcommand; a refusal exits with status 1."""
    args = parse_arguments(argv)
    try:
        if args.action == "record":
            record(args)
        else:
            delete(args)
    except Refusal as refusal:
        print(f"error: {refusal}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

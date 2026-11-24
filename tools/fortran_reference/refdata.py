#!/usr/bin/env python3
"""The Fortran reference data: the MANIFEST.tsv of each reference case, and the comparison of two
reference trees.

Every file a Fortran run generates for the regression tests has one line in the MANIFEST.tsv of
its case directory (tests/reference/<case>/, tests/crplib/reference/<case>/,
tests/scm/reference/<case>/): its path, sha256 and how it was made. The columns are described in
tools/fortran_reference/README.md; ctest ReferenceManifest.FilesMatch checks the hashes.

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

  compare OLD NEW --expected FILE [--no-provenance]
      Compares the reference files (tests/reference/**, tests/crplib/reference/**,
      tests/scm/reference/**) of the repository trees OLD (e.g. the main checkout) and NEW (e.g. a
      worktree after regenerate_all.sh). Each file is classified by the row of the
      expected-differences file FILE whose glob matches it (a file two rows match is refused); a
      file no row matches must be identical. Prints one line per file that is not identical, then a
      summary. NEW's provenance is checked too: every case has a MANIFEST.tsv that lists exactly its
      generated files with their sha256, all made by the install FILE names ("# install: <id>"),
      whose BUILDINFO copy (tools/fortran_reference/buildinfo/<id>.tsv) has the executable the lines
      name, without run-time drift, and every fine analysis read the case's fortran_fort21.bin.
      --no-provenance skips only these MANIFEST and BUILDINFO checks (e.g. for files copied in by
      hand).

Usage: python3 -I -B tools/fortran_reference/refdata.py {manifest {record,delete},compare} ... (-h)
Exit status: 0 done (compare: every difference allowed), 1 refused or a difference not allowed
("error: <what>; <fix>"), 2 bad arguments.
Environment: none read. Run time: manifest under a second; compare 4-7 s on the full trees
(measured 2026-11-12; it reads every BK coefficient of the files that differ).
Docs: tools/fortran_reference/README.md.
"""

import argparse
import cmath
import datetime
import functools
import hashlib
import math
import os
import re
import struct
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

REPO_DIR = Path(__file__).resolve().parents[2]
MANIFEST_NAME = "MANIFEST.tsv"
MANIFEST_COLUMNS = ("path", "sha256", "script", "commit", "fortran_ref", "input_sha", "runtime_drift", "cpu", "date")
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
CASE_INPUT_FILES = ("DVR.inp", "input.json", "pes_cuts.json")  # pes_cuts.json: the write2DPES cuts (gate 1c)
CASE_INPUT_SUFFIX = ".in"  # the CRP and SCM drivers' inputs (point lists, SCM parameter and slab files)
CASE_INPUT_DIR = "dataset"  # a CRP data set
NOT_A_CASE = ("fortran_unformatted",)
# The reference trees, one directory per case. Copies: reference_config.sh FORTRAN_REFERENCE_TREES,
# ReferenceManifest.cmake (also its MANIFEST header and NOT_A_CASE), tools/check_docs.py and the
# exclude of .pre-commit-config.yaml; ctest FortranReference.ReferenceListsAgree checks them.
REFERENCE_TREES = ("tests/reference", "tests/crplib/reference", "tests/scm/reference")


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
                raise Refusal(f"{path}:{number}: the header is not {'<TAB>'.join(columns)}; restore the file from git")
            header_seen = True
        elif len(fields) != len(columns):
            raise Refusal(f"{path}:{number}: {len(fields)} fields, expected {len(columns)}; restore the file from git")
        else:
            rows.append(dict(zip(columns, fields, strict=True)))
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
    """True for a committed input of a case (never generated, never deleted): the case's DVR.inp,
    input.json and pes_cuts.json, the driver inputs (*.in: CRP and SCM point lists, SCM parameter and
    slab files) and the CRP data set (dataset/). Everything else in a case directory is generated,
    including fine_<variant>/DVR.inp."""
    parts = PurePosixPath(path).parts
    at_top = len(parts) == 1
    return (
        parts[0] == CASE_INPUT_DIR
        or (at_top and parts[0] in CASE_INPUT_FILES)
        or (at_top and parts[0].endswith(CASE_INPUT_SUFFIX))
    )


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
        raise Refusal(
            f"{name} is a committed input of {case_dir}, not a generated file; record only the files a run wrote"
        )
    return str(path)


def git(*arguments):
    """The stripped stdout of a git command in this repository."""
    result = subprocess.run(["git", "-C", str(REPO_DIR), *arguments], capture_output=True, text=True)
    if result.returncode != 0:
        raise Refusal(
            f"git {' '.join(arguments)} failed in {REPO_DIR}: {result.stderr.strip()}; run from a git clone or worktree"
        )
    return result.stdout.strip()


def repository_commit():
    """`git describe --always` of this repository, plus "-dirty" when a tracked file outside the
    reference trees differs from it (regenerating changes those trees, so they cannot count)."""
    changed = git(
        "--no-optional-locks", "status", "--porcelain", "--untracked-files=no", "--", ".", OUTSIDE_REFERENCE_TREES
    )
    return git("describe", "--always") + ("-dirty" if changed else "")


def fortran_ref(install, binary):
    """'<install id> <program>:<sha256, 12 digits>': the install (BUILDINFO.tsv `id`) and the
    program that ran (its executable, or the CRP or SCM driver linked against its objects)."""
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
        values = dict(
            commit=EMPTY, fortran_ref=PRE_MANIFEST_REF, input_sha=EMPTY, runtime_drift=EMPTY, cpu=EMPTY, date=EMPTY
        )
    else:
        values = dict(
            commit=repository_commit(),
            fortran_ref=fortran_ref(args.install, args.binary),
            input_sha=input_hashes(args.input),
            runtime_drift=args.runtime_drift,
            cpu=cpu_description(),
            date=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d"),
        )
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


# ---------------------------------------------------------------------------------------------
# compare: a regenerated reference tree against the committed one
# ---------------------------------------------------------------------------------------------
BUILDINFO_COPIES = "tools/fortran_reference/buildinfo"  # the committed BUILDINFO.tsv of each install
FORTRAN_EXECUTABLE = "spodvr.x"  # the program name of the executable runs in fortran_ref
FINE_DIR_PREFIX = "fine_"  # the fine-analysis runs of a case: fine_<variant>/
FINE_INPUT = "fort.22"  # the BK file a fine analysis reads, as named in input_sha
EXPECTED_COLUMNS = ("glob", "allowed", "limit", "reason")
# The '#' line of the expected-differences file that names the install whose output it describes:
# NEW's MANIFEST lines must name that install.
EXPECTED_INSTALL = re.compile(r"# install: (\S+)")
# The closed list of classes of the expected-differences file (README "Comparing two reference
# trees"). The text classes nest: labels allows what volatile allows, and abs, printed_units and
# signs_only allow what labels allows.
CLASSES = (
    "identical",
    "volatile",
    "labels",
    "abs",
    "bk_scaled",
    "printed_units",
    "signs_only",
    "initial_phase",
    "unread",
    "removed",
    "added",
)
TEXT_CLASSES = ("volatile", "labels", "abs", "printed_units", "signs_only")
LIMIT_CLASSES = ("abs", "bk_scaled", "printed_units")  # their limit is a positive number

# The BK coefficients (fort.21): BK_HEADER_RECORDS records, then one record of complex coefficients
# per analysis step. Gate 4's scale is max(max|c| over all records, BK_SCALE_FLOOR)
# (tests/analysis_reference.hpp, compare_bk_with_fortran).
BK_FILE = "fortran_fort21.bin"
BK_HEADER_RECORDS = 10
BK_SCALE_FLOOR = 1e-4

# Lines that differ between two runs of one build: the version line (build and run dates), the
# thread count, the set-up timings, and the "iz = N" lines that the potential set-up prints from a
# parallel loop (their order depends on the threads). Dropped before every text comparison.
VOLATILE_LINE = re.compile(
    r"\*\*\* SPO-DVR program version .*|\* Using a maximum number of \d+ threads in this run\."
    r"| Total time (?:for potential initialisation|setting up the potential) = .* secs"
    r"| +iz = +\d+"
)
TIMING_TABLE = re.compile(r" CPU TIMES *")  # the timing table that ends a stdout: dropped from here on
# A CRP driver line "TAG index HEX decimal": the 17-digit decimal determines the double, so the hex
# word is dropped (its digits would read as numbers such as 4e5123).
HEX_DOUBLE_LINE = re.compile(r"^(\S+ +\d+ +)[0-9A-Fa-f]{16}( +\S+)$")

# The label whitelist of class labels: the stdout lines whose text changed from Fortran 5.19 to
# 5.24 while their numbers did not, as (pattern, canonical form); the canonical form None drops the
# line. Found by comparing every stdout of the r1713 -O3 build with the 5.19 references committed
# then (2026-11-11); expected/r1713.tsv names the files they apply to.
LABELS_519_TO_524 = (
    # The line is gone in 5.24.
    (re.compile(r"\* NAG calls disabled \.\.\."), None),
    # 5.19 prints the first point of each absorber (nzchsp, nschsp, nrchsp), 5.24 the point before
    # it (nzchmi, nschmi, nrchmi): other quantities, so both are dropped.
    (re.compile(r" +(?:[zr]\()?n[zsr]ch(?:sp|mi)\)? *=.*"), None),
    # "lattice constant in X/Y" became "in U/V".
    (re.compile(r"(\* The surface lattice constant in )U( \(bohr\):.*)"), r"\1X\2"),
    (re.compile(r"(\* The surface lattice constant in )V( \(bohr\):.*)"), r"\1Y\2"),
    # The headers of the absorber blocks gained "start" and "end".
    (re.compile(r"( +)(?:start|end) (optical potential in [zr] (?:scat|spec)\.grid:)"), r"\1\2"),
    # The absorber strengths gained "[high end]" and moved to a new column.
    (
        re.compile(
            r"(\* Optical potential in [Zr] )(?:\[high end\] )?((?:specular )?absorbs with coefficient:) +(\S+)"
        ),
        r"\1\2 \3",
    ),
)
# Labels changed after 5.24 (svn r1714-r1727, upstream fixes of 2026-11): the "(nproj = 2)" of the specular
# tolerance line became "(nproj = 1)" in r1725, with the same number. Used with LABELS_519_TO_524.
LABELS_R1725 = ((re.compile(r"(\* Tolerance for transferring the wavefunction \(nproj = )1(\):.*)"), r"\g<1>2\2"),)
# A number as Fortran prints it (I, F, E or D format), and a right-aligned field of one: a sign
# change takes a blank of the field, so only the field width must stay.
NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")
NUMBER_FIELD = re.compile(r" *" + NUMBER.pattern)
# The rounding a sign change may bring along, in units of the last printed digit: 1 of the 2277
# phases that change sign in cu111_dailight_small fine_full/fort.25 (5.19 -> r1713) also changes
# its last digit (-0.666745 -> 0.666744, measured 2026-11-12).
SIGN_CHANGE_ROUNDING = 1
# Relative: |a - b| / unit is computed in floating point ((0.666745 - 0.666744) / 1e-6 = 1.00000000003).
UNIT_SLACK = 1e-9


class FormatError(Exception):
    """A file that is not in the format its name implies (a gfortran record file)."""


@dataclass(frozen=True)
class Rule:
    """One row of the expected-differences file. limit: the number of abs, bk_scaled and
    printed_units, else 0; positions: the numbers of a line (1-based) whose sign signs_only lets
    change."""

    glob: str
    allowed: str
    limit: float
    positions: tuple
    reason: str
    pattern: re.Pattern


UNMATCHED = Rule("(no row)", "identical", 0.0, (), "a file no row matches must be identical", re.compile(""))


@dataclass
class Differences:
    """What differs between two versions of a file, and the first difference its rule does not
    allow ("" if none). For BK files, lines counts records and scale is the gate-4 scale."""

    lines: int = 0
    numbers: int = 0
    max_abs: float = 0.0
    max_units: float = 0.0  # in units of the last printed digit
    scale: float = 0.0
    line_counts: str = ""  # "N lines became M" when the number of lines changed
    violation: str = ""
    kind: str = "text"

    def forbid(self, what):
        """Records the first difference the rule does not allow."""
        self.violation = self.violation or what

    def summary(self):
        """One line: what differs, then what is not allowed."""
        if self.kind == "bk":
            text = (
                f"{self.lines} BK records differ, {self.numbers} coefficients; max |Δc| {self.max_abs:.3g} = "
                f"{self.max_abs / self.scale:.3g} of the gate-4 scale {self.scale:.3g}"
            )
        elif self.kind != "text":
            text = "the files differ"
        elif self.lines == 0 and not self.line_counts:
            text = "differs in volatile lines and whitelisted labels only"
        else:
            text = (
                f"{self.lines} lines differ, {self.numbers} numbers; max abs {self.max_abs:.3g}, "
                f"max {self.max_units:.3g} printed units"
            )
            text += f"; {self.line_counts}" if self.line_counts else ""
        return text + (f"; not allowed: {self.violation}" if self.violation else "")


def glob_pattern(glob):
    """The regular expression of a glob whose * and ? match within one path component."""
    parts = ("[^/]*" if c == "*" else "[^/]" if c == "?" else re.escape(c) for c in glob)
    return re.compile("".join(parts))


def parse_limit(allowed, limit):
    """(limit, positions) of a row of class allowed, or None if the limit does not fit the class."""
    try:
        if allowed in LIMIT_CLASSES:
            value = float(limit)
            return (value, ()) if math.isfinite(value) and value > 0 else None
        if allowed == "signs_only":
            positions = tuple(int(p) for p in limit.split(","))
            return (0.0, positions) if min(positions) >= 1 else None
    except ValueError:
        return None
    return (0.0, ()) if limit == EMPTY else None


def parse_rule(where, row):
    """The Rule of one row of the expected-differences file, refused if malformed."""
    allowed, limit = row["allowed"], row["limit"]
    if allowed not in CLASSES:
        raise Refusal(f"{where}: class {allowed!r} is not one of {', '.join(CLASSES)}; use one of them")
    if EMPTY in (row["glob"], row["reason"]):
        raise Refusal(f"{where}: a row needs a glob and a reason; add them")
    parsed = parse_limit(allowed, limit)
    if parsed is None:
        raise Refusal(
            f"{where}: limit {limit!r} does not fit class {allowed} (a positive number for "
            f"{', '.join(LIMIT_CLASSES)}; positions like 4,6 for signs_only; else -); fix it"
        )
    return Rule(row["glob"], allowed, *parsed, row["reason"], glob_pattern(row["glob"]))


def read_expected(path):
    """(the install id of its '# install:' line, the rules) of the expected-differences file."""
    rows = read_tsv(path, EXPECTED_COLUMNS)
    globs = [row["glob"] for row in rows]
    if globs != sorted(set(globs)):
        raise Refusal(f"{path}: the rows are not sorted by glob, or a glob repeats; sort them (LC_ALL=C sort)")
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    installs = [match.group(1) for match in map(EXPECTED_INSTALL.fullmatch, lines) if match]
    if len(installs) != 1:
        raise Refusal(
            f"{path} has {len(installs)} '# install: <id>' lines, not one; name the install whose "
            "output the file describes"
        )
    return installs[0], [parse_rule(f"{path}, row {row['glob']}", row) for row in rows]


def matching_rule(rules, path):
    """The rule whose glob matches path, and its index (UNMATCHED and -1 if none). Two matching rows
    are refused: the order of the rows would decide silently."""
    matches = [(index, rule) for index, rule in enumerate(rules) if rule.pattern.fullmatch(path)]
    if len(matches) > 1:
        raise Refusal(
            f"{path} matches the rows {' and '.join(rule.glob for _, rule in matches)}; "
            "make the globs of the expected file disjoint"
        )
    return matches[0] if matches else (-1, UNMATCHED)


def reference_files(root):
    """The paths (relative to root, POSIX) of every file of the reference trees but the MANIFESTs,
    which the provenance checks read."""
    found = set()
    for tree in REFERENCE_TREES:
        for path in (Path(root) / tree).rglob("*"):
            if path.is_file() and path.name != MANIFEST_NAME:
                found.add(path.relative_to(root).as_posix())
    return found


def comparable_lines(text, labels):
    """(line number, line) of a Fortran text output without its volatile lines and timing table,
    and without the hex words of the CRP driver; with labels, each whitelisted label line in its
    canonical form (or dropped)."""
    kept = []
    for number, line in enumerate(text.splitlines(), start=1):
        if TIMING_TABLE.fullmatch(line):
            break
        if VOLATILE_LINE.fullmatch(line):
            continue
        line = HEX_DOUBLE_LINE.sub(r"\1\2", line)
        for pattern, canonical in LABELS_519_TO_524 + LABELS_R1725 if labels else ():
            match = pattern.fullmatch(line)
            if match:
                line = None if canonical is None else match.expand(canonical)
                break
        if line is not None:
            kept.append((number, line))
    return kept


def e_format(token):
    """A printed number with an E exponent: D becomes E, and a three-digit exponent without its letter (gfortran
    prints 0.100000-119) gets it back."""
    text = token.upper().replace("D", "E")
    sign = max(text.rfind("+"), text.rfind("-"))
    if "E" not in text and sign > 0:
        text = text[:sign] + "E" + text[sign:]
    return text


def as_float(token):
    """The value of a printed number (D exponents and three-digit exponents without their letter included)."""
    return float(e_format(token))


def printed_unit(token):
    """The value of one unit in the last printed digit of a number."""
    mantissa, _, exponent = e_format(token).partition("E")
    decimals = len(mantissa.partition(".")[2])
    # Not 10.0 ** n, which overflows: hex digits in the CRP outputs read as numbers like 4e5123.
    return float(f"1e{int(exponent or 0) - decimals}")


def printed_digits(token):
    """(the digits of the mantissa as an integer, the exponent) of a number, without its sign:
    two numbers with the same exponent differ by that many units of the last printed digit."""
    mantissa, _, exponent = token.lstrip("+-").upper().replace("D", "E").partition("E")
    return int(mantissa.replace(".", "") or "0"), exponent


def sign_change_allowed(rule, position, old, new, row_length):
    """signs_only: in a table row of at least max(positions) numbers (an S-matrix row of fort.25,
    not a block header with its labels), the number at a listed position changed its sign, and its
    printed digits by at most SIGN_CHANGE_ROUNDING units (the negated value is computed along
    another path)."""
    (digits_old, exponent_old), (digits_new, exponent_new) = printed_digits(old), printed_digits(new)
    return (
        row_length >= max(rule.positions)
        and position in rule.positions
        and old.startswith("-") != new.startswith("-")
        and exponent_old == exponent_new
        and abs(digits_old - digits_new) <= SIGN_CHANGE_ROUNDING
    )


def number_allowed(rule, position, old, new, row_length):
    """Whether the rule lets the finite number at position (1-based) of a line change from old to
    new; row_length: how many numbers the line holds if it holds nothing else, else 0. The
    comparisons are written so that a NaN fails them."""
    difference = abs(as_float(old) - as_float(new))
    if rule.allowed == "abs":
        return difference <= rule.limit
    if rule.allowed == "printed_units":
        return difference / max(printed_unit(old), printed_unit(new)) <= rule.limit * (1 + UNIT_SLACK)
    if rule.allowed == "signs_only":
        return sign_change_allowed(rule, position, old, new, row_length)
    return False


def layout(line):
    """The line with each number field replaced by as many '#' as it is wide."""
    return NUMBER_FIELD.sub(lambda match: "#" * len(match.group()), line)


def compare_line(diff, rule, where, old, new):
    """Adds the differences of two versions of one line that differ."""
    diff.lines += 1
    old_numbers, new_numbers = NUMBER.findall(old), NUMBER.findall(new)
    # The same layout can hold fewer, wider numbers ("   0  0.94E-25" and "   0.00000000000").
    if layout(old) != layout(new) or len(old_numbers) != len(new_numbers):
        diff.forbid(f"{where}: text {old.strip()!r} became {new.strip()!r}")
        return
    row_length = 0 if NUMBER.sub("", old).strip() else len(old_numbers)
    for position, (a, b) in enumerate(zip(old_numbers, new_numbers, strict=True), start=1):
        if a == b:
            continue
        diff.numbers += 1
        x, y = as_float(a), as_float(b)
        if not (math.isfinite(x) and math.isfinite(y)):
            diff.forbid(f"{where}: number {position} {a} became {b}, not finite")
            continue
        difference = abs(x - y)
        diff.max_abs = max(diff.max_abs, difference)
        diff.max_units = max(diff.max_units, difference / max(printed_unit(a), printed_unit(b)))
        if not number_allowed(rule, position, a, b, row_length):
            diff.forbid(f"{where}: number {position} {a} became {b}")


def text_differences(old, new, rule):
    """The Differences of two Fortran text outputs, line by line after comparable_lines (the label
    whitelist for every class but volatile)."""
    labels = rule.allowed != "volatile"
    a = comparable_lines(old.decode(), labels)
    b = comparable_lines(new.decode(), labels)
    diff = Differences()
    for (number, x), (_, y) in zip(a, b, strict=False):  # the common lines; the count is checked below
        if x != y:
            compare_line(diff, rule, f"line {number}", x, y)
    if len(a) != len(b):
        diff.line_counts = f"{len(a)} lines became {len(b)} (without volatile and label lines)"
        diff.forbid(diff.line_counts)
    return diff


# Class initial_phase (an on-the-fly fortran_fort25.txt): upstream r1720 refers the phase of the initial
# rot-dif state's S elements (every v) to the fixed specular analysis line Z_sp instead of ZAnalysis
# (Analyzor.f90 smatrix_calc: S = ... exp(i arg), arg = -skz zRef + ...), so those rows rotate by
# exp(-i k_v (Z_sp - Zana)) and nothing else changes. k_v is the channel's normal momentum: the block's
# "total momentum" for v = VIn (EParallel 0), from the printed gas-phase levels for the other v.
CHAEV = 27.2114  # the Fortran's hartree -> eV (chaev)
CAMAU = 1822.89  # the Fortran's amu -> electron masses (camau)
GAS_LEVEL_DIGITS = 1e-6  # the gas-phase table prints its energies (eV) with 6 decimals
FORT25_HEADER = re.compile(r" v j mj dif\.order are: +(-?\d+) +(-?\d+) +(-?\d+) +(-?\d+) *")
FORT25_ROW = re.compile(r" +(-?\d+) +(-?\d+) +(\S+) +(\S+) +(\S+) +(\S+) *")
FORT25_MOMENTUM = re.compile(r" total momentum = +(\S+) *")
GAS_LEVELS_TITLE = " asymptotic molecular vibrational wave functions"
GAS_LEVEL = re.compile(r" +(\d+) +(\d+) +(-?\d+\.\d+) +(-?\d+\.\d+) .*")


@dataclass(frozen=True, kw_only=True)
class InitialPhaseCase:
    """What class initial_phase reads from a case: the initial state (VIn, JIn, MJIn, DVR.inp), the analysis
    line ZAnalysis and the fixed specular line Z_sp (printed, with half a printed unit each), the total mass
    (electron masses) and the gas-phase levels E(v, j) (eV) of fortran_stdout.txt."""

    v_in: int
    j_in: int
    mj_in: int
    z_ana: float
    z_sp: float
    z_rounding: float
    mass: float
    levels: dict


def dvr_value(text, key):
    """The value of KEY in a DVR.inp ("Key : TYPE : value;" or "Key : value;"), or None."""
    match = re.search(rf"^{key} *:(?:[^:;]*:)? *([^;]+);", text, re.MULTILINE)
    return match.group(1).strip() if match else None


def printed_value(text, label):
    """(value, half a printed unit) of the one stdout line "<label> <number>"; Refusal if not exactly one."""
    values = re.findall(rf"^{re.escape(label)} +(\S+) *$", text, re.MULTILINE)
    if len(values) != 1:
        raise Refusal(f"{len(values)} lines {label!r} in the case's fortran_stdout.txt, not one")
    return as_float(values[0]), 0.5 * printed_unit(values[0])


def gas_levels(text):
    """{(v, j): E (eV)} of the gas-phase table of a stdout (the lines after GAS_LEVELS_TITLE)."""
    lines = text.splitlines()
    if GAS_LEVELS_TITLE not in lines:
        raise Refusal("no gas-phase table in the case's fortran_stdout.txt")
    levels = {}
    for line in lines[lines.index(GAS_LEVELS_TITLE) + 1 :]:
        match = GAS_LEVEL.fullmatch(line)
        if match:
            levels[(int(match.group(1)), int(match.group(2)))] = float(match.group(3))
        elif line.strip() and levels:
            break
    return levels


def initial_phase_case(root, path):
    """The InitialPhaseCase of the case of path in the tree root (its DVR.inp and fortran_stdout.txt)."""
    located = case_and_relative_path(path)
    if located is None:
        raise Refusal("not in a case directory")
    case_dir = Path(root) / located[0]
    dvr = (case_dir / "DVR.inp").read_text(encoding="utf-8")
    stdout = (case_dir / "fortran_stdout.txt").read_text(encoding="utf-8")
    if as_float(dvr_value(dvr, "EParallel") or "0") != 0.0:
        raise Refusal("the case has EParallel /= 0 (the block's total momentum is not k_z then)")
    z_ana, z_ana_rounding = printed_value(stdout, "* Wave function analysed at (bohr):")
    z_sp, z_sp_rounding = printed_value(stdout, "   ==> ZAnalysis on specular channel fixed at Z =")
    mass_1, _ = printed_value(stdout, "* Mass of atom 1 (amu):")
    mass_2, _ = printed_value(stdout, "* Mass of atom 2 (amu):")
    return InitialPhaseCase(
        v_in=int(dvr_value(dvr, "VIn") or 0),
        j_in=int(dvr_value(dvr, "JIn") or 0),
        mj_in=int(dvr_value(dvr, "MJIn") or 0),
        z_ana=z_ana,
        z_sp=z_sp,
        z_rounding=z_ana_rounding + z_sp_rounding,
        mass=CAMAU * (mass_1 + mass_2),
        levels=gas_levels(stdout),
    )


def channel_momentum(case, v, momentum):
    """(k_v, its uncertainty from the printed digits) of the initial state's channel at level v, from the
    block's printed total momentum (k of v = VIn); None if a level is missing or the channel is closed."""
    k_in, dk_in = as_float(momentum), 0.5 * printed_unit(momentum)
    if v == case.v_in:
        return k_in, dk_in
    e_in, e_v = case.levels.get((case.v_in, case.j_in)), case.levels.get((v, case.j_in))
    if e_in is None or e_v is None:
        return None
    k_squared = k_in * k_in + 2.0 * case.mass * (e_in - e_v) / CHAEV
    if not k_squared > 0.0:
        return None
    k = math.sqrt(k_squared)
    dk_squared = 2.0 * k_in * dk_in + 2.0 * case.mass * GAS_LEVEL_DIGITS / CHAEV
    return k, dk_squared / (2.0 * k)


def initial_row_change(diff, case, where, v, momentum, old, new):
    """Records whether one initial-state row (Re S, Im S, P, phase as printed; changed or not) changed only as
    the r1720 rotation does: S_new = S_old exp(-i k_v (Z_sp - Zana)) and |S| within the printed digits of Re S
    and Im S (and the rounding of k_v, Zana and Z_sp), P within one printed unit, and the printed phase rotated
    by the same angle mod 2 pi. The comparisons are written so that a NaN fails them."""
    diff.numbers += sum(a != b for a, b in zip(old, new, strict=True))
    values = [as_float(token) for token in (*old, *new)]
    if not all(math.isfinite(x) for x in values):
        diff.forbid(f"{where}: a number is not finite")
        return
    re_old, im_old, p_old, phase_old, re_new, im_new, p_new, phase_new = values
    s_tolerance = 0.5 * (
        math.hypot(printed_unit(old[0]), printed_unit(old[1])) + math.hypot(printed_unit(new[0]), printed_unit(new[1]))
    )
    if not abs(math.hypot(re_old, im_old) - math.hypot(re_new, im_new)) <= s_tolerance * (1 + UNIT_SLACK):
        diff.forbid(f"{where}: |S| changed ({old[0]} {old[1]} became {new[0]} {new[1]})")
    if not abs(p_old - p_new) <= max(printed_unit(old[2]), printed_unit(new[2])) * (1 + UNIT_SLACK):
        diff.forbid(f"{where}: P {old[2]} became {new[2]}, more than one printed unit")
    if re_old == im_old == re_new == im_new == 0.0:
        return  # S = 0 has no phase to rotate
    if momentum is None:
        diff.forbid(f"{where}: no 'total momentum' line before the row")
        return
    found = channel_momentum(case, v, momentum)
    if found is None:
        diff.forbid(f"{where}: no momentum for v = {v} (a gas-phase level missing, or the channel closed)")
        return
    k, dk = found
    distance = case.z_sp - case.z_ana
    expected = -k * distance
    angle_tolerance = abs(distance) * dk + k * case.z_rounding
    rotated = complex(re_old, im_old) * cmath.exp(1j * expected)
    s_deviation = abs(complex(re_new, im_new) - rotated)
    s_rotation_tolerance = s_tolerance + abs(complex(re_old, im_old)) * angle_tolerance
    if not s_deviation <= s_rotation_tolerance * (1 + UNIT_SLACK):
        diff.forbid(
            f"{where}: S {old[0]} {old[1]} became {new[0]} {new[1]}, not S_old exp(-i k (Z_sp - Zana)) "
            f"(off by {s_deviation:.3g}, tolerance {s_rotation_tolerance:.3g})"
        )
    if re_old == 0.0 or re_new == 0.0:
        return  # the Fortran prints a zero phase when Re S = 0 (smatrix_combine)
    deviation = math.remainder(phase_new - phase_old - expected, 2.0 * math.pi)
    tolerance = 0.5 * (printed_unit(old[3]) + printed_unit(new[3])) + angle_tolerance
    diff.max_abs = max(diff.max_abs, abs(deviation))
    if not abs(deviation) <= tolerance * (1 + UNIT_SLACK):
        diff.forbid(
            f"{where}: the phase changed by {phase_new - phase_old:.6f}, not by -k (Z_sp - Zana) = {expected:.6f} "
            f"mod 2 pi (off by {deviation:.3g} rad, tolerance {tolerance:.3g})"
        )


def initial_phase_differences(old, new, case):
    """The Differences of two versions of an on-the-fly fort.25 under class initial_phase: only the rows
    m1 = m2 = 0 below a block header of the initial state (v any, JIn, MJIn, diffraction order 0) may change,
    and each of them, changed or not, as initial_row_change allows; every other line must be identical."""
    diff = Differences()
    a, b = old.decode().splitlines(), new.decode().splitlines()
    if len(a) != len(b):
        diff.line_counts = f"{len(a)} lines became {len(b)}"
        diff.forbid(diff.line_counts)
        return diff
    state, momentum = None, None
    for number, (x, y) in enumerate(zip(a, b, strict=True), start=1):
        header = FORT25_HEADER.fullmatch(x)
        if header:
            state = tuple(int(group) for group in header.groups())
        found = FORT25_MOMENTUM.fullmatch(x)
        if found:
            momentum = found.group(1)
        row_old, row_new = FORT25_ROW.fullmatch(x), FORT25_ROW.fullmatch(y)
        initial = state is not None and state[1:] == (case.j_in, case.mj_in, 0)
        if x == y:
            if initial and row_old and row_old.group(1, 2) == ("0", "0"):  # unchanged: must be the rotation too
                initial_row_change(
                    diff, case, f"line {number}", state[0], momentum, row_old.groups()[2:], row_old.groups()[2:]
                )
            continue
        diff.lines += 1
        if not (initial and row_old and row_new and row_old.group(1, 2) == ("0", "0") == row_new.group(1, 2)):
            diff.forbid(f"line {number}: {x.strip()!r} became {y.strip()!r}, not a row of the initial state")
            continue
        initial_row_change(diff, case, f"line {number}", state[0], momentum, row_old.groups()[2:], row_new.groups()[2:])
    return diff


def read_records(data):
    """The records of a gfortran unformatted sequential file: a 4-byte little-endian length before
    and after each. Records over 2 GiB are split into subrecords (negative lengths), which no
    reference file has: refused. extract_pes.py reads the .pes files with it too."""
    records, position = [], 0
    while position < len(data):
        if position + 4 > len(data):
            raise FormatError(f"a truncated record marker at byte {position}")
        (length,) = struct.unpack_from("<i", data, position)
        end = position + 4 + length
        if length < 0:
            raise FormatError(f"subrecords (records over 2 GiB) at byte {position} are not supported")
        if end + 4 > len(data) or struct.unpack_from("<i", data, end)[0] != length:
            raise FormatError(f"the record at byte {position} has no matching end marker")
        records.append(data[position + 4 : end])
        position = end + 4
    return records


def coefficients(record):
    """The complex BK coefficients of one record (pairs of little-endian doubles)."""
    if len(record) % 16:
        raise FormatError(f"a BK record of {len(record)} bytes is not a list of complex numbers")
    values = struct.unpack(f"<{len(record) // 8}d", record)
    return [complex(values[k], values[k + 1]) for k in range(0, len(values), 2)]


def gate4_scale(records):
    """max(the largest coefficient modulus of the BK records, BK_SCALE_FLOOR); NaN if a
    coefficient is not finite (Python's max would skip a NaN)."""
    largest = 0.0
    for record in records:
        for modulus in map(abs, coefficients(record)):
            if not math.isfinite(modulus):
                return math.nan
            largest = max(largest, modulus)
    return max(largest, BK_SCALE_FLOOR)


def record_changes(diff, number, old, new):
    """Adds the coefficient changes of two versions of BK record number (1-based) that differ. A
    non-finite coefficient that is not bitwise the old one is never allowed."""
    diff.lines += 1
    for k, (p, q) in enumerate(zip(coefficients(old), coefficients(new), strict=True)):
        if p == q:
            continue
        if not (math.isfinite(abs(p)) and math.isfinite(abs(q))):
            if old[16 * k : 16 * k + 16] != new[16 * k : 16 * k + 16]:
                diff.numbers += 1
                diff.forbid(f"BK record {number}, coefficient {k + 1}: {p} became {q}, not finite")
            continue
        diff.numbers += 1
        diff.max_abs = max(diff.max_abs, abs(p - q))


def bk_differences(old, new, rule):
    """The Differences of two BK coefficient files: the header records must be identical, each
    coefficient |new - old| within limit times the gate-4 scale of old (class bk_scaled)."""
    diff = Differences(kind="bk")
    a, b = read_records(old), read_records(new)
    if [len(r) for r in a] != [len(r) for r in b]:
        raise FormatError(f"the record count or lengths differ ({len(a)} records became {len(b)})")
    if a[:BK_HEADER_RECORDS] != b[:BK_HEADER_RECORDS]:
        diff.forbid("the header records differ")
    diff.scale = gate4_scale(a[BK_HEADER_RECORDS:])
    if not math.isfinite(diff.scale):
        diff.forbid("OLD holds a coefficient that is not finite, so the gate-4 scale is undefined")
    for number, (x, y) in enumerate(zip(a[BK_HEADER_RECORDS:], b[BK_HEADER_RECORDS:], strict=True), start=1):
        if x != y:
            record_changes(diff, number, x, y)
    limit = rule.limit * diff.scale if rule.allowed == "bk_scaled" else 0.0
    if not diff.max_abs <= limit:  # written so that a NaN limit fails
        diff.forbid(f"max |Δc| {diff.max_abs:.3g} exceeds {limit:.3g}")
    return diff


def file_kind(path, data):
    """bk (BK_FILE), text (UTF-8) or binary."""
    if PurePosixPath(path).name == BK_FILE:
        return "bk"
    try:
        data.decode()
        return "text"
    except UnicodeDecodeError:
        return "binary"


def file_differences(path, old, new, rule, case=None):
    """The Differences of the two versions (bytes) of a file under its rule; case: the InitialPhaseCase of
    class initial_phase."""
    kind = file_kind(path, old)
    if file_kind(path, new) != kind:
        return Differences(kind="binary", violation=f"the {kind} file became {file_kind(path, new)}")
    if rule.allowed == "initial_phase":
        if kind != "text":
            return Differences(kind="binary", violation=f"class initial_phase does not apply to a {kind} file")
        return initial_phase_differences(old, new, case)
    if (rule.allowed == "bk_scaled" and kind != "bk") or (rule.allowed in TEXT_CLASSES and kind != "text"):
        return Differences(
            kind="binary", violation=f"class {rule.allowed} does not apply to a {kind} file; fix the row"
        )
    try:
        if kind == "bk":
            diff = bk_differences(old, new, rule)
        elif kind == "text":
            diff = text_differences(old, new, rule)
        else:
            diff = Differences(kind="binary", violation="binary files may only be identical")
    except FormatError as error:
        diff = Differences(kind="binary", violation=f"BK records: {error}")
    if rule.allowed in ("identical", "removed", "added"):
        diff.forbid("no difference is allowed")
    if rule.allowed == "unread":
        diff.violation = ""
    return diff


def classify(old_root, new_root, path, rule):
    """(allowed, text) for one path of either tree: whether its state is allowed, and the line to
    print ("" for an identical file)."""
    old, new = Path(old_root) / path, Path(new_root) / path
    if not new.exists():
        if rule.allowed == "removed":
            return True, "removed"
        return False, "missing in NEW; not allowed: only class removed lets a file go"
    if not old.exists():
        return added(new_root, path, rule)
    old_bytes, new_bytes = old.read_bytes(), new.read_bytes()
    if old_bytes == new_bytes:
        return True, ""
    case = None
    if rule.allowed == "initial_phase":
        try:
            case = initial_phase_case(new_root, path)
        except (Refusal, OSError, ValueError) as error:
            return False, f"class initial_phase needs the case's inputs; not allowed: {error}"
    diff = file_differences(path, old_bytes, new_bytes, rule, case)
    # A text file reports a change of its line count; a binary or BK file its size (e.g. truncated).
    resized = diff.kind != "text" and len(old_bytes) != len(new_bytes)
    sizes = f" ({len(old_bytes)} bytes became {len(new_bytes)})" if resized else ""
    return not diff.violation, diff.summary() + sizes


def added(new_root, path, rule):
    """(allowed, text) for a file that only NEW has: allowed under a row of class added when it is a committed
    input of its case (is_input: inputs have no MANIFEST line) or its case's MANIFEST.tsv in NEW lists it."""
    if rule.allowed != "added":
        return False, "added in NEW; not allowed: only class added lets a file appear"
    located = case_and_relative_path(path)
    if located is None:
        return False, "added in NEW outside a case directory; not allowed: add it to a case"
    case_dir, relative = located
    if is_input(relative):
        return True, "added (input)"
    manifest = Path(new_root) / case_dir / MANIFEST_NAME
    if not manifest.is_file() or not any(row["path"] == relative for row in read_tsv(manifest, MANIFEST_COLUMNS)):
        return False, "added in NEW without a line in its case's MANIFEST.tsv; generate it with its run_*.sh script"
    return True, "added"


def case_and_relative_path(path):
    """(case directory, path inside it) of a path of a reference tree (relative to the repository), or None
    when the path lies in no case directory."""
    for tree in REFERENCE_TREES:
        if path.startswith(tree + "/"):
            case, _, relative = path[len(tree) + 1 :].partition("/")
            return (f"{tree}/{case}", relative) if relative else None
    return None


def compare_trees(old_root, new_root, rules):
    """Prints one line per file that is not identical; returns (counts per outcome, failures,
    rules that matched no file)."""
    counts, failures, used = {"identical": 0}, 0, set()
    for path in sorted(reference_files(old_root) | reference_files(new_root)):
        index, rule = matching_rule(rules, path)
        used.add(index)
        allowed, text = classify(old_root, new_root, path, rule)
        if not text:
            counts["identical"] += 1
            continue
        outcome = rule.allowed if allowed else "not allowed"
        counts[outcome] = counts.get(outcome, 0) + 1
        failures += not allowed
        print(f"  {'ok  ' if allowed else 'FAIL'}  {rule.allowed:<13} {path}: {text}")
    unused = [rule for index, rule in enumerate(rules) if index not in used]
    return counts, failures, unused


def case_dirs(root):
    """The case directories of the reference trees under root (NOT_A_CASE excluded)."""
    return sorted(
        d for tree in REFERENCE_TREES for d in (Path(root) / tree).glob("*") if d.is_dir() and d.name not in NOT_A_CASE
    )


def generated_files(case_dir):
    """The case-relative paths of the files of a case directory that a Fortran run generated."""
    paths = (f.relative_to(case_dir).as_posix() for f in case_dir.rglob("*") if f.is_file())
    return {p for p in paths if p != MANIFEST_NAME and not is_input(p)}


def buildinfo_copy(root, install_id, cache):
    """The rows (key -> value) of the committed BUILDINFO copy of an install, or None if absent."""
    if install_id not in cache:
        path = Path(root) / BUILDINFO_COPIES / f"{install_id}.tsv"
        cache[install_id] = (
            {row["key"]: row["value"] for row in read_tsv(path, BUILDINFO_COLUMNS)} if path.is_file() else None
        )
    return cache[install_id]


@functools.lru_cache(maxsize=None)
def bk_file_sha(case_dir):
    """The sha256 (12 digits) of a case's BK file, "" if it has none; hashed once per case."""
    bk_file = case_dir / BK_FILE
    return sha256_of(bk_file)[:SHORT_HASH] if bk_file.is_file() else ""


def input_problems(case_dir, row):
    """The BK input of one MANIFEST line: a fine-analysis file must name its input fort.22, and
    that must be the case's BK file (a fine analysis of a stale BK file would pass unseen)."""
    inputs = dict(entry.partition(":")[::2] for entry in row["input_sha"].split(","))
    if FINE_INPUT in inputs:
        same = inputs[FINE_INPUT] == bk_file_sha(case_dir)
        return [] if same else [f"its input {FINE_INPUT} is not the case's {BK_FILE}"]
    fine = row["path"].startswith(FINE_DIR_PREFIX)
    return [f"a fine-analysis file without input {FINE_INPUT}"] if fine else []


def program_problems(root, row, buildinfos):
    """(install id or None, problems) of the fortran_ref of one MANIFEST line: an install with a
    committed BUILDINFO copy that has the executable the line names (buildinfos: a cache)."""
    match = re.fullmatch(r"(\S+) (\S+):([0-9a-f]{12})", row["fortran_ref"])
    if not match:
        return None, [f"fortran_ref {row['fortran_ref']!r} names no install and program"]
    install_id, program, program_sha = match.groups()
    buildinfo = buildinfo_copy(root, install_id, buildinfos)
    if buildinfo is None:
        return install_id, [f"install {install_id} has no {BUILDINFO_COPIES}/{install_id}.tsv"]
    if program == FORTRAN_EXECUTABLE and buildinfo.get("exe.sha256", "")[:SHORT_HASH] != program_sha:
        return install_id, [f"{program}:{program_sha} is not the executable of {install_id} (exe.sha256)"]
    return install_id, []


def line_problems(root, case_dir, row, buildinfos):
    """(install id or None, problems) of one MANIFEST line of a regenerated tree."""
    problems, target = [], case_dir / row["path"]
    if not target.is_file() or sha256_of(target) != row["sha256"]:
        problems.append("missing, or its sha256 differs from the MANIFEST")
    if row["runtime_drift"] != EMPTY:
        problems.append(f"made with run-time drift {row['runtime_drift']}")
    install_id, program = program_problems(root, row, buildinfos)
    return install_id, problems + input_problems(case_dir, row) + program


def provenance_problems(root, install):
    """(problems, number of MANIFEST lines) of the MANIFESTs of a regenerated tree, whose lines
    must all name install (the expected file's)."""
    problems, lines, buildinfos, others = [], 0, {}, {}
    for case_dir in case_dirs(root):
        case = case_dir.relative_to(root).as_posix()
        manifest = case_dir / MANIFEST_NAME
        if not manifest.is_file():
            problems.append(f"{case}: no {MANIFEST_NAME}")
            continue
        rows = read_tsv(manifest, MANIFEST_COLUMNS)
        lines += len(rows)
        for path in sorted(generated_files(case_dir) - {row["path"] for row in rows}):
            problems.append(f"{case}/{path}: generated, but not in the MANIFEST")
        for row in rows:
            install_id, found = line_problems(root, case_dir, row, buildinfos)
            problems += [f"{case}/{row['path']}: {p}" for p in found]
            if install_id not in (None, install):
                others[install_id] = others.get(install_id, 0) + 1
    # The expected differences were measured with one install: another one's output is unchecked.
    problems += [
        f"{count} MANIFEST lines name install {other}, not {install} (the expected file's)"
        for other, count in sorted(others.items())
    ]
    return problems, lines


def compare(args):
    """compare: prints the classification and the provenance problems; refuses unless all is
    allowed."""
    for root in (args.old, args.new):
        if not all((Path(root) / tree).is_dir() for tree in REFERENCE_TREES):
            raise Refusal(f"{root} lacks one of {', '.join(REFERENCE_TREES)}; name a repository checkout")
    install, rules = read_expected(args.expected)
    print(f"* Comparing the reference files of {args.old} (OLD) and {args.new} (NEW) under {args.expected}")
    counts, failures, unused = compare_trees(args.old, args.new, rules)
    for rule in unused:
        print(f"  note: the row {rule.glob} ({rule.allowed}) matches no file")
    provenance = "not checked (--no-provenance)"
    if not args.no_provenance:
        problems, lines = provenance_problems(args.new, install)
        for problem in problems:
            print(f"  FAIL  provenance    {problem}")
        failures += len(problems)
        provenance = f"{lines} MANIFEST lines checked against install {install}, {len(problems)} problems"
    files = sum(counts.values())
    outcomes = ", ".join(f"{name} {count}" for name, count in sorted(counts.items()))
    print(f"* {files} files: {outcomes}; NEW's provenance: {provenance}")
    if failures:
        raise Refusal(
            f"{failures} differences or provenance problems are not allowed (FAIL lines above); "
            f"fix the regeneration, or review the differences and extend {args.expected}"
        )
    print("* Every difference is allowed")


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
    add.add_argument(
        "--binary", metavar="FILE", help="the program that ran: the install's executable, or the CRP or SCM driver"
    )
    add.add_argument(
        "--input",
        metavar="PATH",
        action="append",
        default=[],
        help="an input as run (repeatable): a file, or a directory (a data set, hashed as "
        "a tree); named by the path as given if relative, else by its last component",
    )
    add.add_argument(
        "--runtime-drift",
        metavar="TEXT",
        default=EMPTY,
        help="run-time packages that differ from BUILDINFO.tsv (default: -)",
    )
    add.add_argument(
        "--pre-manifest",
        action="store_true",
        help=f'files made before MANIFESTs existed: fortran_ref "{PRE_MANIFEST_REF}", "-" '
        "for the columns not recorded then",
    )

    remove = actions.add_parser("delete", help="delete the listed files and the MANIFESTs")
    remove.add_argument("case_dirs", metavar="CASE_DIR", nargs="+", help="a case directory with a MANIFEST.tsv")

    check = commands.add_parser("compare", help="compare two reference trees against the expected differences")
    check.add_argument("old", metavar="OLD", help="repository checkout with the committed references")
    check.add_argument("new", metavar="NEW", help="repository checkout with the regenerated references")
    check.add_argument(
        "--expected",
        metavar="FILE",
        required=True,
        help="the expected-differences file, e.g. tools/fortran_reference/expected/r1713.tsv",
    )
    check.add_argument("--no-provenance", action="store_true", help="skip the MANIFEST and BUILDINFO checks of NEW")

    args = parser.parse_args(argv)
    if args.command == "manifest" and args.action == "record":
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
        if args.command == "compare":
            compare(args)
        elif args.action == "record":
            record(args)
        else:
            delete(args)
    except Refusal as refusal:
        sys.stdout.flush()  # the report first, also in a log of stdout and stderr
        print(f"error: {refusal}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

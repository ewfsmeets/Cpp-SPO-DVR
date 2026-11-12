#!/usr/bin/env python3
"""Check that the documentation agrees with the code and the repository.

Checks (each finding is printed as `file:line: <check>: <what>; <how to fix>`):
  input-keys    keys read by src/input/json_input.cpp <-> docs/input_schema.md (both directions)
  output-files  file names written by src/app, src/io, src/report <-> docs/output_files.md
  links         tracked *.md: relative links, images and #anchors resolve (GitHub slugs);
                `file.ext:N` / `path/file.ext:N-M` anchors name a repository file with >= N lines
  decision-ids  `flux Qn`, `SCM Dn`, `CRP Dn`, `absorber Dn`, `low-Z Dn` cited in a tracked file
                are defined in docs/fortran_deviations.md (a heading, list item or table row that
                starts with the ID); the unqualified "decision Dn" wording is reported
  local-paths   tracked text files hold no ~/work path other than the Fortran trees and no plan
                vocabulary (the thesis record docs/validation_thesis_fig5_3_dz/ is exempt)
  cli           the options of `spodvr --help` <-> the usage block of README "Running"
  gates         every GateN test name has gate N in docs/validation.md; every case directory
                tests/reference/<case>/ (one with a DVR.inp) is named there
  env-vars      every SPODVR_* environment variable read by code, tests, tools or CMake is in a
                living doc; retired variables appear only in their one allowed file. Report-only:
                variables documented in several living docs
  numbers       report-only: numbers with >= 3 significant digits or an exponent that occur in
                more than one living doc (README.md, docs/*.md, the other tracked README.md files
                outside the thesis record)

Inputs: the tracked files of --root (`git ls-files`; every file below it, without .git and build*
directories, if it is not the top of a git work tree), the allow-list and, for check cli, the
output of `<spodvr> --help`. The allow-list is a TSV file (`#` comments, header
`check file match reason`, rows sorted): a finding is accepted when a row has its check and file
and `match` is a substring of the finding line or of the offending source line. A row that is
the only match of no finding is stale and itself a finding.
Output (stdout): the report-only lists, the findings, a summary line with counts per check.
Exit status: 0 no findings, 1 findings, 2 bad arguments.

Usage: python3 -I -B tools/check_docs.py [--root DIR] [--spodvr PATH] [--allowlist FILE]
                                         [--checks NAME,...] [--no-report]
  --root DIR        repository to check (default: the repository that contains this script)
  --spodvr PATH     built program for check cli (default: none; check cli is skipped, with a note)
  --allowlist FILE  accepted findings (default: <root>/tools/check_docs_allowlist.tsv)
  --checks LIST     comma-separated subset of the checks above (default: all)
  --no-report       do not print the report-only lists (their sizes are still in the summary)
Environment: none read. Run time: under 1 s. ctest: Docs.Consistency. Docs: README.md "Testing".
"""

import argparse
import os
import re
import subprocess
import sys
import time
import unicodedata
import urllib.parse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

CHECKS = ("input-keys", "output-files", "links", "decision-ids", "local-paths", "cli", "gates",
          "env-vars", "numbers")
REFERENCE_DATA = ("tests/reference/", "tests/crplib/reference/", "potentials/")
THESIS_RECORD = "docs/validation_thesis_fig5_3_dz/"
# This checker and its allow-list spell out the patterns they look for.
SELF = ("tools/check_docs.py", "tools/check_docs_allowlist.tsv")
# Retired SPODVR_* environment variables: they may appear only in RETIRED_ENV_HOME.
RETIRED_ENV_VARS: tuple[str, ...] = ("SPODVR_FORTRAN_EXE", "SPODVR_FORTRAN_OBJS")
RETIRED_ENV_HOME = "tools/fortran_reference/reference_config.sh"
ALLOWLIST_HEADER = ["check", "file", "match", "reason"]
ALLOWABLE_CHECKS = tuple(c for c in CHECKS if c != "numbers")  # numbers is report-only


@dataclass(frozen=True)
class Finding:
    """One finding: where, which check, what is wrong (and the offending source line)."""
    check: str
    file: str
    line: int
    what: str
    text: str = ""

    def __str__(self):
        return f"{self.file}:{self.line}: {self.check}: {self.what}"


class Repo:
    """The tracked files of the repository and a cache of their text lines."""

    def __init__(self, root: Path):
        self.root = root
        self.files = list_files(root)
        self.file_set = set(self.files)
        self.dirs = {str(p) for f in self.files for p in PurePosixPath(f).parents if str(p) != "."}
        self._lines: dict[str, list[str] | None] = {}

    def lines(self, rel: str) -> list[str] | None:
        """The lines of a file, or None if it is missing or binary."""
        if rel not in self._lines:
            self._lines[rel] = read_text_lines(self.root / rel)
        return self._lines[rel]

    def text_files(self, skip: tuple[str, ...] = ()) -> list[str]:
        """Tracked text files outside the reference data, this checker and the `skip` prefixes."""
        excluded = REFERENCE_DATA + skip
        return [f for f in self.files if not f.startswith(excluded) and f not in SELF
                and self.lines(f) is not None]

    def living_docs(self) -> list[str]:
        """README.md, docs/*.md and the other tracked README.md files outside the thesis record."""
        return [f for f in self.files
                if f == "README.md" or (f.startswith("docs/") and f.count("/") == 1 and f.endswith(".md"))
                or (f.endswith("/README.md") and not f.startswith((THESIS_RECORD,) + REFERENCE_DATA))]


def git(root: Path, *args: str) -> str | None:
    """Output of a git command in root, or None if it fails."""
    try:
        return subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True).stdout.decode()
    except (OSError, subprocess.CalledProcessError):
        return None


def list_files(root: Path) -> list[str]:
    """`git ls-files` if root is the top of a git work tree, else every file below root
    (without .git and build* directories), as sorted relative POSIX paths."""
    top = git(root, "rev-parse", "--show-toplevel")
    if top is not None and Path(top.strip()).resolve() == root.resolve():
        return sorted(f for f in git(root, "ls-files", "-z").split("\0") if f and (root / f).is_file())
    files = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d != ".git" and not d.startswith("build")]
        files += [(Path(dirpath) / n).relative_to(root).as_posix() for n in filenames]
    return sorted(files)


def read_text_lines(path: Path) -> list[str] | None:
    try:
        data = path.read_bytes()
    except OSError:
        return None
    if b"\0" in data:
        return None
    try:
        return data.decode("utf-8").splitlines()
    except UnicodeDecodeError:
        return None


def line_of(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def markdown_lines(lines: list[str]):
    """(line number, line) of the lines outside fenced code blocks."""
    fence = None
    for number, line in enumerate(lines, 1):
        marker = re.match(r"\s*(`{3,}|~{3,})", line)
        if marker:
            char = marker.group(1)[0]
            if fence is None:
                fence = char
            elif fence == char:
                fence = None
            continue
        if fence is None:
            yield number, line


CODE_SPAN = re.compile(r"`([^`]+)`")


def without_code_spans(line: str) -> str:
    return CODE_SPAN.sub(lambda m: " " * len(m.group(0)), line)


# ---------------------------------------------------------------------------------------------
# 1. Input keys
# ---------------------------------------------------------------------------------------------
INPUT_SOURCE = "src/input/json_input.cpp"
INPUT_SCHEMA = "docs/input_schema.md"
KEY_FORMS = (re.compile(r'\bin\.\w+\(\s*"(\w+)"'),        # in.read("Key", …), in.contains("Key")
             re.compile(r'\bInputError\(\s*"(\w+)"'))     # errors name the key they concern
KEY_LIST = re.compile(r'\{\s*"\w+"(?:\s*,\s*"\w+")*\s*\}')  # {"KeyA", "KeyB"} lists of keys
STRING = re.compile(r'"(\w+)"')


def keys_read(text: str) -> dict[str, int]:
    """Every key json_input.cpp reads, with the line of its first occurrence."""
    keys: dict[str, int] = {}
    for form in KEY_FORMS:
        for m in form.finditer(text):
            keys.setdefault(m.group(1), line_of(text, m.start(1)))
    for block in KEY_LIST.finditer(text):
        for m in STRING.finditer(block.group(0)):
            keys.setdefault(m.group(1), line_of(text, block.start() + m.start(1)))
    return dict(sorted(keys.items(), key=lambda item: item[1]))


def table_rows(lines: list[str]):
    """(line number, cells, header cells) of every markdown table body row."""
    header = None
    for number, line in markdown_lines(lines):
        if not line.lstrip().startswith("|"):
            header = None
            continue
        cells = [c.strip() for c in re.split(r"(?<!\\)\|", line.strip())[1:-1]]
        if header is None:
            header = cells
        elif not all(re.fullmatch(r":?-+:?", c) for c in cells):
            yield number, cells, header


def check_input_keys(repo: Repo) -> list[Finding]:
    source, schema = repo.lines(INPUT_SOURCE), repo.lines(INPUT_SCHEMA)
    missing = [f for f, lines in ((INPUT_SOURCE, source), (INPUT_SCHEMA, schema)) if lines is None]
    if missing:
        return [Finding("input-keys", f, 1, "file not found; restore it or update tools/check_docs.py")
                for f in missing]
    read = keys_read("\n".join(source))
    mentioned = {m.group(1).strip() for line in schema for m in CODE_SPAN.finditer(line)}
    findings = [Finding("input-keys", INPUT_SOURCE, line,
                        f'key "{key}" is read but not documented in {INPUT_SCHEMA}; add it to a key table there',
                        source[line - 1])
                for key, line in read.items() if key not in mentioned]
    for number, cells, header in table_rows(schema):
        if header[0] != "Key":
            continue
        for m in CODE_SPAN.finditer(cells[0]):
            if m.group(1) not in read:
                findings.append(Finding("input-keys", INPUT_SCHEMA, number,
                                        f'key "{m.group(1)}" is documented but not read by {INPUT_SOURCE}; '
                                        "remove the row or read the key", schema[number - 1]))
    return findings


# ---------------------------------------------------------------------------------------------
# 2. Output files
# ---------------------------------------------------------------------------------------------
OUTPUT_DOC = "docs/output_files.md"
WRITER_DIRS = ("src/app/", "src/io/", "src/report/")
OUTPUT_EXTENSIONS = "dat|tsv|json|bin|out|txt"
WRITTEN_NAME = re.compile(r'"((?:[\w-]+/)*[\w-]+\.(?:' + OUTPUT_EXTENSIONS + r')|fort\.\d+)"')
DOC_NAME = re.compile(r"(?:[\w*-]+/)*(?:[\w*.-]+\.(?:" + OUTPUT_EXTENSIONS + r"|tmp)|fort\.\d+)")
DOC_DIR = re.compile(r"(?:[\w-]+/)+")
ANY_STRING = re.compile(r'"([^"\\\n]*)"')


def written_names(repo: Repo) -> tuple[dict[str, tuple[str, int, str]], set[str]]:
    """Output file names in the writer sources (with a location) and all their string literals.

    Lines that look up restart input (`first_existing(`) name files that are read, not written.
    """
    names: dict[str, tuple[str, int, str]] = {}
    literals: set[str] = set()
    for f in repo.files:
        if not f.startswith(WRITER_DIRS) or not f.endswith((".cpp", ".hpp")):
            continue
        for number, line in enumerate(repo.lines(f) or [], 1):
            literals.update(ANY_STRING.findall(line))
            if "first_existing(" in line:
                continue
            for m in WRITTEN_NAME.finditer(line):
                names.setdefault(m.group(1), (f, number, line))
    return names, literals


def check_output_files(repo: Repo) -> list[Finding]:
    doc = repo.lines(OUTPUT_DOC)
    if doc is None:
        return [Finding("output-files", OUTPUT_DOC, 1, "file not found; restore it or update tools/check_docs.py")]
    names, literals = written_names(repo)
    text = "\n".join(doc)
    findings = []
    for name, (f, number, line) in names.items():
        if not re.search(r"(?<![\w.-])" + re.escape(name) + r"(?![\w-])", text):
            findings.append(Finding("output-files", f, number,
                                    f'output file "{name}" is written but not named in {OUTPUT_DOC}; document it there',
                                    line))
    written = {PurePosixPath(n).name for n in names}
    for number, line in enumerate(doc, 1):
        for m in CODE_SPAN.finditer(line):
            problem = documented_name_problem(m.group(1).strip(), written, literals)
            if problem:
                findings.append(Finding("output-files", OUTPUT_DOC, number,
                                        f"{problem}; remove it from the doc or fix the name", line))
    return findings


def documented_name_problem(span: str, written: set[str], literals: set[str]) -> str | None:
    """Why a file or directory named in output_files.md is not written by the program, if it is not."""
    if DOC_DIR.fullmatch(span):
        absent = [d for d in span.strip("/").split("/") if d not in literals]
        return f'directory "{span}" is documented but "{absent[0]}" is not in the writer sources' if absent else None
    if not DOC_NAME.fullmatch(span) or "*" in span:
        return None
    name = PurePosixPath(span).name
    if name.endswith(".tmp"):
        if ".tmp" in literals and name[:-4] in written:
            return None
    elif name in written:
        return None
    return f'file "{span}" is documented but not written by {", ".join(WRITER_DIRS)}'


# ---------------------------------------------------------------------------------------------
# 3. Links and line anchors in tracked markdown files
# ---------------------------------------------------------------------------------------------
LINK = re.compile(r'!?\[[^\]]*\]\(\s*<?([^)\s>]*)>?(?:\s+"[^"]*")?\s*\)')
HTML_LINK = re.compile(r'\b(?:src|href)\s*=\s*"([^"]*)"')
REFERENCE_LINK = re.compile(r"^\s{0,3}\[[^\]]+\]:\s*<?([^\s>]+)>?")
LINE_ANCHOR = re.compile(r"(?<![\w./~:-])((?:[\w.~-]+/)*[\w-][\w.-]*\."
                         r"(?:cpp|hpp|h|cu|cuh|py|sh|cmake|md|f90|F90|f|inp|json|txt)):(\d+)(?:-(\d+))?(?!\d)")
HEADING = re.compile(r"^\s{0,3}#{1,6}\s+(.*?)(?:\s+#+)?\s*$")
HTML_ANCHOR = re.compile(r'<a\s+(?:name|id)\s*=\s*"([^"]+)"')


def github_slug(heading: str) -> str:
    """The anchor GitHub gives a heading: lower case, punctuation dropped, spaces to hyphens."""
    text = re.sub(r"!?\[([^\]]*)\]\([^)]*\)", r"\1", heading)
    text = re.sub(r"<[^>]+>", "", text).replace("`", "").strip().lower()
    kept = (c for c in text if c.isalnum() or c in "-_ " or unicodedata.category(c).startswith("M"))
    return "".join(kept).replace(" ", "-")


def heading_anchors(lines: list[str]) -> set[str]:
    anchors: set[str] = set()
    counts: dict[str, int] = defaultdict(int)
    for _, line in markdown_lines(lines):
        anchors.update(HTML_ANCHOR.findall(line))
        m = HEADING.match(line)
        if m:
            slug = github_slug(m.group(1))
            anchors.add(slug if counts[slug] == 0 else f"{slug}-{counts[slug]}")
            counts[slug] += 1
    return anchors


def link_targets(lines: list[str]):
    """(line number, target, line) of every link, image and reference definition outside code."""
    for number, line in markdown_lines(lines):
        bare = without_code_spans(line)
        for pattern in (LINK, HTML_LINK, REFERENCE_LINK):
            for m in pattern.finditer(bare):
                yield number, m.group(1), line


def link_problem(repo: Repo, md: str, target: str, anchors: dict[str, set[str]]) -> str | None:
    """Why a link target of `md` does not resolve, if it does not."""
    if not target or re.match(r"[A-Za-z][\w+.-]*:|//", target):
        return None
    path, _, anchor = target.partition("#")
    path = urllib.parse.unquote(path)
    if path:
        base = PurePosixPath() if path.startswith("/") else PurePosixPath(md).parent
        resolved = os.path.normpath((base / path.lstrip("/")).as_posix())
        if resolved not in repo.file_set and resolved not in repo.dirs:
            return f'link target "{target}" is not a tracked file or directory'
    else:
        resolved = md
    if anchor and resolved.endswith(".md") and resolved in repo.file_set:
        if resolved not in anchors:
            anchors[resolved] = heading_anchors(repo.lines(resolved) or [])
        if anchor not in anchors[resolved]:
            return f'anchor "#{anchor}" is not a heading of {resolved}'
    return None


def anchor_problem(repo: Repo, md: str, path: str, first: int, last: int | None) -> str | None:
    """Why a `file.ext:N[-M]` anchor does not name an existing repository file with enough lines."""
    if "/" in path:
        candidates = [os.path.normpath(p) for p in (path, (PurePosixPath(md).parent / path).as_posix())]
        matches = [p for p in candidates if p in repo.file_set][:1]
    else:
        matches = [f for f in repo.files if PurePosixPath(f).name == path]
    if not matches:
        return (f'"{path}" is not a file of the repository; drop the line number '
                "(Fortran anchors are file + routine) or fix the name")
    if len(matches) > 1:
        return f'"{path}" matches {len(matches)} files ({", ".join(matches[:3])}, …); write its path'
    count = len(repo.lines(matches[0]) or [])
    needed = max(first, last or 0)
    if count < needed:
        return f"{matches[0]} has {count} lines, fewer than {needed}; update the line number"
    return None


def check_links(repo: Repo) -> list[Finding]:
    findings = []
    anchors: dict[str, set[str]] = {}
    for md in (f for f in repo.files if f.endswith(".md")):
        lines = repo.lines(md) or []
        for number, target, line in link_targets(lines):
            problem = link_problem(repo, md, target, anchors)
            if problem:
                findings.append(Finding("links", md, number, f"{problem}; fix the link", line))
        for number, line in markdown_lines(lines):
            for m in LINE_ANCHOR.finditer(line):
                last = int(m.group(3)) if m.group(3) else None
                problem = anchor_problem(repo, md, m.group(1), int(m.group(2)), last)
                if problem:
                    findings.append(Finding("links", md, number, f"line anchor {m.group(0)}: {problem}", line))
    return findings


# ---------------------------------------------------------------------------------------------
# 4. Decision IDs
# ---------------------------------------------------------------------------------------------
DEVIATIONS = "docs/fortran_deviations.md"
ID_PREFIX = r"(flux Q|SCM D|CRP D|absorber D|low-Z D)"
ID_CITATION = re.compile(r"\b" + ID_PREFIX + r"(\d+)((?:/[DQ]?\d+)*)")
ID_DEFINITION = re.compile(r"^\s*(?:#{1,6}\s+|[-*]\s+|\|\s*)?(?:\*\*)?" + ID_PREFIX + r"(\d+)(?!\d)")
UNQUALIFIED_ID = re.compile(r"\bdecision (D\d+(?:/D?\d+)*)")


def defined_ids(repo: Repo) -> tuple[set[str], list[Finding]]:
    """IDs defined in fortran_deviations.md, and findings for IDs defined twice."""
    defined: dict[str, int] = {}
    findings = []
    for number, line in enumerate(repo.lines(DEVIATIONS) or [], 1):
        m = ID_DEFINITION.match(line)
        if not m:
            continue
        decision = m.group(1) + m.group(2)
        if decision in defined:
            findings.append(Finding("decision-ids", DEVIATIONS, number,
                                    f'"{decision}" is defined twice (first at line {defined[decision]}); '
                                    "give one of them a new ID", line))
        defined.setdefault(decision, number)
    return set(defined), findings


def check_decision_ids(repo: Repo) -> list[Finding]:
    defined, findings = defined_ids(repo)
    for f in repo.text_files():
        for number, line in enumerate(repo.lines(f), 1):
            for m in ID_CITATION.finditer(line):
                numbers = [m.group(2)] + re.findall(r"\d+", m.group(3))
                for decision in (m.group(1) + n for n in numbers):
                    if decision not in defined:
                        findings.append(Finding("decision-ids", f, number,
                                                f'"{decision}" is not defined in {DEVIATIONS}; define it there '
                                                "(a heading, list item or table row starting with the ID)", line))
            for m in UNQUALIFIED_ID.finditer(line):
                findings.append(Finding("decision-ids", f, number,
                                        f'unqualified decision ID "{m.group(0)}"; name its area (flux Qn, SCM Dn, '
                                        f"CRP Dn, absorber Dn, low-Z Dn) as defined in {DEVIATIONS}", line))
    return findings


# ---------------------------------------------------------------------------------------------
# 5. Local paths and plan vocabulary
# ---------------------------------------------------------------------------------------------
LOCAL_PATH = re.compile(r"(?:~|\$HOME|\$\{HOME\}|/home/[\w.-]+)/work(?:/([^\s`'\")\]>,;]*))?")
FORTRAN_TREE = re.compile(r"SPO-DVR(?:_v5\.19|_ref)?(?:/.*)?")
PLAN_WORDS = re.compile(r"(?<![\w.])\(new\)|\b[Ss]tage 0\b|\b[Rr]ound \d+\b"
                        r"|\b[Pp]hase [A-Z]?\d+[a-z]?(?!\w|\.\d)|\bR-D\d+")


def check_local_paths(repo: Repo) -> list[Finding]:
    findings = []
    for f in repo.text_files(skip=(THESIS_RECORD,)):
        for number, line in enumerate(repo.lines(f), 1):
            for m in LOCAL_PATH.finditer(line):
                if not FORTRAN_TREE.fullmatch((m.group(1) or "").rstrip(".:")):
                    findings.append(Finding("local-paths", f, number,
                                            f'local path "{m.group(0)}"; only the Fortran trees ~/work/SPO-DVR, '
                                            "~/work/SPO-DVR_v5.19, ~/work/SPO-DVR_ref may be named: drop it", line))
            for m in PLAN_WORDS.finditer(line):
                findings.append(Finding("local-paths", f, number,
                                        f'plan vocabulary "{m.group(0)}"; describe the content instead', line))
    return findings


# ---------------------------------------------------------------------------------------------
# 6. Command-line options
# ---------------------------------------------------------------------------------------------
USAGE_ITEM = re.compile(r"\[(--[\w-]+)([^\]]*)\]")


def usage_items(text: str) -> dict[str, str]:
    """Option name -> its usage item, whitespace collapsed (e.g. "--threads N")."""
    return {m.group(1): " ".join((m.group(1) + m.group(2)).split()) for m in USAGE_ITEM.finditer(text)}


def readme_usage(lines: list[str]) -> tuple[int, str]:
    """(line, text) of the first fenced block in README "Running" that starts with `spodvr `."""
    in_running, block, start = False, None, 0
    for number, line in enumerate(lines, 1):
        if line.startswith("## "):
            in_running = line.strip() == "## Running"
        if not in_running:
            continue
        if line.lstrip().startswith("```"):
            if block is not None and block and block[0].lstrip().startswith("spodvr "):
                return start, "\n".join(block)
            block, start = ([], number + 1) if block is None else (None, 0)
        elif block is not None:
            block.append(line)
    return 0, ""


def help_usage(spodvr: str) -> str:
    """The usage lines of `spodvr --help` (from "Usage:" to the next blank line)."""
    result = subprocess.run([spodvr, "--help"], capture_output=True, text=True, timeout=60)
    if result.returncode != 0:
        raise RuntimeError(f"{spodvr} --help exited with {result.returncode}")
    lines = result.stdout.splitlines()
    start = next((i for i, line in enumerate(lines) if line.lstrip().startswith("Usage:")), None)
    if start is None:
        raise RuntimeError(f'{spodvr} --help prints no "Usage:" line')
    usage = []
    for line in lines[start:]:
        if not line.strip():
            break
        usage.append(line)
    return "\n".join(usage)


def check_cli(repo: Repo, spodvr: str) -> list[Finding]:
    readme = repo.lines("README.md") or []
    number, block = readme_usage(readme)
    if not block:
        return [Finding("cli", "README.md", 1, 'no usage block (fenced, starting with "spodvr ") in "## Running"; '
                        "add it")]
    documented, printed = usage_items(block), usage_items(help_usage(spodvr))
    findings = []
    for name in printed.keys() - documented.keys():
        findings.append(Finding("cli", "README.md", number,
                                f'option "{printed[name]}" of spodvr --help is missing from the usage block; add it'))
    for name in documented.keys() - printed.keys():
        findings.append(Finding("cli", "README.md", number,
                                f'option "{documented[name]}" is not an option of spodvr --help; remove it'))
    for name in sorted(documented.keys() & printed.keys()):
        if documented[name] != printed[name]:
            findings.append(Finding("cli", "README.md", number,
                                    f'usage block has "[{documented[name]}]", spodvr --help "[{printed[name]}]"; '
                                    "copy the program's form"))
    return sorted(findings, key=lambda f: f.what)


# ---------------------------------------------------------------------------------------------
# 7. Gates
# ---------------------------------------------------------------------------------------------
VALIDATION = "docs/validation.md"
TEST_NAME = re.compile(r"\bTEST(?:_F|_P)?\(\s*(\w+)\s*,\s*(\w+)\s*\)")


def gate_entries(lines: list[str]) -> set[int]:
    """Gate numbers with a numbered entry ("N. ") in the "## Gates" section."""
    gates, in_gates = set(), False
    for _, line in markdown_lines(lines):
        if line.startswith("## "):
            in_gates = line.strip() == "## Gates"
        m = re.match(r"(\d+)\.\s", line)
        if in_gates and m:
            gates.add(int(m.group(1)))
    return gates


def check_gates(repo: Repo) -> list[Finding]:
    doc = repo.lines(VALIDATION)
    if doc is None:
        return [Finding("gates", VALIDATION, 1, "file not found; restore it or update tools/check_docs.py")]
    entries, text = gate_entries(doc), "\n".join(doc)
    findings = []
    for f in (f for f in repo.files if f.startswith("tests/") and f.endswith(".cpp")):
        for number, line in enumerate(repo.lines(f) or [], 1):
            for m in TEST_NAME.finditer(line):
                for name in m.groups():
                    gate = re.match(r"Gate(\d+)", name)
                    if gate and int(gate.group(1)) not in entries:
                        findings.append(Finding("gates", f, number,
                                                f'test {name} names gate {gate.group(1)}, which has no entry in '
                                                f'"## Gates" of {VALIDATION}; add it or rename the test', line))
    # A case directory holds the Fortran input of its case (DVR.inp).
    cases = sorted({f.split("/")[2] for f in repo.files if re.fullmatch(r"tests/reference/[^/]+/DVR\.inp", f)})
    for case in cases:
        if not re.search(r"(?<![\w-])" + re.escape(case) + r"(?![\w-])", text):
            findings.append(Finding("gates", f"tests/reference/{case}", 1,
                                    f'reference directory "{case}" is not named in {VALIDATION}; describe it there'))
    return findings


# ---------------------------------------------------------------------------------------------
# 8. Environment variables
# ---------------------------------------------------------------------------------------------
ENV_READS = (
    ((".cpp", ".hpp", ".h", ".cu", ".cuh"), re.compile(r'getenv\(\s*"(SPODVR_\w+)"')),
    (("CMakeLists.txt", ".cmake"), re.compile(r"\$ENV\{(SPODVR_\w+)\}")),
    ((".sh",), re.compile(r"\$\{?(SPODVR_\w+)")),
    ((".py",), re.compile(r"""(?:environ(?:\.get)?\s*[\[(]|getenv\()\s*["'](SPODVR_\w+)""")),
)


def env_vars_read(repo: Repo) -> dict[str, tuple[str, int, str]]:
    """SPODVR_* variables read from the environment, with their first location."""
    found: dict[str, tuple[str, int, str]] = {}
    for f in repo.text_files():
        for suffixes, pattern in ENV_READS:
            if not f.endswith(suffixes):
                continue
            for number, line in enumerate(repo.lines(f), 1):
                for name in pattern.findall(line):
                    found.setdefault(name, (f, number, line))
    return found


def check_env_vars(repo: Repo) -> tuple[list[Finding], list[str]]:
    findings, report = [], []
    docs = {d: "\n".join(repo.lines(d) or []) for d in repo.living_docs()}
    for name, (f, number, line) in sorted(env_vars_read(repo).items()):
        where = [d for d, text in docs.items() if re.search(r"\b" + name + r"\b", text)]
        if not where:
            findings.append(Finding("env-vars", f, number,
                                    f"environment variable {name} is read but not documented in a living doc "
                                    "(README.md, docs/*.md, a README.md); document it", line))
        elif len(where) > 1:
            report.append(f"env-vars: {name} is documented in {len(where)} docs: {', '.join(where)}")
    for name in RETIRED_ENV_VARS:
        for f in repo.text_files():
            if f == RETIRED_ENV_HOME:
                continue
            for number, line in enumerate(repo.lines(f), 1):
                if re.search(r"\b" + name + r"\b", line):
                    findings.append(Finding("env-vars", f, number,
                                            f"retired variable {name} may appear only in {RETIRED_ENV_HOME}; "
                                            "drop it", line))
    return findings, report


# ---------------------------------------------------------------------------------------------
# 9. Repeated numbers (report-only)
# ---------------------------------------------------------------------------------------------
NUMBER = re.compile(r"(?<![\w.\-/−])(\d+(?:\.\d+)?)([eE][-+−]?\d+)?(?![\w]|\.\d|-\d)")


def significant_digits(mantissa: str) -> int:
    digits = mantissa.replace(".", "").lstrip("0")
    return len(digits if "." in mantissa else digits.rstrip("0"))


def check_numbers(repo: Repo) -> list[str]:
    places: dict[str, dict[str, list[int]]] = defaultdict(lambda: defaultdict(list))
    for doc in repo.living_docs():
        for number, line in markdown_lines(repo.lines(doc) or []):
            line = re.sub(r"\d{4}-\d\d-\d\d|\b[0-9a-f]{7,40}\b|https?://\S+", " ", line)
            for m in NUMBER.finditer(line):
                if m.group(2) or significant_digits(m.group(1)) >= 3:
                    places[m.group(0).replace("−", "-").lower()][doc].append(number)
    report = []
    for value, docs in sorted(places.items()):
        if len(docs) > 1:
            where = ", ".join(f"{d}:{','.join(map(str, lines))}" for d, lines in sorted(docs.items()))
            report.append(f"numbers: {value} in {where}")
    return report


# ---------------------------------------------------------------------------------------------
# Allow-list
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class AllowEntry:
    """One allow-list row: accept findings of `check` in `file` whose line or source contains `match`."""
    check: str
    file: str
    match: str
    reason: str
    line: int

    def accepts(self, finding: Finding) -> bool:
        return (self.check == finding.check and self.file == finding.file
                and (self.match in str(finding) or self.match in finding.text))


def read_allowlist(path: Path, label: str) -> tuple[list[AllowEntry], list[Finding]]:
    """The allow-list rows, and findings for malformed or unsorted rows."""
    entries, problems = [], []
    lines = read_text_lines(path) or []
    rows = [(n, line) for n, line in enumerate(lines, 1) if line.strip() and not line.startswith("#")]
    if not rows or rows[0][1].split("\t") != ALLOWLIST_HEADER:
        problems.append(Finding("allowlist", label, rows[0][0] if rows else 1,
                                "the first row must be the header check<TAB>file<TAB>match<TAB>reason"))
    for number, line in rows[1:]:
        cells = line.split("\t")
        if len(cells) != 4 or not all(c.strip() for c in cells) or cells[0] not in ALLOWABLE_CHECKS:
            problems.append(Finding("allowlist", label, number,
                                    "row needs 4 non-empty tab-separated cells and a check from "
                                    f"{', '.join(ALLOWABLE_CHECKS)}", line))
            continue
        entries.append(AllowEntry(*cells, line=number))
    keys = [(e.check, e.file, e.match) for e in entries]
    for previous, entry, key in zip(keys, entries[1:], keys[1:]):
        if key < previous:
            problems.append(Finding("allowlist", label, entry.line,
                                    "rows are not sorted (by check, file, match); sort them", ""))
            break
    return entries, problems


def apply_allowlist(findings: list[Finding], entries: list[AllowEntry], label: str):
    """(remaining findings, number accepted, findings for stale rows).

    A row is stale unless it is the only row that accepts some finding: then deleting any row
    makes the check fail, and rows that suppress nothing cannot accumulate.
    """
    remaining, accepted, needed = [], 0, set()
    for finding in findings:
        matching = [e for e in entries if e.accepts(finding)]
        if not matching:
            remaining.append(finding)
            continue
        accepted += 1
        if len(matching) == 1:
            needed.add(matching[0])
    stale = [Finding("allowlist", label, e.line,
                     f'row "{e.check}  {e.file}  {e.match}" suppresses nothing (or only findings another row '
                     "also suppresses); delete it")
             for e in entries if e not in needed]
    return remaining, accepted, stale


# ---------------------------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------------------------
def parse_arguments(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent,
                        help="repository to check (default: the repository of this script)")
    parser.add_argument("--spodvr", help="built spodvr program for check cli (default: skip check cli)")
    parser.add_argument("--allowlist", type=Path,
                        help="accepted findings (default: <root>/tools/check_docs_allowlist.tsv)")
    parser.add_argument("--checks", default=",".join(CHECKS), help="comma-separated checks (default: all)")
    parser.add_argument("--no-report", action="store_true", help="do not print the report-only lists")
    args = parser.parse_args(argv)
    args.checks = [c.strip() for c in args.checks.split(",") if c.strip()]
    unknown = [c for c in args.checks if c not in CHECKS]
    if unknown:
        parser.error(f"unknown check(s) {', '.join(unknown)}; choose from {', '.join(CHECKS)}")
    if not args.root.is_dir():
        parser.error(f"--root {args.root} is not a directory")
    args.root = args.root.resolve()
    args.allowlist = args.allowlist or args.root / "tools/check_docs_allowlist.tsv"
    if not args.allowlist.is_file():
        parser.error(f"--allowlist {args.allowlist} is not a file")
    if args.spodvr and not os.access(args.spodvr, os.X_OK):
        parser.error(f"--spodvr {args.spodvr} is not an executable file")
    return args


def run_checks(repo: Repo, args: argparse.Namespace) -> tuple[dict[str, list[Finding]], list[str], list[str]]:
    """Findings per check, report-only lines and notes."""
    findings: dict[str, list[Finding]] = {}
    report, notes = [], []
    simple = {"input-keys": check_input_keys, "output-files": check_output_files, "links": check_links,
              "decision-ids": check_decision_ids, "local-paths": check_local_paths, "gates": check_gates}
    for check in args.checks:
        if check in simple:
            findings[check] = simple[check](repo)
        elif check == "cli" and not args.spodvr:
            notes.append("note: check cli skipped: no --spodvr given")
        elif check == "cli":
            findings[check] = check_cli(repo, args.spodvr)
        elif check == "env-vars":
            findings[check], env_report = check_env_vars(repo)
            report += env_report
        elif check == "numbers":
            numbers = check_numbers(repo)
            report += numbers
            notes.append(f"note: {len(numbers)} numbers repeated across living docs (report-only)")
    return findings, report, notes


def main(argv: list[str]) -> int:
    args = parse_arguments(argv)
    start = time.monotonic()
    repo = Repo(args.root)
    try:
        findings, report, notes = run_checks(repo, args)
    except (OSError, RuntimeError, subprocess.TimeoutExpired) as error:
        print(f"error: check cli could not run the program: {error}; build spodvr or pass another --spodvr")
        return 1
    label = (args.allowlist.resolve().relative_to(args.root).as_posix()
             if args.allowlist.resolve().is_relative_to(args.root) else str(args.allowlist))
    entries, problems = read_allowlist(args.allowlist, label)
    entries = [e for e in entries if e.check in findings]  # rows of skipped checks are not judged
    counts, remaining_all = [], []
    for check, found in findings.items():
        remaining, accepted, _ = apply_allowlist(found, entries, label)
        remaining_all += remaining
        counts.append(f"{check} {len(remaining)}" + (f" (+{accepted} allowed)" if accepted else ""))
    _, _, stale = apply_allowlist([f for found in findings.values() for f in found], entries, label)
    remaining_all += problems + stale
    if report and not args.no_report:
        print("Report-only (never fails):")
        print("\n".join(f"  {line}" for line in report))
    for note in notes:
        print(note)
    for finding in remaining_all:
        print(finding)
    print(f"check_docs: {len(remaining_all)} finding(s): {', '.join(counts)}, allowlist {len(problems + stale)}; "
          f"report-only {len(report)}; {time.monotonic() - start:.1f} s")
    return 1 if remaining_all else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

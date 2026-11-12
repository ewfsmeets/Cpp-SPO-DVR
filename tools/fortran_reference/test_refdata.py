#!/usr/bin/env python3
"""Tests of refdata.py compare with expected/r1713.tsv (ctest FortranReference.Compare*).

Each case builds an OLD tree in the 5.19 style and a NEW tree in the 5.24 style under WORK, at
the paths of real reference cases (small synthetic files), runs compare as a program and checks
its exit status and the FAIL line it prints:
  expected-file        the repository against itself: r1713.tsv is valid, every row matches a
                       file (a removed row may match none: its files are gone after the switch)
  allowed              every class of r1713.tsv at work (labels, a sign change of Im S and of the
                       phase with the one-digit rounding, BK coefficients within the limit, tiny
                       stdout probabilities, unread, removed), with complete MANIFESTs: passes
  edited-number        one number of a stdout of class labels edited: fails
  deleted-file         a reference file missing from NEW: fails
  unknown-label        a stdout label changed that is not on the whitelist: fails
  s-magnitude          a fine fort.25: |S| changed (Re S edited), Re S's sign changed, the sign of
                       the diffraction order in a block header changed: fails (each)
  bk-limit             a BK coefficient changed by more than 1e-12 of the gate-4 scale: fails
  bk-nan               a BK coefficient of NEW NaN, all of them NaN, OLD's scale NaN: fails (each)
  non-finite-number    a printed number of NEW out of the double range (abs class): fails
  number-count         two printed numbers that became one of the same width (abs class): fails
  provenance           NEW's MANIFESTs: a generated file without a line, a file changed after its
                       line, a line of an executable other than BUILDINFO's: fails (each)
  install-mismatch     NEW made by an install other than the expected file's: fails
  fine-without-fort22  a fine-analysis MANIFEST line without its input fort.22: fails
  multi-match          an expected file with a row that overlaps another: refused
  refusals             an expected file with unsorted rows, an unknown class, a limit that does
                       not fit its class (inf included), no install line; a tree without the
                       reference directories: refused (each)
The content cases run with --no-provenance, so that only the classification can fail them.

Usage: python3 -I -B tools/fortran_reference/test_refdata.py CASE WORK
  CASE  one of the cases above;  WORK  scratch directory (recreated)
Exit status: 0 the case behaved as expected, 1 it did not, 2 bad arguments.
Environment: none read. Run time: about a second (expected-file: about 4 s, it reads both trees).
Docs: tools/fortran_reference/README.md "Comparing two reference trees".
"""

import argparse
import hashlib
import re
import shutil
import struct
import subprocess
import sys
from pathlib import Path

TOOLS_DIR = Path(__file__).resolve().parent
REPO_DIR = TOOLS_DIR.parents[1]
REFDATA = TOOLS_DIR / "refdata.py"
EXPECTED = TOOLS_DIR / "expected" / "r1713.tsv"
# The install r1713.tsv names; the stub install of these tests takes its id.
INSTALL_ID = re.search(r"^# install: (\S+)$", EXPECTED.read_text(encoding="utf-8"), re.MULTILINE).group(1)
CASES = ("expected-file", "allowed", "edited-number", "deleted-file", "unknown-label", "s-magnitude", "bk-limit",
         "bk-nan", "non-finite-number", "number-count", "provenance", "install-mismatch", "fine-without-fort22",
         "multi-match", "refusals")
EXPECTED_HEADER = "# install: x\nglob\tallowed\tlimit\treason\n"

HARDWALL = "tests/reference/hardwall_small"
DAILIGHT_SHORT = "tests/reference/cu111_dailight_short"
DAILIGHT_SMALL = "tests/reference/cu111_dailight_small"
KERNELS = "tests/crplib/reference/kernels"

STDOUT_519 = """\
*** SPO-DVR program version 5.19 (Compiled using GFORTRAN for X86 on Oct 26 2026 at 10:27:24) run on 20261026 at 1507 ***
* NAG calls disabled ...
* Using a maximum number of 16 threads in this run.
* The surface lattice constant in X (bohr):                 0.62396600E+01
* The surface lattice constant in Y (bohr):                 0.52396600E+01
* Number of time steps taken:                               7000
     optical potential in z scat.grid:
     nzchsp =                                         137
     z(nzchsp)=                                        8.880000
 iz =   1
 iz =   2
 Total time for potential initialisation =  0.36015001E-04 secs
* Optical potential in Z absorbs with coefficient:          0.19938827E-01
     probability of channel 3                         {probability}
   0  0.9440341E-25
 CPU TIMES
 Total computation                         127.8329 (a)
"""
STDOUT_524 = """\
*** SPO-DVR program version 5.24 (Compiled using GNU for X86 on Nov 11 2026 at 18:16:07) run on 20261111 at 1857 ***
* Using a maximum number of 4 threads in this run.
* The surface lattice constant in U (bohr):                 0.62396600E+01
* The surface lattice constant in V (bohr):                 0.52396600E+01
* Number of time steps taken:                               7000
     start optical potential in z scat.grid:
     nzchmi =                                         136
     z(nzchmi)=                                        8.800000
 iz =   2
 iz =   1
 Total time for potential initialisation =  0.24450000E-04 secs
* Optical potential in Z [high end] absorbs with coefficient:                   0.19938827E-01
     probability of channel 3                         {probability}
   0  0.9440341E-25
 CPU TIMES
 Total computation                         151.6756 (a)
"""
# A fine-analysis fort.25 before and after the r1713 fix: Im S (number 4) and the phase (6) change
# sign, and the second phase also rounds differently in its last digit (as measured on 5.24).
FORT25_519 = """\
 v j mj dif.order are:   0   0   0   1
  m1  m2        s-matrix element          prob          phase
   0   0   0.926238D+00 -0.381233D+00  0.100325D+01  -0.390460
   4   3   0.314358D-06 -0.247390D-06  0.160023D-12  -0.666745
"""
FORT25_524 = """\
 v j mj dif.order are:   0   0   0   1
  m1  m2        s-matrix element          prob          phase
   0   0   0.926238D+00  0.381233D+00  0.100325D+01   0.390460
   4   3   0.314358D-06  0.247390D-06  0.160023D-12   0.666744
"""
# The BK coefficients (re, im) of two records of two coefficients each; max|c| = 0.5.
BK_VALUES = ((0.5, 0.0), (0.1, -0.2), (0.3, 0.1), (0.0, 0.25))
NAN = float("nan")


def bk_file(values=BK_VALUES):
    """A fort.21 of 10 header records and 2 BK records of the given coefficients (re, im)."""
    def record(payload):
        return struct.pack("<i", len(payload)) + payload + struct.pack("<i", len(payload))
    header = b"".join(record(struct.pack("<2i", k, 2)) for k in range(10))
    return header + b"".join(record(struct.pack("<4d", *values[k], *values[k + 1])) for k in (0, 2))


def changed_bk_file(index, re_part, im_part):
    """bk_file with coefficient index replaced."""
    values = list(BK_VALUES)
    values[index] = (re_part, im_part)
    return bk_file(values)


def old_files():
    """The OLD tree (5.19 style): path -> content."""
    return {
        f"{HARDWALL}/DVR.inp": "nStep : 7000;\n",
        f"{HARDWALL}/fortran_stdout.txt": STDOUT_519.format(probability="0.12345678E-03"),
        f"{HARDWALL}/fortran_averages.out": "     0.10000000E+02  0.99999986E+00\n",
        f"{HARDWALL}/fortran_fort25.txt": " 0.5 -0.1\n",
        f"{DAILIGHT_SHORT}/fortran_stdout.txt": STDOUT_519.format(probability="0.10000000E-20"),
        f"{DAILIGHT_SMALL}/fortran_fort21.bin": bk_file(),
        f"{DAILIGHT_SMALL}/fine_full/fort.25": FORT25_519,
        f"{KERNELS}/VERSIONS": "driver 5.19\n",
        f"{KERNELS}/spline.in": "1 2 3\n",
        f"{KERNELS}/spline.out": "cx 1 BFF0F7E3C95C02B0  -1.06051996856710318E+000\n",
    }


def new_files():
    """The NEW tree (5.24 style, every difference allowed by r1713.tsv): path -> content."""
    files = old_files()
    files.update({
        f"{HARDWALL}/fortran_stdout.txt": STDOUT_524.format(probability="0.12345678E-03"),
        f"{HARDWALL}/fortran_fort25.txt": " 0.5 -0.135\n",
        f"{DAILIGHT_SHORT}/fortran_stdout.txt": STDOUT_524.format(probability="0.11000000E-20"),
        f"{DAILIGHT_SMALL}/fortran_fort21.bin": changed_bk_file(1, 0.1, -0.2 + 1e-13),
        f"{DAILIGHT_SMALL}/fine_full/fort.25": FORT25_524,
    })
    del files[f"{KERNELS}/VERSIONS"]
    return files


def write_tree(root, files):
    """Writes path -> content (str or bytes) below root."""
    for path, content in files.items():
        target = Path(root) / path
        target.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(content, bytes):
            target.write_bytes(content)
        else:
            target.write_text(content, encoding="utf-8")


def run(*arguments):
    """(exit status, stdout + stderr) of refdata.py with arguments."""
    result = subprocess.run([sys.executable, "-I", "-B", str(REFDATA), *map(str, arguments)],
                            capture_output=True, text=True)
    return result.returncode, result.stdout + result.stderr


def make_stub_install(work, install_id):
    """A stub install (bin/spodvr.x, BUILDINFO.tsv with install_id) and its directory."""
    install = work / "install"
    (install / "bin").mkdir(parents=True)
    (install / "bin" / "spodvr.x").write_text("stub executable\n")
    sha = hashlib.sha256(b"stub executable\n").hexdigest()
    (install / "BUILDINFO.tsv").write_text(f"key\tvalue\nexe.sha256\t{sha}\nid\t{install_id}\n")
    return install


def record_manifests(work, new, install_id, fine_input=True):
    """MANIFEST.tsv of every case of NEW, from a stub install whose BUILDINFO copy NEW holds; the
    fine-analysis line names its input fort.22 if fine_input."""
    install = make_stub_install(work, install_id)
    copies = new / "tools" / "fortran_reference" / "buildinfo"
    copies.mkdir(parents=True)
    shutil.copy(install / "BUILDINFO.tsv", copies / f"{install_id}.tsv")
    fort22 = work / "fort.22"
    shutil.copy(new / DAILIGHT_SMALL / "fortran_fort21.bin", fort22)
    calls = [(HARDWALL, ["fortran_stdout.txt", "fortran_averages.out", "fortran_fort25.txt"], []),
             (DAILIGHT_SHORT, ["fortran_stdout.txt"], []),
             (DAILIGHT_SMALL, ["fortran_fort21.bin"], []),
             (DAILIGHT_SMALL, ["fine_full/fort.25"], ["--input", fort22] if fine_input else []),
             (KERNELS, ["spline.out"], [])]
    for case, files, inputs in calls:
        status, output = run("manifest", "record", new / case, *files, "--script", "test_refdata.py",
                             "--install", install, "--binary", install / "bin" / "spodvr.x", *inputs)
        if status:
            raise SystemExit(f"FAIL: manifest record for {case}: {output}")


def make_trees(work, install_id=None, fine_input=True):
    """(OLD, NEW) under work, recreated; with install_id, NEW gets MANIFESTs of that install."""
    shutil.rmtree(work, ignore_errors=True)
    old, new = work / "old", work / "new"
    write_tree(old, old_files())
    write_tree(new, new_files())
    if install_id:
        record_manifests(work, new, install_id, fine_input)
    return old, new


def compare(old, new, *options, expected=EXPECTED):
    """(exit status, output) of compare, with r1713.tsv by default."""
    return run("compare", old, new, "--expected", expected, *options)


def expect_pass(status, output):
    """compare exited 0 with no FAIL line."""
    if status != 0 or "FAIL" in output:
        raise SystemExit(f"FAIL: compare should pass (exit {status}):\n{output}")


def expect_failure(status, output, pattern):
    """compare exited 1 with a FAIL line matching pattern."""
    lines = [line for line in output.splitlines() if line.lstrip().startswith("FAIL")]
    if status != 1 or not any(re.search(pattern, line) for line in lines) or "error: " not in output:
        raise SystemExit(f"FAIL: compare should fail with a FAIL line matching {pattern!r} (exit {status}):\n{output}")
    print(f"failed as expected: {next(line for line in lines if re.search(pattern, line)).strip()[:200]}")


def expect_refusal(status, output, pattern):
    """refdata.py exited 1 with "error: <what>; <fix>" matching pattern."""
    if status != 1 or not re.search(rf"^error: .*{pattern}.*; [a-z]", output, re.MULTILINE):
        raise SystemExit(f"FAIL: should refuse with {pattern!r} (exit {status}):\n{output}")
    print(f"refused as expected: {output.strip().splitlines()[-1][:200]}")


def edit(root, path, old, new):
    """Replaces the one occurrence of old in the file root/path by new."""
    target = Path(root) / path
    text = target.read_text(encoding="utf-8")
    if text.count(old) != 1:
        raise SystemExit(f"FAIL: test set-up: {old!r} is not once in {target}")
    target.write_text(text.replace(old, new), encoding="utf-8")


def stdout_cases(name, old, new):
    """The cases that edit a stdout of NEW (classes labels and abs)."""
    if name == "edited-number":
        edit(new, f"{HARDWALL}/fortran_stdout.txt", "0.12345678E-03", "0.12345679E-03")
        expect_failure(*compare(old, new, "--no-provenance"), r"labels .*hardwall_small/fortran_stdout.txt.*0\.12345678E-03 became")
    elif name == "unknown-label":
        edit(new, f"{HARDWALL}/fortran_stdout.txt", "Number of time steps taken", "Number of steps taken")
        expect_failure(*compare(old, new, "--no-provenance"), r"hardwall_small/fortran_stdout.txt.*text .*Number of time steps")
    elif name == "non-finite-number":
        edit(new, f"{DAILIGHT_SHORT}/fortran_stdout.txt", "0.11000000E-20", "0.1100000E+999")
        expect_failure(*compare(old, new, "--no-provenance"), r"abs .*cu111_dailight_short/fortran_stdout.txt.*not finite")
    elif name == "number-count":
        edit(new, f"{DAILIGHT_SHORT}/fortran_stdout.txt", "   0  0.9440341E-25", "   0.00000000000000")
        expect_failure(*compare(old, new, "--no-provenance"), r"abs .*cu111_dailight_short/fortran_stdout.txt.*0\.9440341E-25")


def fort25_cases(old, new):
    """s-magnitude: three edits of the fine fort.25 of NEW, each refused by signs_only."""
    path = f"{DAILIGHT_SMALL}/fine_full/fort.25"
    pattern = r"signs_only .*fine_full/fort.25.*"
    edit(new, path, "0.926238D+00", "0.926239D+00")
    expect_failure(*compare(old, new, "--no-provenance"), pattern + r"number 3 0.926238D\+00 became")
    edit(new, path, " 0.926239D+00", "-0.926238D+00")
    expect_failure(*compare(old, new, "--no-provenance"), pattern + r"number 3 0.926238D\+00 became -")
    write_tree(new, {path: FORT25_524.replace("   0   0   0   1", "   0   0   0  -1")})
    expect_failure(*compare(old, new, "--no-provenance"), pattern + r"line 1: number 4 1 became -1")


def bk_nan_case(work):
    """A NaN coefficient in NEW (one, or all), and a NaN in OLD (the scale), each fail."""
    old, new = make_trees(work)
    target = new / DAILIGHT_SMALL / "fortran_fort21.bin"
    pattern = r"bk_scaled .*fortran_fort21.bin.*"
    target.write_bytes(changed_bk_file(1, 0.1, NAN))
    expect_failure(*compare(old, new, "--no-provenance"), pattern + "not finite")
    target.write_bytes(bk_file(((NAN, NAN),) * 4))
    expect_failure(*compare(old, new, "--no-provenance"), pattern + "not finite")
    (old / DAILIGHT_SMALL / "fortran_fort21.bin").write_bytes(changed_bk_file(0, NAN, 0.0))
    target.write_bytes(changed_bk_file(0, NAN, 1e-13))
    expect_failure(*compare(old, new, "--no-provenance"), pattern + "gate-4 scale is undefined")


def content_case(name, work):
    """The cases that edit NEW's content and compare without provenance."""
    old, new = make_trees(work)
    if name == "deleted-file":
        (new / HARDWALL / "fortran_averages.out").unlink()
        expect_failure(*compare(old, new, "--no-provenance"), r"hardwall_small/fortran_averages.out: missing in NEW")
    elif name == "s-magnitude":
        fort25_cases(old, new)
    elif name == "bk-limit":
        (new / DAILIGHT_SMALL / "fortran_fort21.bin").write_bytes(changed_bk_file(1, 0.1, -0.2 + 1e-11))
        expect_failure(*compare(old, new, "--no-provenance"), r"bk_scaled .*fortran_fort21.bin.*exceeds")
    else:
        stdout_cases(name, old, new)


def provenance_case(work):
    """NEW's MANIFESTs: three faults, each on a fresh copy of the trees."""
    old, new = make_trees(work, INSTALL_ID)
    expect_pass(*compare(old, new))
    (new / HARDWALL / "stray.txt").write_text("not recorded\n")
    expect_failure(*compare(old, new), r"provenance .*hardwall_small/stray.txt: generated, but not in the MANIFEST")
    old, new = make_trees(work, INSTALL_ID)
    (new / HARDWALL / "fortran_averages.out").write_text("     0.10000000E+02  0.99999986E+00\n\n")
    expect_failure(*compare(old, new), r"provenance .*fortran_averages.out: missing, or its sha256 differs")
    old, new = make_trees(work, INSTALL_ID)
    copy = new / "tools" / "fortran_reference" / "buildinfo" / f"{INSTALL_ID}.tsv"
    copy.write_text(re.sub(r"exe\.sha256\t\w+", "exe.sha256\t" + "0" * 64, copy.read_text()))
    expect_failure(*compare(old, new), rf"provenance .*spodvr.x:\w+ is not the executable of {INSTALL_ID}")


def refusals_case(work):
    """The refusals of compare before it compares anything."""
    old, new = make_trees(work)
    expected = work / "expected.tsv"
    bad_files = {
        "not sorted": EXPECTED_HEADER + "tests/b\tlabels\t-\tr\ntests/a\tlabels\t-\tr\n",
        "is not one of": EXPECTED_HEADER + "tests/a\tsimilar\t-\tr\n",
        "does not fit class abs": EXPECTED_HEADER + "tests/a\tabs\t-\tr\n",
        "does not fit class bk_scaled": EXPECTED_HEADER + "tests/a\tbk_scaled\tinf\tr\n",
        "does not fit class signs_only": EXPECTED_HEADER + "tests/a\tsigns_only\t0,4\tr\n",
        "does not fit class labels": EXPECTED_HEADER + "tests/a\tlabels\t1e-12\tr\n",
        "has 0 '# install: <id>' lines": "glob\tallowed\tlimit\treason\ntests/a\tlabels\t-\tr\n",
    }
    for pattern, text in bad_files.items():
        expected.write_text(text, encoding="utf-8")
        expect_refusal(*compare(old, new, expected=expected), pattern)
    shutil.rmtree(new / "tests" / "crplib")
    expect_refusal(*compare(old, new, "--no-provenance"), "has no tests/reference and tests/crplib/reference")


def multi_match_case(work):
    """r1713.tsv plus a row whose glob overlaps another one's: refused."""
    old, new = make_trees(work)
    lines = EXPECTED.read_text(encoding="utf-8").splitlines()
    comments = [line for line in lines if line.startswith("#")]
    rows = [line for line in lines if not line.startswith("#")]
    overlap = "tests/reference/*_small/fortran_stdout.txt\tabs\t1e-20\toverlaps the hardwall_* row"
    expected = work / "expected.tsv"
    expected.write_text("\n".join(comments + rows[:1] + sorted(rows[1:] + [overlap])) + "\n", encoding="utf-8")
    expect_refusal(*compare(old, new, "--no-provenance", expected=expected),
                   r"hardwall_small/fortran_stdout.txt matches the rows")


def run_case(name, work):
    """Runs one case; raises SystemExit with a FAIL message if it does not behave as expected."""
    if name == "expected-file":
        status, output = compare(REPO_DIR, REPO_DIR, "--no-provenance")
        expect_pass(status, output)
        # A glob typo; a removed row matches nothing once the switch has deleted its files.
        if re.search(r"note: .*\((?!removed\))", output):
            raise SystemExit(f"FAIL: a row of {EXPECTED} matches no committed file:\n{output}")
    elif name == "allowed":
        old, new = make_trees(work, INSTALL_ID)
        status, output = compare(old, new)
        expect_pass(status, output)
        for allowed in ("labels", "abs", "signs_only", "bk_scaled", "unread", "removed"):
            if not re.search(rf"^  ok    {allowed} ", output, re.MULTILINE):
                raise SystemExit(f"FAIL: no file of class {allowed} was compared:\n{output}")
    elif name == "provenance":
        provenance_case(work)
    elif name == "install-mismatch":
        expect_failure(*compare(*make_trees(work, "other-install")), rf"MANIFEST lines name install other-install, not {INSTALL_ID}")
    elif name == "fine-without-fort22":
        expect_failure(*compare(*make_trees(work, INSTALL_ID, fine_input=False)), r"fine_full/fort.25: a fine-analysis file without input fort.22")
    elif name == "bk-nan":
        bk_nan_case(work)
    elif name == "multi-match":
        multi_match_case(work)
    elif name == "refusals":
        refusals_case(work)
    else:
        content_case(name, work)


def main(argv):
    """Runs the named case in WORK."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case", metavar="CASE", choices=CASES, help="the case to run")
    parser.add_argument("work", metavar="WORK", type=Path, help="scratch directory (recreated)")
    args = parser.parse_args(argv)
    run_case(args.case, args.work.resolve())
    print(f"PASS: {args.case}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

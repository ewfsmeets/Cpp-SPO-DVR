"""Tests of refdata.py compare with expected/r1713.tsv (ctest FortranReference.Compare*), and of the
copies of refdata.py's lists (ctest FortranReference.ReferenceListsAgree).

Most tests build an OLD tree in the 5.19 style and a NEW tree in the 5.24 style (test_initial_phase_*
and test_nproj_label: r1713 and r1727 style) in their scratch directory `work`, at the paths of real
reference cases (small synthetic files), run compare as a program and check its exit status and the
FAIL line it prints:
  test_expected_file         the repository against itself: r1727.tsv is valid, every row matches a
                             file (a removed row may match none: its files are gone after the
                             switch)
  test_allowed               every class of r1713.tsv at work (labels, a sign change of Im S and of
                             the phase with the one-digit rounding, BK coefficients within the
                             limit, tiny stdout probabilities, unread, removed, added), with complete
                             MANIFESTs: passes
  test_added_file            a file of class added without its MANIFEST line, a new file no added
                             row matches (even with a MANIFEST line), a file of class added present
                             in both trees but changed: fails (each)
  test_edited_number         one number of a stdout of class labels edited: fails
  test_deleted_file          a reference file missing from NEW: fails
  test_unknown_label         a stdout label changed that is not on the whitelist: fails
  test_s_magnitude           a fine fort.25: |S| changed (Re S edited), Re S's sign changed, the
                             sign of the diffraction order in a block header changed: fails (each)
  test_bk_limit              a BK coefficient changed by more than 1e-12 of the gate-4 scale: fails
  test_bk_nan                a BK coefficient of NEW NaN, all of them NaN, OLD's scale NaN: fails
                             (each)
  test_non_finite_number     a printed number of NEW out of the double range (abs class): fails
  test_number_count          two printed numbers that became one of the same width (abs class):
                             fails
  test_provenance            NEW's MANIFESTs: a generated file without a line, a file changed after
                             its line, a line of an executable other than BUILDINFO's: fails (each)
  test_install_mismatch      NEW made by an install other than the expected file's: fails
  test_fine_without_fort22   a fine-analysis MANIFEST line without its input fort.22: fails
  test_multi_match           an expected file with a row that overlaps another: refused
  test_refusals              an expected file with unsorted rows, an unknown class, a limit that
                             does not fit its class (inf included), a row without reason, no install
                             line; a tree without the reference directories: refused (each)
  test_initial_phase_*       class initial_phase (an on-the-fly fort.25 of a case whose specular line is
                             fixed): the initial state's rows rotated by -k_v (Z_sp - Zana) for v = 0 and 1,
                             |S| and P kept: passes; a changed row of another state, a changed |S|, a wrong
                             rotation, P changed by 2 printed units, a changed block header: fails (each)
  test_nproj_label           the stdout label "(nproj = 2)" that r1725 prints as "(nproj = 1)": passes under
                             labels; with its number changed: fails
  test_reference_lists       the MANIFEST header and column count, the reference trees and
                             NOT_A_CASE in ReferenceManifest.cmake, tools/check_docs.py and
                             reference_config.sh equal refdata.py's, and the exclude of .pre-commit-
                             config.yaml covers every tree; copies with a tree or NOT_A_CASE changed
                             are reported (5 lines)
The tests that edit NEW's content run compare with --no-provenance, so that only the classification can fail
them.

Usage: uv run python -I -B -m pytest -c tools/fortran_reference/pytest.ini \
           --rootdir=tools/fortran_reference tools/fortran_reference/test_refdata.py
Environment: none read. Run time: about a second (test_expected_file: about 4 s, it reads both trees).
Docs: tools/fortran_reference/README.md "Comparing two reference trees".
"""

import cmath
import hashlib
import importlib.util
import math
import re
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest

TOOLS_DIR = Path(__file__).resolve().parent
REPO_DIR = TOOLS_DIR.parents[1]
REFDATA = TOOLS_DIR / "refdata.py"
EXPECTED = TOOLS_DIR / "expected" / "r1713.tsv"  # the synthetic 5.19 -> 5.24 trees of test_allowed
CURRENT_EXPECTED = TOOLS_DIR / "expected" / "r1727.tsv"  # the switch that made the committed references
# The install r1713.tsv names; the stub install of these tests takes its id.
INSTALL_ID = re.search(r"^# install: (\S+)$", EXPECTED.read_text(encoding="utf-8"), re.MULTILINE).group(1)
EXPECTED_HEADER = "# install: x\nglob\tallowed\tlimit\treason\n"

HARDWALL = "tests/reference/hardwall_small"
DAILIGHT_SHORT = "tests/reference/cu111_dailight_short"
DAILIGHT_SMALL = "tests/reference/cu111_dailight_small"
KERNELS = "tests/crplib/reference/kernels"
SCM_EXAMPLE = "tests/scm/reference/example"

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
# A write2DPES cut (pes_cut_<k>.dat, gate 1c): new in NEW, class added.
PES_CUT = "#u=         0.0000000000  v=         0.0000000000 theta=         1.5707963268 phi=         0.0000000000\n"
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
        f"{SCM_EXAMPLE}/points.in": "0.0 0.0 6.0 1.4 1.5707963267948966 0.0\n",
        f"{SCM_EXAMPLE}/rescaling_on.out": "SCM_DRIVER_DONE\n",
    }


def new_files():
    """The NEW tree (5.24 style, every difference allowed by r1713.tsv): path -> content."""
    files = old_files()
    files.update(
        {
            f"{HARDWALL}/fortran_stdout.txt": STDOUT_524.format(probability="0.12345678E-03"),
            f"{HARDWALL}/fortran_fort25.txt": " 0.5 -0.135\n",
            f"{DAILIGHT_SHORT}/fortran_stdout.txt": STDOUT_524.format(probability="0.11000000E-20"),
            f"{DAILIGHT_SMALL}/fortran_fort21.bin": changed_bk_file(1, 0.1, -0.2 + 1e-13),
            f"{DAILIGHT_SMALL}/fine_full/fort.25": FORT25_524,
            f"{HARDWALL}/pes_cut_1.dat": PES_CUT,
            f"{HARDWALL}/pes_cuts.json": '{"cuts": [{}]}\n',
        }
    )
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
    result = subprocess.run(
        [sys.executable, "-I", "-B", str(REFDATA), *map(str, arguments)], capture_output=True, text=True
    )
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
    calls = [
        (HARDWALL, ["fortran_stdout.txt", "fortran_averages.out", "fortran_fort25.txt", "pes_cut_1.dat"], []),
        (DAILIGHT_SHORT, ["fortran_stdout.txt"], []),
        (DAILIGHT_SMALL, ["fortran_fort21.bin"], []),
        (DAILIGHT_SMALL, ["fine_full/fort.25"], ["--input", fort22] if fine_input else []),
        (KERNELS, ["spline.out"], []),
        (SCM_EXAMPLE, ["rescaling_on.out"], []),
    ]
    for case, files, inputs in calls:
        status, output = run(
            "manifest",
            "record",
            new / case,
            *files,
            "--script",
            "test_refdata.py",
            "--install",
            install,
            "--binary",
            install / "bin" / "spodvr.x",
            *inputs,
        )
        if status:
            pytest.fail(f"manifest record for {case}: {output}")


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
        pytest.fail(f"compare should pass (exit {status}):\n{output}")


def expect_failure(status, output, pattern):
    """compare exited 1 with a FAIL line matching pattern."""
    lines = [line for line in output.splitlines() if line.lstrip().startswith("FAIL")]
    if status != 1 or not any(re.search(pattern, line) for line in lines) or "error: " not in output:
        pytest.fail(f"compare should fail with a FAIL line matching {pattern!r} (exit {status}):\n{output}")
    print(f"failed as expected: {next(line for line in lines if re.search(pattern, line)).strip()[:200]}")


def expect_refusal(status, output, pattern):
    """refdata.py exited 1 with "error: <what>; <fix>" matching pattern."""
    if status != 1 or not re.search(rf"^error: .*{pattern}.*; [a-z]", output, re.MULTILINE):
        pytest.fail(f"should refuse with {pattern!r} (exit {status}):\n{output}")
    print(f"refused as expected: {output.strip().splitlines()[-1][:200]}")


def edit(root, path, old, new):
    """Replaces the one occurrence of old in the file root/path by new."""
    target = Path(root) / path
    text = target.read_text(encoding="utf-8")
    if text.count(old) != 1:
        pytest.fail(f"test set-up: {old!r} is not once in {target}")
    target.write_text(text.replace(old, new), encoding="utf-8")


# Class initial_phase: a case of DAILIGHT_SHORT with VIn = JIn = MJIn = 0, Zana 8, the specular line fixed at
# Z_sp = 10.88, the gas-phase levels E(0, 0) and E(1, 0), and one energy block with the initial state's rows for
# v = 0 and 1 and a row of another state.
INITIAL_PHASE_DVR = "VIn : INTEGER : 0;\nJIn : INTEGER : 0;\nMJIn : INTEGER : 0;\nEParallel : 0.0;\n"
INITIAL_PHASE_STDOUT = """\
* Mass of atom 1 (amu):                                     0.10078000E+01
* Mass of atom 2 (amu):                                     0.10078000E+01
* Wave function analysed at (bohr):                         0.80000000E+01
* Tolerance for transferring the wavefunction (nproj = {nproj}):  0.10000000E-06
 asymptotic molecular vibrational wave functions

   v   j       e(ev)    erel(ev)         psi(1)        psi(nr)       norm

   0   0   -4.467141    0.000000  0.1853026D-01  0.1513393D-02  1.0000000 1.4425817

   1   0   -3.963306    0.503835  0.1793600D-01  0.1500000D-02  1.0000000 1.4800000

   ==> ZAnalysis on specular channel fixed at Z =  0.10880000E+02
"""
# The initial rows (v, S) and a row of another state (m1, m2) = (1, 0) of v = 0, as OLD prints them.
INITIAL_PHASE_ROWS = ((0, complex(-0.470660e-03, -0.850680e-03)), (1, complex(0.212345e-02, -0.123456e-02)))
OTHER_ROW = complex(0.100000e-02, 0.200000e-02)
# eV; the six fixed-line cases of the r1727 switch have EParallel 0 (class initial_phase refuses EParallel /= 0),
# so the block's momentum is k_z of v = VIn
TRANS_ENERGY = 0.7
MASS = 1822.89 * (1.0078 + 1.0078)
CHAEV = 27.2114


def fortran_d(x):
    """x in the Fortran's D13.6 format (0.dddddd and the exponent)."""
    if x == 0.0:
        return "0.000000D+00"
    exponent = math.floor(math.log10(abs(x))) + 1
    mantissa = round(x / 10.0**exponent, 6)
    if abs(mantissa) >= 1.0:
        mantissa, exponent = mantissa / 10.0, exponent + 1
    return f"{mantissa:.6f}D{exponent:+03d}"


def from_fortran(token):
    """The value of a printed D number."""
    return float(token.replace("D", "E"))


def fort25_row(m1, m2, s):
    """A row m1 m2 Re S Im S P phase as smatrix_combine prints it, from the printed Re S and Im S."""
    re_s, im_s = fortran_d(s.real), fortran_d(s.imag)
    printed = complex(from_fortran(re_s), from_fortran(im_s))
    phase = math.atan2(printed.imag, printed.real) if printed.real != 0.0 else 0.0
    return f"{m1:4d}{m2:4d}  {re_s:>13} {im_s:>13} {fortran_d(abs(printed) ** 2):>13} {phase:10.6f}"


def initial_phase_fort25(rotation=lambda v: 1.0, other=OTHER_ROW, order=0, j2_row=None):
    """A fort.25 block: each initial row multiplied by rotation(v), the other row as given; the header of v = 0
    with diffraction order `order`; with j2_row, also a block of the state (v, j, mJ, order) = (0, 2, 0, 0) with
    that S as its row (0, 0)."""
    k = math.sqrt(2.0 * MASS * TRANS_ENERGY / CHAEV)
    lines = [
        " total energy = " + fortran_d(-4.467141 + TRANS_ENERGY).replace("D", "D", 1),
        " trans energy = " + fortran_d(TRANS_ENERGY),
        " trans energy in z = " + fortran_d(TRANS_ENERGY),
        f" total momentum =  {k / 100.0:.8f}D+02",
    ]
    for v, s in INITIAL_PHASE_ROWS:
        lines.append(f" v j mj dif.order are:   {v}   0   0   {order if v == 0 else 0}")
        lines.append("  m1  m2        s-matrix element          prob          phase")
        lines.append(fort25_row(0, 0, s * rotation(v)))
        if v == 0:
            lines.append(fort25_row(1, 0, other))
    if j2_row is not None:
        lines.append(" v j mj dif.order are:   0   2   0   0")
        lines.append("  m1  m2        s-matrix element          prob          phase")
        lines.append(fort25_row(0, 0, j2_row))
    return "\n".join(lines) + "\n"


def momentum(v):
    """k_v of the initial state's channel at level v, as class initial_phase computes it."""
    k_in = math.sqrt(2.0 * MASS * TRANS_ENERGY / CHAEV)
    return k_in if v == 0 else math.sqrt(k_in * k_in + 2.0 * MASS * (-4.467141 - -3.963306) / CHAEV)


def r1720_rotation(v, extra=0.0):
    """exp(-i k_v (Z_sp - Zana)) times exp(i extra)."""
    return cmath.exp(1j * (-momentum(v) * (10.88 - 8.0) + extra))


def initial_phase_trees(work, new_fort25, new_stdout=None, old_fort25=None):
    """(OLD, NEW, expected file) for class initial_phase: the case as above, NEW's fort.25 as given (OLD's too,
    if given)."""
    shutil.rmtree(work, ignore_errors=True)
    old, new = work / "old", work / "new"
    for root, fort25, stdout in (
        (old, old_fort25 or initial_phase_fort25(), INITIAL_PHASE_STDOUT.format(nproj=2)),
        (new, new_fort25, new_stdout or INITIAL_PHASE_STDOUT.format(nproj=2)),
    ):
        for tree in ("tests/reference", "tests/crplib/reference", "tests/scm/reference"):
            (root / tree).mkdir(parents=True)
        write_tree(
            root,
            {
                f"{DAILIGHT_SHORT}/DVR.inp": INITIAL_PHASE_DVR,
                f"{DAILIGHT_SHORT}/fortran_stdout.txt": stdout,
                f"{DAILIGHT_SHORT}/fortran_fort25.txt": fort25,
            },
        )
    expected = work / "expected.tsv"
    expected.write_text(
        EXPECTED_HEADER
        + f"{DAILIGHT_SHORT}/fortran_fort25.txt\tinitial_phase\t-\tr1720\n"
        + f"{DAILIGHT_SHORT}/fortran_stdout.txt\tlabels\t-\tr1725\n",
        encoding="utf-8",
    )
    return old, new, expected


def compare_initial_phase(work, new_fort25, new_stdout=None, old_fort25=None):
    """(exit status, output) of compare on initial_phase_trees."""
    old, new, expected = initial_phase_trees(work, new_fort25, new_stdout, old_fort25)
    return compare(old, new, "--no-provenance", expected=expected)


def test_initial_phase_rotation_passes(work):
    """The initial rows of v = 0 and 1 rotated by -k_v (Z_sp - Zana), the other row unchanged: passes."""
    status, output = compare_initial_phase(work, initial_phase_fort25(rotation=r1720_rotation))
    expect_pass(status, output)
    if not re.search(r"^  ok    initial_phase .*fortran_fort25.txt: 2 lines differ", output, re.MULTILINE):
        pytest.fail(f"the two initial rows should differ and be allowed:\n{output}")


def test_initial_phase_other_row_fails(work):
    """A row of another state changed (as if rotated too): fails."""
    fort25 = initial_phase_fort25(rotation=r1720_rotation, other=OTHER_ROW * r1720_rotation(0))
    expect_failure(*compare_initial_phase(work, fort25), r"initial_phase .*not a row of the initial state")


def test_initial_phase_s_magnitude_fails(work):
    """The initial row of v = 0 rotated and scaled by 1.001: fails on |S|."""
    fort25 = initial_phase_fort25(rotation=lambda v: r1720_rotation(v) * (1.001 if v == 0 else 1.0))
    expect_failure(*compare_initial_phase(work, fort25), r"initial_phase .*\|S\| changed")


def test_initial_phase_wrong_rotation_fails(work):
    """The initial row of v = 1 rotated 0.01 rad too far: fails on S (and the phase)."""
    fort25 = initial_phase_fort25(rotation=lambda v: r1720_rotation(v, extra=0.01 if v == 1 else 0.0))
    expect_failure(*compare_initial_phase(work, fort25), r"initial_phase .*not S_old exp")


def test_initial_phase_phase_column_fails(work):
    """The rows rotated as r1720 does, but the printed phase of v = 1 0.01 rad off: fails on the phase."""
    fort25 = initial_phase_fort25(rotation=r1720_rotation)
    row = [line for line in fort25.splitlines() if line.startswith("   0   0")][1]
    phase = row.split()[5]
    fort25 = fort25.replace(row, row[: -len(phase)] + f"{float(phase) + 0.01:.6f}")
    expect_failure(*compare_initial_phase(work, fort25), r"initial_phase .*the phase changed by")


def test_initial_phase_p_change_fails(work):
    """The rotated rows with P of v = 0 raised by 2 printed units: fails."""
    fort25 = initial_phase_fort25(rotation=r1720_rotation)
    row = next(line for line in fort25.splitlines() if line.startswith("   0   0"))
    p = row.split()[4]
    raised = fortran_d(from_fortran(p) + 2e-12)
    fort25 = fort25.replace(row, row.replace(p, raised))
    expect_failure(*compare_initial_phase(work, fort25), r"initial_phase .*P \S+ became \S+, more than one printed")


def test_initial_phase_header_fails(work):
    """The block header of v = 0 with diffraction order 1 instead of 0: fails."""
    fort25 = initial_phase_fort25(rotation=r1720_rotation, order=1)
    expect_failure(*compare_initial_phase(work, fort25), r"initial_phase .*dif.order are.*not a row of the initial")


def test_initial_phase_unrotated_row_fails(work):
    """v = 0 rotated, the initial row of v = 1 left as it was (identical): fails, an unchanged initial row must be
    the rotation too."""
    fort25 = initial_phase_fort25(rotation=lambda v: r1720_rotation(v) if v == 0 else 1.0)
    expect_failure(*compare_initial_phase(work, fort25), r"initial_phase .*not S_old exp\(-i k \(Z_sp - Zana\)\)")


def test_initial_phase_row_of_another_state_fails(work):
    """The row (0, 0) of the state j = 2 (not JIn) rotated as r1720 rotates the initial state: fails."""
    other = complex(0.300000e-02, -0.100000e-02)
    fort25 = initial_phase_fort25(rotation=r1720_rotation, j2_row=other * r1720_rotation(0))
    status, output = compare_initial_phase(work, fort25, old_fort25=initial_phase_fort25(j2_row=other))
    expect_failure(status, output, r"initial_phase .*not a row of the initial state")


def test_initial_phase_im_sign_flipped_fails(work):
    """The v = 0 row rotated, but its Im S printed with the wrong sign (|S|, P and the phase column as rotated):
    fails on S."""
    fort25 = initial_phase_fort25(rotation=r1720_rotation)
    row = next(line for line in fort25.splitlines() if line.startswith("   0   0"))
    fields = row.split()
    flipped = fortran_d(-from_fortran(fields[3]))
    fort25 = fort25.replace(row, row.replace(fields[3], flipped))
    expect_failure(*compare_initial_phase(work, fort25), r"initial_phase .*not S_old exp")


def test_three_digit_exponent():
    """gfortran prints a three-digit exponent without its letter (0.100000-119): read with it."""
    refdata = load_module("refdata", REFDATA)
    assert refdata.as_float("0.100000-119") == 0.1e-119
    assert refdata.as_float("-0.250000+123") == -0.25e123
    assert refdata.as_float("0.500000D+00") == 0.5
    assert refdata.printed_unit("0.100000-119") == pytest.approx(1e-125, rel=1e-12)


def test_nproj_label(work):
    """r1725's "(nproj = 1)" label passes under labels; with its number changed it fails."""
    fort25 = initial_phase_fort25()
    expect_pass(*compare_initial_phase(work, fort25, INITIAL_PHASE_STDOUT.format(nproj=1)))
    changed = INITIAL_PHASE_STDOUT.format(nproj=1).replace("0.10000000E-06", "0.20000000E-06")
    expect_failure(*compare_initial_phase(work, fort25, changed), r"labels .*fortran_stdout.txt.*became")


def test_bk_nan(work):
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


def test_provenance(work):
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


def test_refusals(work):
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
        "a row needs a glob and a reason": EXPECTED_HEADER + "tests/a\tlabels\t-\t-\n",
        "has 0 '# install: <id>' lines": "glob\tallowed\tlimit\treason\ntests/a\tlabels\t-\tr\n",
    }
    for pattern, text in bad_files.items():
        expected.write_text(text, encoding="utf-8")
        expect_refusal(*compare(old, new, expected=expected), pattern)
    shutil.rmtree(new / "tests" / "crplib")
    expect_refusal(
        *compare(old, new, "--no-provenance"),
        "lacks one of tests/reference, tests/crplib/reference, tests/scm/reference",
    )


def test_multi_match(work):
    """r1713.tsv plus a row whose glob overlaps another one's: refused."""
    old, new = make_trees(work)
    lines = EXPECTED.read_text(encoding="utf-8").splitlines()
    comments = [line for line in lines if line.startswith("#")]
    rows = [line for line in lines if not line.startswith("#")]
    overlap = "tests/reference/*_small/fortran_stdout.txt\tabs\t1e-20\toverlaps the hardwall_* row"
    expected = work / "expected.tsv"
    expected.write_text("\n".join(comments + rows[:1] + sorted(rows[1:] + [overlap])) + "\n", encoding="utf-8")
    expect_refusal(
        *compare(old, new, "--no-provenance", expected=expected), r"hardwall_small/fortran_stdout.txt matches the rows"
    )


def load_module(name, path):
    """The Python tool at path as a module (python3 -I keeps the tools' directories off sys.path)."""
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def cmake_lists(text):
    """name -> values of the set(name value...) lines of a CMake script; a quoted value is one
    string, with its \\t escapes as tabs."""
    values = {}
    for name, value in re.findall(r"^set\((\w+) (.*)\)$", text, re.MULTILINE):
        values[name] = (value.strip('"').replace("\\t", "\t"),) if value.startswith('"') else tuple(value.split())
    return values


def list_mismatches(cmake_text, check_docs_trees, config_trees, precommit_exclude):
    """One line per copy of a refdata.py list that differs from it: ReferenceManifest.cmake (the text
    of the script), check_docs.py REFERENCE_TREES, reference_config.sh FORTRAN_REFERENCE_TREES; and
    the trees whose files the pre-commit hooks would edit (precommit_exclude: the exclude regex)."""
    refdata = load_module("refdata", REFDATA)
    cmake = cmake_lists(cmake_text)
    edited = tuple(tree for tree in refdata.REFERENCE_TREES if not re.match(precommit_exclude, f"{tree}/MANIFEST.tsv"))
    copies = {
        "ReferenceManifest.cmake manifest_header": (
            cmake.get("manifest_header"),
            ("\t".join(refdata.MANIFEST_COLUMNS),),
        ),
        "ReferenceManifest.cmake manifest_columns": (
            cmake.get("manifest_columns"),
            (str(len(refdata.MANIFEST_COLUMNS)),),
        ),
        "ReferenceManifest.cmake reference_trees": (cmake.get("reference_trees"), refdata.REFERENCE_TREES),
        "ReferenceManifest.cmake not_a_case": (cmake.get("not_a_case"), refdata.NOT_A_CASE),
        "check_docs.py REFERENCE_TREES": (
            tuple(tree.rstrip("/") for tree in check_docs_trees),
            refdata.REFERENCE_TREES,
        ),
        "reference_config.sh FORTRAN_REFERENCE_TREES": (tuple(config_trees), refdata.REFERENCE_TREES),
        "the trees not excluded in .pre-commit-config.yaml": (edited, ()),
    }
    return [f"{name} is {copy}, not {original}" for name, (copy, original) in copies.items() if copy != original]


def test_reference_lists():
    """The copies agree with refdata.py; copies with a tree or NOT_A_CASE changed do not."""
    cmake_text = (TOOLS_DIR / "ReferenceManifest.cmake").read_text(encoding="utf-8")
    check_docs_trees = load_module("check_docs", REPO_DIR / "tools" / "check_docs.py").REFERENCE_TREES
    config_trees = subprocess.run(
        [
            "bash",
            "-c",
            'source "$0"; printf "%s\\n" "${FORTRAN_REFERENCE_TREES[@]}"',
            str(TOOLS_DIR / "reference_config.sh"),
        ],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    precommit = (REPO_DIR / ".pre-commit-config.yaml").read_text(encoding="utf-8")
    exclude = re.search(r"^exclude: (\S+)", precommit, re.MULTILINE).group(1)
    mismatches = list_mismatches(cmake_text, check_docs_trees, config_trees, exclude)
    if mismatches:
        pytest.fail("a copy of a refdata.py list differs:\n" + "\n".join(mismatches))
    planted = cmake_text.replace(" tests/scm/reference)", ")").replace(
        "set(not_a_case fortran_unformatted)", "set(not_a_case other)"
    )
    found = list_mismatches(planted, check_docs_trees[:2], config_trees[:2], exclude.replace("|scm", ""))
    if len(found) != 5:
        pytest.fail(f"5 planted differences, {len(found)} reported:\n" + "\n".join(found))
    print("planted differences reported as expected:\n" + "\n".join(found))


def test_expected_file():
    """The repository against itself: r1727.tsv (the latest switch) is valid and every row matches a committed
    file."""
    status, output = compare(REPO_DIR, REPO_DIR, "--no-provenance", expected=CURRENT_EXPECTED)
    expect_pass(status, output)
    # A glob typo; a removed row matches nothing once the switch has deleted its files.
    if re.search(r"note: .*\((?!removed\))", output):
        pytest.fail(f"a row of {CURRENT_EXPECTED} matches no committed file:\n{output}")


def test_allowed(work):
    """Every class of r1713.tsv at work, with complete MANIFESTs: passes, and each class was used."""
    old, new = make_trees(work, INSTALL_ID)
    status, output = compare(old, new)
    expect_pass(status, output)
    for allowed in ("labels", "abs", "signs_only", "bk_scaled", "unread", "removed", "added"):
        if not re.search(rf"^  ok    {allowed} ", output, re.MULTILINE):
            pytest.fail(f"no file of class {allowed} was compared:\n{output}")


def test_added_file(work):
    """A file of class added (a pes_cut file) without its MANIFEST line fails; a new file that no added row
    matches fails even with its MANIFEST line; a file of class added in both trees fails when it changed."""
    old, new = make_trees(work, INSTALL_ID)
    manifest = new / HARDWALL / "MANIFEST.tsv"
    lines = manifest.read_text(encoding="utf-8").splitlines(keepends=True)
    manifest.write_text("".join(line for line in lines if not line.startswith("pes_cut_1.dat\t")), encoding="utf-8")
    expect_failure(
        *compare(old, new, "--no-provenance"), r"added .*hardwall_small/pes_cut_1.dat: added in NEW without a line"
    )
    old, new = make_trees(work, INSTALL_ID)
    (new / HARDWALL / "pes_cut_1.dat").rename(new / HARDWALL / "pes_slice_1.dat")
    edit(new, f"{HARDWALL}/MANIFEST.tsv", "pes_cut_1.dat\t", "pes_slice_1.dat\t")
    expect_failure(*compare(old, new, "--no-provenance"), r"pes_slice_1.dat: added in NEW; not allowed")
    old, new = make_trees(work, INSTALL_ID)
    write_tree(old, {f"{HARDWALL}/pes_cut_1.dat": PES_CUT.replace("1.5707963268", "1.5707963267")})
    expect_failure(
        *compare(old, new, "--no-provenance"), r"added .*hardwall_small/pes_cut_1.dat: .*not allowed: line 1: number 3"
    )


def test_install_mismatch(work):
    expect_failure(
        *compare(*make_trees(work, "other-install")),
        rf"MANIFEST lines name install other-install, not {INSTALL_ID}",
    )


def test_fine_without_fort22(work):
    expect_failure(
        *compare(*make_trees(work, INSTALL_ID, fine_input=False)),
        r"fine_full/fort.25: a fine-analysis file without input fort.22",
    )


def test_edited_number(work):
    """One number of a stdout of class labels edited: fails."""
    old, new = make_trees(work)
    edit(new, f"{HARDWALL}/fortran_stdout.txt", "0.12345678E-03", "0.12345679E-03")
    expect_failure(
        *compare(old, new, "--no-provenance"), r"labels .*hardwall_small/fortran_stdout.txt.*0\.12345678E-03 became"
    )


def test_deleted_file(work):
    """A reference file missing from NEW: fails."""
    old, new = make_trees(work)
    (new / HARDWALL / "fortran_averages.out").unlink()
    expect_failure(*compare(old, new, "--no-provenance"), r"hardwall_small/fortran_averages.out: missing in NEW")


def test_unknown_label(work):
    """A stdout label changed that is not on the whitelist: fails."""
    old, new = make_trees(work)
    edit(new, f"{HARDWALL}/fortran_stdout.txt", "Number of time steps taken", "Number of steps taken")
    expect_failure(
        *compare(old, new, "--no-provenance"), r"hardwall_small/fortran_stdout.txt.*text .*Number of time steps"
    )


def test_s_magnitude(work):
    """Three edits of the fine fort.25 of NEW, each refused by signs_only: |S| changed, Re S's sign changed, the
    sign of the diffraction order in a block header changed."""
    old, new = make_trees(work)
    path = f"{DAILIGHT_SMALL}/fine_full/fort.25"
    pattern = r"signs_only .*fine_full/fort.25.*"
    edit(new, path, "0.926238D+00", "0.926239D+00")
    expect_failure(*compare(old, new, "--no-provenance"), pattern + r"number 3 0.926238D\+00 became")
    edit(new, path, " 0.926239D+00", "-0.926238D+00")
    expect_failure(*compare(old, new, "--no-provenance"), pattern + r"number 3 0.926238D\+00 became -")
    write_tree(new, {path: FORT25_524.replace("   0   0   0   1", "   0   0   0  -1")})
    expect_failure(*compare(old, new, "--no-provenance"), pattern + r"line 1: number 4 1 became -1")


def test_bk_limit(work):
    """A BK coefficient changed by more than 1e-12 of the gate-4 scale: fails."""
    old, new = make_trees(work)
    (new / DAILIGHT_SMALL / "fortran_fort21.bin").write_bytes(changed_bk_file(1, 0.1, -0.2 + 1e-11))
    expect_failure(*compare(old, new, "--no-provenance"), r"bk_scaled .*fortran_fort21.bin.*exceeds")


def test_non_finite_number(work):
    """A printed number of NEW out of the double range (abs class): fails."""
    old, new = make_trees(work)
    edit(new, f"{DAILIGHT_SHORT}/fortran_stdout.txt", "0.11000000E-20", "0.1100000E+999")
    expect_failure(*compare(old, new, "--no-provenance"), r"abs .*cu111_dailight_short/fortran_stdout.txt.*not finite")


def test_number_count(work):
    """Two printed numbers that became one of the same width (abs class): fails."""
    old, new = make_trees(work)
    edit(new, f"{DAILIGHT_SHORT}/fortran_stdout.txt", "   0  0.9440341E-25", "   0.00000000000000")
    expect_failure(
        *compare(old, new, "--no-provenance"), r"abs .*cu111_dailight_short/fortran_stdout.txt.*0\.9440341E-25"
    )

#!/usr/bin/env bash
# Generates the Fortran reference data of gate 1c (write2DPES: V(r, Z) at a fixed x, y, theta, phi).
#
# For each cut k (from 1) of the case's pes_cuts.json (a list "cuts" of objects whose keys are the input
# keys x2DPES, y2DPES (bohr), theta2DPES, phi2DPES (degrees); a key left out keeps the Fortran default),
# runs the reference executable (reference_config.sh) in a fresh directory with the case's DVR.inp plus
# `write2DPES : LOGICAL : TRUE;` and those keys as `Key : REAL : value;` (the Fortran reads these keys only
# with their type: InputNew.f90 GetRealField). The Fortran writes 2dpes.dat right after the PES is set up,
# before the 6D potential and the gas phase, and stops (Main.f90). After all cuts have run, the script copies
# them into tests/reference/<case>/:
#   pes_cut_<k>.dat  2dpes.dat of cut k: a header '#u= x v= y theta= t phi= p' (x, y in bohr, theta,
#                    phi in radians), then r, Z (bohr) and V (eV; the PES itself: no switching, no
#                    VMin/VMax) as 3(4x,e15.8), Z fastest, a blank line after each r block
#                    (Potential.f90 Write2DPESToFile)
# and their lines in the case's MANIFEST.tsv (inputs: the DVR.inp as run and the staged data files).
# Refuses, before the runs, a bad pes_cuts.json, a case DVR.inp that sets one of these keys itself, a
# pes_cut_<n + 1>.dat left from more cuts, and to overwrite a file of these that differs from git
# (require_unmodified).
#
# Usage: tools/fortran_reference/run_pes_cut_reference.sh [case] [scratch_dir]
#   case         reference case directory name (default: hardwall_small); it needs a pes_cuts.json
#   scratch_dir  where to run, a directory that does not exist yet (default: a new mktemp -d); one
#                subdirectory cut_<k> per cut
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: seconds per case. Docs: tools/fortran_reference/README.md.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

(($# <= 2)) || { echo "usage: run_pes_cut_reference.sh [case] [scratch_dir]" >&2; exit 2; }
case_dir="${FORTRAN_REPO_DIR}/tests/reference/${1:-hardwall_small}"
cuts_file="${case_dir}/pes_cuts.json"
[[ -f "${cuts_file}" ]] || refuse "${cuts_file} not found; name a case with write2DPES cuts"

# cut_lines FILE K: the DVR.inp lines of cut K (from 1) of FILE, or the number of cuts when K is 0;
# refuses (exit 1) a file that is not JSON with a list "cuts", a key that is not a write2DPES coordinate
# and a value that is not a finite number. Python's repr of a float is the shortest decimal that reads back
# to the same double, so gfortran and the C++ (which reads the JSON) get the same values.
cut_lines() {
  python3 -I -B -c '
import json, math, sys
path, k = sys.argv[1], int(sys.argv[2])
try:
    cuts = json.load(open(path, encoding="utf-8"))["cuts"]
    assert isinstance(cuts, list) and all(isinstance(cut, dict) for cut in cuts)
except (ValueError, KeyError, TypeError, AssertionError):
    sys.exit(f"error: {path} is not JSON with a list \"cuts\" of objects; fix the file")
if k == 0:
    print(len(cuts))
    sys.exit(0)
for key, value in cuts[k - 1].items():
    if key not in ("x2DPES", "y2DPES", "theta2DPES", "phi2DPES"):
        sys.exit(f"error: {path}: cut {k}: {key} is not a write2DPES coordinate; use x2DPES, y2DPES, theta2DPES"
                 " or phi2DPES")
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        sys.exit(f"error: {path}: cut {k}: {key} = {value!r} is not a finite number; give a number")
    print(f"{key} : REAL : {float(value)!r};")
' "$1" "$2"
}

n_cuts="$(cut_lines "${cuts_file}" 0)"
((n_cuts > 0)) || refuse "${cuts_file} lists no cut; add one"
stale="${case_dir}/pes_cut_$((n_cuts + 1)).dat"
[[ ! -e "${stale}" ]] ||
  refuse "${stale} is left from more cuts than ${cuts_file} lists; delete it and its MANIFEST.tsv line"
# The appended keys must be the only ones: a second line of a key in DVR.inp would be ambiguous.
if [[ -f "${case_dir}/DVR.inp" ]] &&
  grep -qiE '^[[:space:]]*(write2DPES|x2DPES|y2DPES|theta2DPES|phi2DPES)[[:space:]]*:' "${case_dir}/DVR.inp"; then
  refuse "${case_dir}/DVR.inp sets a write2DPES key itself; remove it (pes_cuts.json gives the cuts)"
fi
# Every cut is read before the first run, so that a bad one stops the script before anything is written.
cut_inputs=()
outputs=()
for ((k = 1; k <= n_cuts; ++k)); do
  lines="$(cut_lines "${cuts_file}" "${k}")"
  cut_inputs+=("${lines}")
  outputs+=("pes_cut_${k}.dat")
done
require_unmodified "${outputs[@]/#/${case_dir}/}"
scratch_dir="$(new_scratch_dir "${2:-}")"

# All runs first (a cut the Fortran rejects stops it with exit 0 and no 2dpes.dat: run_fortran refuses
# that), then the copies and their MANIFEST lines, so that a refused cut leaves the case as it was.
for ((k = 1; k <= n_cuts; ++k)); do
  run_dir="${scratch_dir}/cut_${k}"
  mkdir "${run_dir}"
  stage_inputs "${case_dir}" "${run_dir}"
  {
    echo
    echo "# Appended by run_pes_cut_reference.sh: cut ${k} of pes_cuts.json."
    echo "write2DPES : LOGICAL : TRUE;"
    [[ -z "${cut_inputs[k - 1]}" ]] || echo "${cut_inputs[k - 1]}"
  } >> "${run_dir}/DVR.inp"
  run_fortran "${run_dir}" fortran_pes_cut_stdout.txt 2dpes.dat
done
for ((k = 1; k <= n_cuts; ++k)); do
  run_dir="${scratch_dir}/cut_${k}"
  cp "${run_dir}/2dpes.dat" "${case_dir}/pes_cut_${k}.dat"
  mapfile -t inputs < <(staged_inputs "${run_dir}")
  record_manifest "${run_dir}" "${case_dir}" "${FORTRAN_REF_EXE}" "${inputs[@]}" -- "pes_cut_${k}.dat"
done
echo "* Reference files written to ${case_dir}"

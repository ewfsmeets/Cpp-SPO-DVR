#!/usr/bin/env bash
# Generates the Fortran reference data for the propagation gate (gate 3).
#
# Runs the reference executable on tests/reference/<case>/DVR.inp unchanged (a full new run:
# setup, propagation, analysis) in a scratch directory and copies into tests/reference/<case>/:
#   fortran_stdout.txt     stdout of the run (initial wave packet info, averages vs time, ...)
#   fortran_averages.out   averages.out (time, norm, <Z>, <r>, rest norm)
#   fortran_fort21.bin     fort.21, the BK coefficients of every analysis step (gate 4)
#   fortran_fort25.txt     fort.25, the state-resolved S-matrix elements and probabilities
# All other Fortran output files stay in the scratch directory.
#
# Usage: tools/fortran_reference/run_propagation_reference.sh [case] [scratch_dir]
#   case         reference case directory name (default: hardwall_short)
#   scratch_dir  where to run (default: a new mktemp -d directory); never the Fortran tree.
# OMP_NUM_THREADS is passed through (default 4).
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
case_name="${1:-hardwall_short}"
case_dir="${repo_dir}/tests/reference/${case_name}"
scratch_dir="${2:-$(mktemp -d)}"
fortran_exe="${SPODVR_FORTRAN_EXE:-${HOME}/work/SPO-DVR/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x}"

[[ -f "${case_dir}/DVR.inp" ]] || { echo "error: ${case_dir}/DVR.inp not found" >&2; exit 1; }
[[ -x "${fortran_exe}" ]] || { echo "error: Fortran executable ${fortran_exe} not found" >&2; exit 1; }

mkdir -p "${scratch_dir}"
cp "${case_dir}/DVR.inp" "${scratch_dir}/DVR.inp"

echo "* Running ${fortran_exe} in ${scratch_dir}"
start=$(date +%s)
(cd "${scratch_dir}" && OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" "${fortran_exe}" > fortran_stdout.txt 2>&1) || {
  echo "error: Fortran run failed; see ${scratch_dir}/fortran_stdout.txt" >&2; exit 1; }
echo "* Run time: $(( $(date +%s) - start )) s"

[[ -f "${scratch_dir}/averages.out" ]] || { echo "error: no averages.out written" >&2; exit 1; }
cp "${scratch_dir}/fortran_stdout.txt" "${case_dir}/fortran_stdout.txt"
cp "${scratch_dir}/averages.out" "${case_dir}/fortran_averages.out"
[[ -f "${scratch_dir}/fort.21" ]] && cp "${scratch_dir}/fort.21" "${case_dir}/fortran_fort21.bin"
[[ -f "${scratch_dir}/fort.25" ]] && cp "${scratch_dir}/fort.25" "${case_dir}/fortran_fort25.txt"
echo "* Reference files written to ${case_dir}"

#!/usr/bin/env bash
# Generates the Fortran reference for the restart gate: the hardwall_small case split into two
# runs of 3500 steps, the second a Fortran restart (NewCalc = FALSE) from the first
# (fort.9 -> fort.8, fort.21 -> fort.22, as the Fortran code expects).
# Writes into tests/reference/hardwall_small/:
#   fortran_restart_stdout.txt   stdout of the restarted run (its probabilities cover 0 .. 14000)
# The wave function files stay in the scratch directory (they are large); a C++ restart from
# them can be tried by hand with "RestartFormat": "fortran" and "RestartDir": <scratch>/part1.
#
# Usage: tools/fortran_reference/run_restart_reference.sh [scratch_dir]
# OMP_NUM_THREADS is passed through (default 4).
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
case_dir="${repo_dir}/tests/reference/hardwall_small"
scratch_dir="${1:-$(mktemp -d)}"
fortran_exe="${SPODVR_FORTRAN_EXE:-${HOME}/work/SPO-DVR/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x}"
[[ -x "${fortran_exe}" ]] || { echo "error: Fortran executable ${fortran_exe} not found" >&2; exit 1; }

# set_field <file> <name> <value>: replaces the value of "name : TYPE : value;".
set_field() {
  sed -i -E "s/^($2[[:space:]]*:[[:space:]]*[A-Z]+[[:space:]]*:).*;/\1 $3;/" "$1"
  grep -q -E "^$2[[:space:]]*:.*: $3;" "$1" || { echo "error: could not set $2 in $1" >&2; exit 1; }
}

run_part() {
  local dir="$1"
  echo "* Running ${fortran_exe} in ${dir}"
  (cd "${dir}" && OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" "${fortran_exe}" > fortran_stdout.txt 2>&1) || {
    echo "error: Fortran run failed; see ${dir}/fortran_stdout.txt" >&2; exit 1; }
}

mkdir -p "${scratch_dir}/part1" "${scratch_dir}/part2"
cp "${case_dir}/DVR.inp" "${scratch_dir}/part1/DVR.inp"
set_field "${scratch_dir}/part1/DVR.inp" tEnd 7000.0
set_field "${scratch_dir}/part1/DVR.inp" nStep 3500
run_part "${scratch_dir}/part1"

cp "${case_dir}/DVR.inp" "${scratch_dir}/part2/DVR.inp"
set_field "${scratch_dir}/part2/DVR.inp" NewCalc FALSE
set_field "${scratch_dir}/part2/DVR.inp" tStart 7000.0
set_field "${scratch_dir}/part2/DVR.inp" nStep 3500
cp "${scratch_dir}/part1/fort.9" "${scratch_dir}/part2/fort.8"
cp "${scratch_dir}/part1/fort.21" "${scratch_dir}/part2/fort.22"
run_part "${scratch_dir}/part2"

cp "${scratch_dir}/part2/fortran_stdout.txt" "${case_dir}/fortran_restart_stdout.txt"
echo "* Reference file written to ${case_dir}/fortran_restart_stdout.txt (scratch: ${scratch_dir})"

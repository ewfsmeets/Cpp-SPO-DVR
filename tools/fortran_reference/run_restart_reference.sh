#!/usr/bin/env bash
# Generates the Fortran reference for the restart gate: the hardwall_small case split into two
# runs of 3500 steps, the second a Fortran restart (NewCalc = FALSE) from the first
# (fort.9 -> fort.8, fort.21 -> fort.22, as the Fortran code expects).
# Writes into tests/reference/hardwall_small/:
#   fortran_restart_stdout.txt   stdout of the restarted run (its probabilities cover 0 .. 14000)
# and its line in the case's MANIFEST.tsv (inputs: the DVR.inp of both runs).
# The wave function files stay in the scratch directory (they are large); a C++ restart from
# them can be tried by hand with "RestartFormat": "fortran" and "RestartDir": <scratch>/part1.
#
# Usage: tools/fortran_reference/run_restart_reference.sh [scratch_dir]
#   scratch_dir  where to run, a directory that does not exist yet (default: a new mktemp -d)
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: minutes. Docs: tools/fortran_reference/README.md.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

(($# <= 1)) || { echo "usage: run_restart_reference.sh [scratch_dir]" >&2; exit 2; }
case_dir="${FORTRAN_REPO_DIR}/tests/reference/hardwall_small"
scratch_dir="$(new_scratch_dir "${1:-}")"
mkdir "${scratch_dir}/part1" "${scratch_dir}/part2"

stage_inputs "${case_dir}" "${scratch_dir}/part1"
set_field "${scratch_dir}/part1/DVR.inp" tEnd 7000.0
set_field "${scratch_dir}/part1/DVR.inp" nStep 3500
run_fortran "${scratch_dir}/part1" fortran_stdout.txt fort.9 fort.21

stage_inputs "${case_dir}" "${scratch_dir}/part2"
set_field "${scratch_dir}/part2/DVR.inp" NewCalc FALSE
set_field "${scratch_dir}/part2/DVR.inp" tStart 7000.0
set_field "${scratch_dir}/part2/DVR.inp" nStep 3500
cp "${scratch_dir}/part1/fort.9" "${scratch_dir}/part2/fort.8"
cp "${scratch_dir}/part1/fort.21" "${scratch_dir}/part2/fort.22"
run_fortran "${scratch_dir}/part2" fortran_stdout.txt averages.out

cp "${scratch_dir}/part2/fortran_stdout.txt" "${case_dir}/fortran_restart_stdout.txt"
record_manifest "${scratch_dir}" "${case_dir}" "${FORTRAN_REF_EXE}" part1/DVR.inp part2/DVR.inp -- \
  fortran_restart_stdout.txt
echo "* Reference file written to ${case_dir}/fortran_restart_stdout.txt (scratch: ${scratch_dir})"

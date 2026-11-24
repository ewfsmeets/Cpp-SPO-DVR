#!/usr/bin/env bash
# Generates the Fortran reference data for the propagation and on-the-fly analysis gates (gates 3 and 4).
#
# Runs the reference executable (reference_config.sh) on tests/reference/<case>/DVR.inp, unchanged
# except the data-file paths (reference_config.sh stage_inputs), as a full new run (setup,
# propagation, analysis) in a scratch directory and copies into tests/reference/<case>/:
#   fortran_stdout.txt     stdout of the run (initial wave packet info, averages vs time, ...)
#   fortran_averages.out   averages.out (time, norm, <Z>, <r>, rest norm)
#   fortran_fort21.bin     fort.21, the BK coefficients of every analysis step (gate 4)
#   fortran_fort25.txt     fort.25, the state-resolved S-matrix elements and probabilities
# and their lines in the case's MANIFEST.tsv (inputs: the DVR.inp and the staged data files of CRP
# and SCM cases, reference_config.sh stage_inputs). All other Fortran output files stay in the
# scratch directory. Refuses, before the run, to overwrite a file of these that differs from git
# (require_unmodified).
#
# Usage: tools/fortran_reference/run_propagation_reference.sh [case] [scratch_dir]
#   case         reference case directory name (default: hardwall_short)
#   scratch_dir  where to run, a directory that does not exist yet (default: a new mktemp -d)
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: minutes; cu111_dailight_small (7000 steps) 1570 s, measured 2026-11-11 (README.md).
# Docs: tools/fortran_reference/README.md.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

(($# <= 2)) || { echo "usage: run_propagation_reference.sh [case] [scratch_dir]" >&2; exit 2; }
case_dir="${FORTRAN_REPO_DIR}/tests/reference/${1:-hardwall_short}"
outputs=(fortran_stdout.txt fortran_averages.out fortran_fort21.bin fortran_fort25.txt)
require_unmodified "${outputs[@]/#/${case_dir}/}"
scratch_dir="$(new_scratch_dir "${2:-}")"

stage_inputs "${case_dir}" "${scratch_dir}"
start=$(date +%s)
run_fortran "${scratch_dir}" fortran_stdout.txt averages.out fort.21 fort.25
echo "* Run time: $(($(date +%s) - start)) s"

cp "${scratch_dir}/fortran_stdout.txt" "${case_dir}/fortran_stdout.txt"
cp "${scratch_dir}/averages.out" "${case_dir}/fortran_averages.out"
cp "${scratch_dir}/fort.21" "${case_dir}/fortran_fort21.bin"
cp "${scratch_dir}/fort.25" "${case_dir}/fortran_fort25.txt"
mapfile -t inputs < <(staged_inputs "${scratch_dir}")
record_manifest "${scratch_dir}" "${case_dir}" "${FORTRAN_REF_EXE}" "${inputs[@]}" -- "${outputs[@]}"
echo "* Reference files written to ${case_dir}"

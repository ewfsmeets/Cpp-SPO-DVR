#!/usr/bin/env bash
# Generates the Fortran reference data for the setup / potential / gas-phase gates.
#
# Runs the reference executable (reference_config.sh) in potential-only mode (nread = 6: stops
# after SetupAsympWF and the rot-dif basis, momentum-shift and absorber set-ups) with DoDump
# enabled, in a scratch directory, and copies the results
# into tests/reference/<case>/:
#   fortran_setup_stdout.txt  stdout of the run (vasy table, grid info, e(v,j) table)
#   debug_asyp.dat            full-precision rovibrational energies and eigenvectors
#   pes_reference.txt         extracted from the binary .pes file (see extract_pes.py)
# and their lines in the case's MANIFEST.tsv (inputs: the DVR.inp as run and the staged data files).
# Refuses, before the run, to overwrite a file of these that differs from git (require_unmodified).
# The data files of CRP and SCM cases are staged under short names (reference_config.sh stage_inputs).
#
# Usage: tools/fortran_reference/run_setup_reference.sh [case] [scratch_dir]
#   case         reference case directory name (default: hardwall_small)
#   scratch_dir  where to run, a directory that does not exist yet (default: a new mktemp -d)
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: seconds (the small cases). Docs: tools/fortran_reference/README.md.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

(($# <= 2)) || { echo "usage: run_setup_reference.sh [case] [scratch_dir]" >&2; exit 2; }
case_dir="${FORTRAN_REPO_DIR}/tests/reference/${1:-hardwall_small}"
outputs=(fortran_setup_stdout.txt debug_asyp.dat pes_reference.txt)
require_unmodified "${outputs[@]/#/${case_dir}/}"
scratch_dir="$(new_scratch_dir "${2:-}")"

stage_inputs "${case_dir}" "${scratch_dir}"
cat >> "${scratch_dir}/DVR.inp" <<'INP'

# Appended by run_setup_reference.sh: potential-only run with full-precision gas-phase dump.
Potential   : LOGICAL : TRUE;
DoDump      : LOGICAL : TRUE;
INP
run_fortran "${scratch_dir}" fortran_setup_stdout.txt '*.pes' debug_asyp.dat

pes_files=("${scratch_dir}"/*.pes)
python3 -I -B "${FORTRAN_REPO_DIR}/tools/fortran_reference/extract_pes.py" "${pes_files[0]}" "${case_dir}/pes_reference.txt"
cp "${scratch_dir}/fortran_setup_stdout.txt" "${scratch_dir}/debug_asyp.dat" "${case_dir}/"
mapfile -t inputs < <(staged_inputs "${scratch_dir}")
record_manifest "${scratch_dir}" "${case_dir}" "${FORTRAN_REF_EXE}" "${inputs[@]}" -- "${outputs[@]}"
echo "* Reference files written to ${case_dir}"

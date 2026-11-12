#!/usr/bin/env bash
# Generates the Fortran reference data for the setup / potential / gas-phase gates.
#
# Runs the reference executable (reference_config.sh) in potential-only mode (nread = 6: stops
# right after SetupAsympWF) with DoDump enabled, in a scratch directory, and copies the results
# into tests/reference/<case>/:
#   fortran_setup_stdout.txt  stdout of the run (vasy table, grid info, e(v,j) table)
#   debug_asyp.dat            full-precision rovibrational energies and eigenvectors
#   pes_reference.txt         extracted from the binary .pes file (see extract_pes.py)
# and their lines in the case's MANIFEST.tsv (inputs: the DVR.inp as run and a staged data set).
# CRP library cases (SystemID CRPlibrary) get their data set copied into the scratch directory.
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
inputs=(DVR.inp)
[[ ! -d "${scratch_dir}/dataset" ]] || inputs+=(dataset)
record_manifest "${scratch_dir}" "${case_dir}" "${FORTRAN_REF_EXE}" "${inputs[@]}" -- \
  fortran_setup_stdout.txt debug_asyp.dat pes_reference.txt
echo "* Reference files written to ${case_dir}"

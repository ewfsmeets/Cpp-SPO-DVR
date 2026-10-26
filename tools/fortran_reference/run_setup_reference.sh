#!/usr/bin/env bash
# Generates the Fortran reference data for the setup / potential / gas-phase gates.
#
# Runs the reference executable in potential-only mode (nread = 6: stops right after
# SetupAsympWF) with DoDump enabled, in a scratch directory, and copies the results into
# tests/reference/<case>/:
#   fortran_setup_stdout.txt  stdout of the run (vasy table, grid info, e(v,j) table)
#   debug_asyp.dat            full-precision rovibrational energies and eigenvectors
#   pes_reference.txt         extracted from the binary .pes file (see extract_pes.py)
#
# Usage: tools/fortran_reference/run_setup_reference.sh [case] [scratch_dir]
#   case         reference case directory name (default: hardwall_small)
#   scratch_dir  where to run (default: a new mktemp -d directory); never the Fortran tree.
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
case_name="${1:-hardwall_small}"
case_dir="${repo_dir}/tests/reference/${case_name}"
scratch_dir="${2:-$(mktemp -d)}"
fortran_exe="${SPODVR_FORTRAN_EXE:-${HOME}/work/SPO-DVR/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x}"

[[ -f "${case_dir}/DVR.inp" ]] || { echo "error: ${case_dir}/DVR.inp not found" >&2; exit 1; }
[[ -x "${fortran_exe}" ]] || { echo "error: Fortran executable ${fortran_exe} not found" >&2; exit 1; }

mkdir -p "${scratch_dir}"
cp "${case_dir}/DVR.inp" "${scratch_dir}/DVR.inp"
cat >> "${scratch_dir}/DVR.inp" <<'INP'

# Appended by run_setup_reference.sh: potential-only run with full-precision gas-phase dump.
Potential   : LOGICAL : TRUE;
DoDump      : LOGICAL : TRUE;
INP

echo "* Running ${fortran_exe} in ${scratch_dir}"
(cd "${scratch_dir}" && OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" "${fortran_exe}" > fortran_setup_stdout.txt 2>&1) || {
  echo "error: Fortran run failed; see ${scratch_dir}/fortran_setup_stdout.txt" >&2; exit 1; }

pes_file="$(ls "${scratch_dir}"/*.pes 2>/dev/null | head -n 1 || true)"
[[ -n "${pes_file}" ]] || { echo "error: no .pes file written in ${scratch_dir}" >&2; exit 1; }
[[ -f "${scratch_dir}/debug_asyp.dat" ]] || { echo "error: no debug_asyp.dat written" >&2; exit 1; }

python3 -I "${repo_dir}/tools/fortran_reference/extract_pes.py" "${pes_file}" "${case_dir}/pes_reference.txt"
cp "${scratch_dir}/fortran_setup_stdout.txt" "${scratch_dir}/debug_asyp.dat" "${case_dir}/"
echo "* Reference files written to ${case_dir}"

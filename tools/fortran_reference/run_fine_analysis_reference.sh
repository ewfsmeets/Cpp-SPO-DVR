#!/usr/bin/env bash
# Generates the Fortran reference data for the fine-analysis gate (gate 5, run mode 2).
#
# Runs the unmodified reference executable as an analysis-only run (onlyAnalysis, nread = 2)
# on the BK coefficients of the full Fortran run (tests/reference/<case>/fortran_fort21.bin,
# copied to fort.22). Only DVR.inp is changed: onlyAnalysis and writeSMatrix are switched on,
# and some variants change tEnd / nStep / NrOfEs:
#   full       tEnd unchanged (1400 records; ndused = 1386 < ndatol)
#   truncated  tEnd 13860, nStep 6930 (only the first 1386 records are used; ndused = ndatol)
#   nen95      as truncated, with NrOfEs = -95: the analysis energies fall on the fine grid
#   small      tEnd 800, nStep 400 (80 records, a few energies; coefficients are still noise)
#   empty      tEnd 100, nStep 50 (10 records: no energy in the confidence interval, nen95 = 0)
# For each variant the files go to tests/reference/<case>/fine_<variant>/:
#   stdout.txt  stdout (asymptotic info, probaf table, interpolated per-energy blocks)
#   fort.10     fil10 (spectrum of the initial channel, reaction probability, averages, sums)
#   fort.3      unit 3 (probability of every state on the fine grid; writeSMatrix)
#   fort.25     the interpolated S-matrix elements, probabilities and phases (writeProbs)
#
# Usage: tools/fortran_reference/run_fine_analysis_reference.sh [case] [variant ...]
#   case     reference case directory name (default: hardwall_small)
#   variant  one or more of the variants above (default: all)
# Every variant runs in its own new mktemp -d directory (never the Fortran tree), in seconds:
# an analysis-only run computes no 6D potential and does not propagate.
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
case_name="${1:-hardwall_small}"
shift || true
variants=("$@")
[[ ${#variants[@]} -gt 0 ]] || variants=(full truncated nen95 small empty)
case_dir="${repo_dir}/tests/reference/${case_name}"
fortran_exe="${SPODVR_FORTRAN_EXE:-${HOME}/work/SPO-DVR/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x}"

[[ -f "${case_dir}/DVR.inp" ]] || { echo "error: ${case_dir}/DVR.inp not found" >&2; exit 1; }
[[ -f "${case_dir}/fortran_fort21.bin" ]] || {
  echo "error: ${case_dir}/fortran_fort21.bin not found (run run_propagation_reference.sh ${case_name})" >&2; exit 1; }
[[ -x "${fortran_exe}" ]] || { echo "error: Fortran executable ${fortran_exe} not found" >&2; exit 1; }

# set_field FILE NAME VALUE: replaces the value of an existing "Name : TYPE : value;" line.
set_field() {
  grep -q "^$2 *:" "$1" || { echo "error: field $2 not found in $1" >&2; exit 1; }
  sed -i -E "s/^($2 *:[^:]*: *)[^;]*;/\1$3;/" "$1"
}

for variant in "${variants[@]}"; do
  scratch_dir="$(mktemp -d)"
  inp="${scratch_dir}/DVR.inp"
  cp "${case_dir}/DVR.inp" "${inp}"
  cp "${case_dir}/fortran_fort21.bin" "${scratch_dir}/fort.22"
  set_field "${inp}" writeSMatrix TRUE
  echo 'onlyAnalysis: LOGICAL : TRUE;' >> "${inp}"
  case "${variant}" in
    full) ;;
    truncated) set_field "${inp}" tEnd 13860.0; set_field "${inp}" nStep 6930 ;;
    nen95) set_field "${inp}" tEnd 13860.0; set_field "${inp}" nStep 6930; set_field "${inp}" NrOfEs -95 ;;
    small) set_field "${inp}" tEnd 800.0; set_field "${inp}" nStep 400 ;;
    empty) set_field "${inp}" tEnd 100.0; set_field "${inp}" nStep 50 ;;
    *) echo "error: unknown variant ${variant}" >&2; exit 1 ;;
  esac

  echo "* ${variant}: running in ${scratch_dir}"
  (cd "${scratch_dir}" && OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" "${fortran_exe}" > stdout.txt 2>&1) || {
    echo "error: Fortran run failed; see ${scratch_dir}/stdout.txt" >&2; exit 1; }

  out_dir="${case_dir}/fine_${variant}"
  mkdir -p "${out_dir}"
  cp "${scratch_dir}/DVR.inp" "${scratch_dir}/stdout.txt" "${out_dir}/"
  for f in fort.10 fort.3 fort.25; do
    [[ -f "${scratch_dir}/${f}" ]] && cp "${scratch_dir}/${f}" "${out_dir}/"
  done
  echo "  written: $(cd "${out_dir}" && ls | tr '\n' ' ')"
done

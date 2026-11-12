#!/usr/bin/env bash
# Generates the Fortran reference data for the fine-analysis gate (gate 5, run mode 2).
#
# Runs the unmodified reference executable (reference_config.sh) as an analysis-only run
# (onlyAnalysis, nread = 2) on the BK coefficients of the full Fortran run
# (tests/reference/<case>/fortran_fort21.bin, copied to fort.22). Only DVR.inp is changed:
# onlyAnalysis and writeSMatrix are switched on, and some variants change tEnd / nStep / NrOfEs:
#   full       tEnd unchanged (1400 records; ndused = 1386 < ndatol)
#   truncated  tEnd 13860, nStep 6930 (only the first 1386 records are used; ndused = ndatol)
#   nen95      as truncated, with NrOfEs = -95: the analysis energies fall on the fine grid
#   small      tEnd 800, nStep 400 (80 records, a few energies; coefficients are still noise)
#   empty      tEnd 100, nStep 50 (10 records: no energy in the confidence interval, nen95 = 0;
#              so no fort.3 and no fort.25)
# For each variant the files go to tests/reference/<case>/fine_<variant>/:
#   DVR.inp     the input as run
#   stdout.txt  stdout (asymptotic info, probaf table, interpolated per-energy blocks)
#   fort.10     fil10 (spectrum of the initial channel, reaction probability, averages, sums)
#   fort.3      unit 3 (probability of every state on the fine grid; writeSMatrix)
#   fort.25     the interpolated S-matrix elements, probabilities and phases (writeProbs)
# and their lines in the case's MANIFEST.tsv (inputs: the DVR.inp as run and fort.22, whose hash
# is that of fortran_fort21.bin).
#
# Usage: tools/fortran_reference/run_fine_analysis_reference.sh [case] [variant ...]
#   case     reference case directory name (default: hardwall_small)
#   variant  one or more of the variants above (default: all)
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: seconds per variant (an analysis-only run computes no 6D potential and does not
# propagate); each runs in its own new mktemp -d directory. Docs: tools/fortran_reference/README.md.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

case_name="${1:-hardwall_small}"
shift || true
variants=("$@")
[[ ${#variants[@]} -gt 0 ]] || variants=(full truncated nen95 small empty)
case_dir="${FORTRAN_REPO_DIR}/tests/reference/${case_name}"
for variant in "${variants[@]}"; do
  [[ "${variant}" =~ ^(full|truncated|nen95|small|empty)$ ]] || {
    echo "usage: run_fine_analysis_reference.sh [case] [full|truncated|nen95|small|empty ...]" >&2
    exit 2
  }
done
[[ -f "${case_dir}/fortran_fort21.bin" ]] ||
  refuse "${case_dir}/fortran_fort21.bin not found; run run_propagation_reference.sh ${case_name} first"

for variant in "${variants[@]}"; do
  scratch_dir="$(new_scratch_dir)"
  inp="${scratch_dir}/DVR.inp"
  stage_inputs "${case_dir}" "${scratch_dir}"
  cp "${case_dir}/fortran_fort21.bin" "${scratch_dir}/fort.22"
  set_field "${inp}" writeSMatrix TRUE
  echo 'onlyAnalysis: LOGICAL : TRUE;' >> "${inp}"
  outputs=(fort.10 fort.3 fort.25)
  case "${variant}" in
    full) ;;
    truncated) set_field "${inp}" tEnd 13860.0; set_field "${inp}" nStep 6930 ;;
    nen95) set_field "${inp}" tEnd 13860.0; set_field "${inp}" nStep 6930; set_field "${inp}" NrOfEs -95 ;;
    small) set_field "${inp}" tEnd 800.0; set_field "${inp}" nStep 400 ;;
    empty) set_field "${inp}" tEnd 100.0; set_field "${inp}" nStep 50; outputs=(fort.10) ;;
  esac

  echo "* ${variant}"
  run_fortran "${scratch_dir}" stdout.txt "${outputs[@]}"
  out_dir="${case_dir}/fine_${variant}"
  mkdir -p "${out_dir}"
  files=(DVR.inp stdout.txt "${outputs[@]}")
  for f in "${files[@]}"; do cp "${scratch_dir}/${f}" "${out_dir}/"; done
  record_manifest "${scratch_dir}" "${case_dir}" "${FORTRAN_REF_EXE}" DVR.inp fort.22 -- \
    "${files[@]/#/fine_${variant}/}"
  echo "  written: DVR.inp stdout.txt ${outputs[*]}"
done

#!/usr/bin/env bash
# The Fortran part of the flux replay (docs/validation_flux_replay/README.md): the reference
# executable runs hardwall_small with the flux keys in two parts, split by a restart, and a control.
#   part1         t = 0 .. 7000 (3500 steps), writes fort.9, fort.21 and the flux coefficient files;
#   part2         t = 7000 .. 7010 (one analysis step), restarted from part1, with part1's flux
#                 coefficient files zero-filled (zero_payload.py): its flux_coeffs.dat and
#                 flux_coeffs_deriv.dat then hold c_n psi_c and c_n d psi_c / dr of that one step;
#   part2_noflux  part2 without the flux keys: its fort.21 and fort.9 must equal part2's (checked).
# Output: the three run directories in SCRATCH_DIR (a new directory; never inside the repository).
# Usage: run_fortran.sh SCRATCH_DIR
# Environment (read by tools/fortran_reference/reference_config.sh): SPODVR_FORTRAN_REF_DIR (the
# reference install, default the pinned one), SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT.
# Run time: about 2 minutes (4 threads, the thread count of every reference run).
# Docs: docs/validation_flux_replay/README.md, tools/fortran_reference/README.md.
set -euo pipefail

if (($# != 1)); then
  echo "usage: run_fortran.sh SCRATCH_DIR" >&2
  exit 2
fi
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/../../tools/fortran_reference/reference_config.sh"
case_dir="${FORTRAN_REPO_DIR}/tests/reference/hardwall_small"
scratch="$(new_scratch_dir "$1")"
mkdir "${scratch}/part1" "${scratch}/part2" "${scratch}/part2_noflux"

# add_flux_keys DVR.inp: the five flux keys, typed (InputNew.f90 ReadInput pre-registers none).
# replay_flux.cpp (part2_config) uses the same values.
add_flux_keys() {
  printf '%s\n' 'doFluxAnalysis : LOGICAL : TRUE;' 'rFluxAnalysis : REAL : 2.0;' 'EMinFluxAna : REAL : 0.2;' \
    'EMaxFluxAna : REAL : 0.6;' 'DeltaEFluxAna : REAL : 0.15;' >> "$1"
}

stage_inputs "${case_dir}" "${scratch}/part1"
set_field "${scratch}/part1/DVR.inp" tEnd 7000.0
set_field "${scratch}/part1/DVR.inp" nStep 3500
add_flux_keys "${scratch}/part1/DVR.inp"
run_fortran "${scratch}/part1" fortran_stdout.txt fort.9 fort.21 flux_coeffs.dat flux_coeffs_deriv.dat

for part in part2 part2_noflux; do
  stage_inputs "${case_dir}" "${scratch}/${part}"
  set_field "${scratch}/${part}/DVR.inp" NewCalc FALSE
  set_field "${scratch}/${part}/DVR.inp" tStart 7000.0
  set_field "${scratch}/${part}/DVR.inp" tEnd 7010.0
  set_field "${scratch}/${part}/DVR.inp" nStep 5
  cp "${scratch}/part1/fort.9" "${scratch}/${part}/fort.8"
  cp "${scratch}/part1/fort.21" "${scratch}/${part}/fort.22"
done
add_flux_keys "${scratch}/part2/DVR.inp"
for name in flux_coeffs.dat flux_coeffs_deriv.dat; do
  python3 -I -B "${here}/zero_payload.py" "${scratch}/part1/${name}" "${scratch}/part2/${name}"
done
run_fortran "${scratch}/part2" fortran_stdout.txt fort.9 fort.21 flux_coeffs.dat flux_coeffs_deriv.dat
run_fortran "${scratch}/part2_noflux" fortran_stdout.txt fort.9 fort.21

# The flux analysis is passive: the step with zero-filled flux files equals the step without flux.
for name in fort.21 fort.9; do
  cmp -s "${scratch}/part2/${name}" "${scratch}/part2_noflux/${name}" ||
    refuse "part2 ${name} differs from part2_noflux ${name}: the flux run changed the propagation; investigate"
  echo "check: part2 ${name} with zero-filled flux files == part2_noflux ${name}: identical"
done
echo "runs: ${scratch}"

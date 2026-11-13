#!/usr/bin/env bash
# Generates the Fortran reference data for the propagation gates (3, 4, and 6 for a flux case).
#
# Runs the reference executable (reference_config.sh) on tests/reference/<case>/DVR.inp unchanged
# (a full new run: setup, propagation, analysis) in a scratch directory and copies into
# tests/reference/<case>/:
#   fortran_stdout.txt     stdout of the run (initial wave packet info, averages vs time, ...)
#   fortran_averages.out   averages.out (time, norm, <Z>, <r>, rest norm)
#   fortran_fort21.bin     fort.21, the BK coefficients of every analysis step (gate 4)
#   fortran_fort25.txt     fort.25, the state-resolved S-matrix elements and probabilities
# and, for a case with doFluxAnalysis (gate 6; the Fortran writes these files only then):
#   fortran_flux_coeffs.bin        flux_coeffs.dat, phi(E) at the flux surface (gfortran unformatted)
#   fortran_flux_coeffs_deriv.bin  flux_coeffs_deriv.dat, phi'(E), its r-derivative
#   fortran_flux_reaction.dat      flux_reaction.dat, the reaction probability P(E)
#   fortran_flux_density.dat       flux_density.dat, the probability density P(x, y; E)
# and their lines in the case's MANIFEST.tsv. All other Fortran output files stay in the scratch
# directory.
#
# Usage: tools/fortran_reference/run_propagation_reference.sh [case] [scratch_dir]
#   case         reference case directory name (default: hardwall_short)
#   scratch_dir  where to run, a directory that does not exist yet (default: a new mktemp -d)
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: minutes; the 7000-step Dai-Light case about half an hour.
# Docs: tools/fortran_reference/README.md.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

(($# <= 2)) || { echo "usage: run_propagation_reference.sh [case] [scratch_dir]" >&2; exit 2; }
case_dir="${FORTRAN_REPO_DIR}/tests/reference/${1:-hardwall_short}"
scratch_dir="$(new_scratch_dir "${2:-}")"

stage_inputs "${case_dir}" "${scratch_dir}"
start=$(date +%s)
run_fortran "${scratch_dir}" fortran_stdout.txt averages.out fort.21 fort.25
echo "* Run time: $(($(date +%s) - start)) s"

cp "${scratch_dir}/fortran_stdout.txt" "${case_dir}/fortran_stdout.txt"
cp "${scratch_dir}/averages.out" "${case_dir}/fortran_averages.out"
cp "${scratch_dir}/fort.21" "${case_dir}/fortran_fort21.bin"
cp "${scratch_dir}/fort.25" "${case_dir}/fortran_fort25.txt"
# The flux files, if the run wrote them; regenerate_all.sh requires them of each flux case.
flux_files=()
if [[ -e "${scratch_dir}/flux_reaction.dat" ]]; then
  cp "${scratch_dir}/flux_coeffs.dat" "${case_dir}/fortran_flux_coeffs.bin"
  cp "${scratch_dir}/flux_coeffs_deriv.dat" "${case_dir}/fortran_flux_coeffs_deriv.bin"
  cp "${scratch_dir}/flux_reaction.dat" "${case_dir}/fortran_flux_reaction.dat"
  cp "${scratch_dir}/flux_density.dat" "${case_dir}/fortran_flux_density.dat"
  flux_files=(fortran_flux_coeffs.bin fortran_flux_coeffs_deriv.bin fortran_flux_reaction.dat fortran_flux_density.dat)
fi
record_manifest "${scratch_dir}" "${case_dir}" "${FORTRAN_REF_EXE}" DVR.inp -- \
  fortran_stdout.txt fortran_averages.out fortran_fort21.bin fortran_fort25.txt "${flux_files[@]}"
echo "* Reference files written to ${case_dir}"

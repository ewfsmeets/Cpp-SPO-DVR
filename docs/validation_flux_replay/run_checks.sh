#!/usr/bin/env bash
# Reruns the flux analysis checks against the reference Fortran (docs/validation_flux_replay/README.md):
#   check 1  b(k), the r-derivative weights D and the energy count (flux_f1_check.{f90,cpp});
#   check 2  the accumulation of phi, phi' and the probability on a synthetic wave function
#            (flux_f2_check.{f90,cpp});
#   replay   one analysis step of hardwall_small (run_fortran.sh, replay_flux.cpp).
# The Fortran drivers are linked against the reference install's objects (all but Main.o); the
# C++ helpers against BUILD_DIR's libraries (build_helper.py). Every output goes to SCRATCH_DIR.
# Usage: run_checks.sh BUILD_DIR SCRATCH_DIR   (BUILD_DIR: a configured and built build, e.g. build;
#        SCRATCH_DIR: a new directory outside the repository)
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Exit status: 0 every check passes, 1 a check fails or a refusal, 2 bad arguments.
# Run time: about 3 minutes (the Fortran run of 3500 steps).
set -euo pipefail

if (($# != 2)); then
  echo "usage: run_checks.sh BUILD_DIR SCRATCH_DIR" >&2
  exit 2
fi
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${here}/../../tools/fortran_reference/reference_config.sh"
build="$(cd "$1" && pwd)" || refuse "build directory $1 not found; name a configured and built one"
scratch="$(new_scratch_dir "$2")"
check_install

# The reference objects a driver links: every object of the install but the main program.
fortran_objects=()
for object in "${FORTRAN_REF_OBJS}"/*.o; do
  [[ "$(basename "${object}")" == Main.o ]] || fortran_objects+=("${object}")
done

# build_fortran_check NAME: compiles NAME.f90 with the reference build's flags and links it.
build_fortran_check() {
  # shellcheck disable=SC2086  # FORTRAN_BUILD_FLAGS is a list of flags
  run_clean gfortran ${FORTRAN_BUILD_FLAGS} -I"${FORTRAN_REF_OBJS}" -J"${scratch}/fortran" \
    -c "${here}/$1.f90" -o "${scratch}/fortran/$1.o"
  run_clean gfortran -fopenmp -o "${scratch}/fortran/$1" "${scratch}/fortran/$1.o" "${fortran_objects[@]}" \
    -lfftw3 -lfftw3_omp
}

mkdir "${scratch}/fortran" "${scratch}/f2_run"
for check in flux_f1_check flux_f2_check; do
  build_fortran_check "${check}"
  python3 -I -B "${here}/build_helper.py" "${build}" "${here}/${check}.cpp" "${scratch}/helpers"
done
python3 -I -B "${here}/build_helper.py" "${build}" "${here}/replay_flux.cpp" "${scratch}/helpers"

status=0
echo "--- check 1: b(k), D, energy count"
run_clean "${scratch}/fortran/flux_f1_check" > "${scratch}/fortran_f1.txt"
grep -i "energies" "${scratch}/fortran_f1.txt" || true
"${scratch}/helpers/flux_f1_check" < "${scratch}/fortran_f1.txt" || status=1
echo "--- check 2: phi, phi', P on a synthetic wave function"
(cd "${scratch}/f2_run" && run_clean "${scratch}/fortran/flux_f2_check" > "${scratch}/fortran_f2.txt")
"${scratch}/helpers/flux_f2_check" "${scratch}/f2_run" "${FORTRAN_REPO_DIR}/tests/reference" \
  < "${scratch}/fortran_f2.txt" || status=1
echo "--- replay: one analysis step of hardwall_small"
"${here}/run_fortran.sh" "${scratch}/runs"
"${scratch}/helpers/replay_flux" "${FORTRAN_REPO_DIR}/tests/reference/hardwall_small/input.json" \
  "${scratch}/runs/part2" "${scratch}/replay_out" || status=1
((status == 0)) && echo "all checks pass" || echo "a check FAILED"
exit "${status}"

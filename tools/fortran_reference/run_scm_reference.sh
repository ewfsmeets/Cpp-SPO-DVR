#!/usr/bin/env bash
# Generates the Fortran reference data of the SCM distorted-surface correction (tests/scm/reference/).
#
# Builds the SCM reference driver (scm_driver/build.sh, linked against the objects of the reference
# install, reference_config.sh) in a scratch directory and runs it, one process per run, on the
# committed point list points.in of a case, with the CRP data set
# potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v (SRP48):
#   example    the example files potentials/H2Cu111SRP_C3v.dir/{params.dat,examplePosfile.dat}
#              (used where they are), configuration 1, cut-off 16 bohr:
#                rescaling_on.out (EnableSurfExpan T), rescaling_off.out (F)
#   synthetic  the case's params.in (layer 2 differs from layer 1) and synthetic slab files:
#                v3_slab3_rescaling_on.out   slab_v3.in, configuration 3 of 3, cut-off 9.5 bohr, T
#                v2_slab2_rescaling_off.out  slab_v2.in, configuration 2 of 3, cut-off 4 bohr, F
# Each <run>.out has its line in the case's MANIFEST.tsv (inputs: points.in, params.dat, slab.dat
# and the data set, as staged for the run). The output format is in scm_driver/scm_reference_driver.f90.
# --generate first rewrites the case's inputs (scm_driver/make_scm_inputs.py).
#
# Usage: tools/fortran_reference/run_scm_reference.sh [--generate] example|synthetic [scratch_dir]
#   scratch_dir  where to build and run, a directory that does not exist yet (default: mktemp -d)
# Every driver run must exit with 0 and end its output with SCM_DRIVER_DONE (a Fortran STOP inside
# the library exits with status 0); the script fails otherwise and writes no output of that run.
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: seconds per case. Docs: tools/fortran_reference/README.md ("The SCM driver").
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

usage() {
  echo "usage: run_scm_reference.sh [--generate] example|synthetic [scratch_dir]" >&2
  exit 2
}
generate=false
if [[ "${1:-}" == --generate ]]; then
  generate=true
  shift
fi
(($# >= 1 && $# <= 2)) && [[ "$1" == example || "$1" == synthetic ]] || usage
target="$1"
case_dir="${FORTRAN_REPO_DIR}/tests/scm/reference/${target}"
example_dir="${FORTRAN_REPO_DIR}/potentials/H2Cu111SRP_C3v.dir"
data_set="${example_dir}/H2Cu111SRP_C3v"
scratch_dir="$(new_scratch_dir "${2:-}")"

if [[ "${generate}" == true ]]; then
  python3 -I -B "${FORTRAN_REPO_DIR}/tools/fortran_reference/scm_driver/make_scm_inputs.py" \
    "${target}" "${case_dir}" "${example_dir}"
fi
[[ -f "${case_dir}/points.in" ]] ||
  refuse "no point list ${case_dir}/points.in; restore it, or create it with --generate"
[[ -f "${data_set}/CRPData" ]] || refuse "data set ${data_set} not found; restore it from git"

driver="$("${FORTRAN_REPO_DIR}/tools/fortran_reference/scm_driver/build.sh" "${scratch_dir}/driver" | tail -n 1)"
# The CRP library opens the data set by a short path: one copy, ../dataset from every run directory.
cp -r "${data_set}" "${scratch_dir}/dataset"

# run_driver NAME PARAMS SLAB SLAB_NUMBER CUTOFF T|F: the driver on points.in with PARAMS and SLAB
# staged as params.dat and slab.dat (InitialiseDistPESModule takes paths of at most 100 characters),
# in its own directory; then <NAME>.out goes into the case directory with its MANIFEST line.
run_driver() {
  local name="$1" params="$2" slab="$3" run="${scratch_dir}/run_$1" status=0
  shift 3
  [[ -f "${params}" && -f "${slab}" ]] || refuse "input ${params} or ${slab} not found; restore it, or use --generate"
  mkdir "${run}"
  cp "${params}" "${run}/params.dat"
  cp "${slab}" "${run}/slab.dat"
  (cd "${run}" && run_clean "${driver}" ../dataset ./params.dat ./slab.dat "$@" "${name}.out" \
    < "${case_dir}/points.in" > "${run}/stdout.txt" 2>&1) || status=$?
  ((status == 0)) || refuse "the driver failed (exit ${status}) in ${run}; see ${run}/stdout.txt"
  [[ -f "${run}/${name}.out" ]] && tail -n 1 "${run}/${name}.out" | grep -q '^SCM_DRIVER_DONE$' ||
    refuse "the driver did not finish in ${run}; see ${run}/stdout.txt and ${run}/${name}.out"
  cp "${run}/${name}.out" "${case_dir}/${name}.out"
  record_manifest "${case_dir}" "${case_dir}" "${driver}" points.in "${run}/params.dat" "${run}/slab.dat" \
    "${scratch_dir}/dataset" -- "${name}.out"
}

if [[ "${target}" == example ]]; then
  run_driver rescaling_on "${example_dir}/params.dat" "${example_dir}/examplePosfile.dat" 1 16.0 T
  run_driver rescaling_off "${example_dir}/params.dat" "${example_dir}/examplePosfile.dat" 1 16.0 F
else
  # The cut-off 9.5 bohr: scm_driver/make_scm_inputs.py places three atoms of slab_v3.in at it.
  run_driver v3_slab3_rescaling_on "${case_dir}/params.in" "${case_dir}/slab_v3.in" 3 9.5 T
  # Version 2 has no cut-off in the Fortran: at 4 bohr only the atom at the origin would be kept, so
  # a reader that applies it in version 2 fails against this run.
  run_driver v2_slab2_rescaling_off "${case_dir}/params.in" "${case_dir}/slab_v2.in" 2 4.0 F
fi
echo "* References of ${target} written to ${case_dir}"

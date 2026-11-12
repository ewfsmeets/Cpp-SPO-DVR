#!/usr/bin/env bash
# Generates the Fortran reference data of the CRP library tests (tests/crplib/reference/).
#
# Builds the CRP reference driver (crp_driver/build.sh, linked against the objects of the
# reference install, reference_config.sh) in a scratch directory and runs it on the committed
# inputs of one case directory:
#   kernels   Fourier tables (tables.out) and spline kernels on the point lists kernels/*.in
#   <case>    the point lists <case>/*.in on the data set <case>/dataset (see tests/crplib/README.md);
#             cu111_srp48: on the real data set potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v
#             (not copied into the case)
# Each <name>.in gives <name>.out, with its line in the case's MANIFEST.tsv (inputs: <name>.in and
# the data set).
# --generate first rewrites the inputs: the kernel point lists (crp_driver/make_kernel_inputs.py),
# or a case's data set and point lists (build/tests/crplib/crplib_make_dataset: needs the C++ build).
#
# Usage: tools/fortran_reference/run_crp_reference.sh [--generate] kernels|<case> [scratch_dir]
#   scratch_dir  where to build and run, a directory that does not exist yet (default: mktemp -d)
# Every driver run must end with CRP_DRIVER_DONE (a Fortran STOP inside the library exits with
# status 0); the script fails otherwise.
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: seconds per case. Docs: tests/crplib/README.md, tools/fortran_reference/README.md.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

usage() {
  echo "usage: run_crp_reference.sh [--generate] kernels|<case> [scratch_dir]" >&2
  exit 2
}
generate=false
if [[ "${1:-}" == --generate ]]; then
  generate=true
  shift
fi
(($# >= 1 && $# <= 2)) && [[ "$1" != -* ]] || usage
target="$1"
case_dir="${FORTRAN_REPO_DIR}/tests/crplib/reference/${target}"
scratch_dir="$(new_scratch_dir "${2:-}")"

driver="$("${FORTRAN_REPO_DIR}/tools/fortran_reference/crp_driver/build.sh" "${scratch_dir}/driver" | tail -n 1)"

# run_driver INPUT OUTPUT MODE [ARGUMENT...]: the driver in MODE, reading INPUT, into OUTPUT.
run_driver() {
  local input="$1" output="$2"
  shift 2
  run_clean "${driver}" "$@" < "${input}" > "${output}"
  tail -n 1 "${output}" | grep -q '^CRP_DRIVER_DONE$' ||
    refuse "the driver did not finish ($*); see ${output}"
}

# require_point_lists: the case has point lists (*.in); without them no output would be written.
require_point_lists() {
  compgen -G "${case_dir}/*.in" > /dev/null ||
    refuse "no point lists (*.in) in ${case_dir}; name a case of tests/crplib/reference, or create it with --generate"
}

if [[ "${target}" == "kernels" ]]; then
  if [[ "${generate}" == true ]]; then
    mkdir -p "${case_dir}"
    python3 -I -B "${FORTRAN_REPO_DIR}/tools/fortran_reference/crp_driver/make_kernel_inputs.py" "${case_dir}"
  fi
  require_point_lists
  run_driver /dev/null "${case_dir}/tables.out" tables
  record_manifest "${case_dir}" "${case_dir}" "${driver}" -- tables.out
  for input in "${case_dir}"/*.in; do
    name="$(basename "${input}" .in)"
    run_driver "${input}" "${case_dir}/${name}.out" "${name%%_*}"
    record_manifest "${case_dir}" "${case_dir}" "${driver}" "${name}.in" -- "${name}.out"
  done
  echo "* Kernel references written to ${case_dir}"
  exit 0
fi

# Cases on a real data set of the repository: the generator writes only the point lists; the data
# set stays where it is (not copied into the case).
declare -A real_data_sets=([cu111_srp48]="potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v")
data_set="${case_dir}/dataset"
[[ -z "${real_data_sets[${target}]:-}" ]] || data_set="${FORTRAN_REPO_DIR}/${real_data_sets[${target}]}"
if [[ "${generate}" == true ]]; then
  generator="${FORTRAN_REPO_DIR}/build/tests/crplib/crplib_make_dataset"
  [[ -x "${generator}" ]] || refuse "${generator} not built; cmake --build build"
  if [[ -n "${real_data_sets[${target}]:-}" ]]; then
    "${generator}" "${target}" "${case_dir}" "${data_set}"
  else
    "${generator}" "${target}" "${case_dir}"
  fi
fi
require_point_lists
[[ -f "${data_set}/CRPData" ]] || refuse "data set ${data_set} not found; restore it, or create it with --generate"
# The Fortran opens the data set with a path of at most 50 characters: run with a short relative one.
work="${scratch_dir}/run"
mkdir "${work}" && cp -r "${data_set}" "${work}/dataset"
for input in "${case_dir}"/*.in; do
  name="$(basename "${input}" .in)"   # <mode>[_<length code>_<energy code>][_<tag>]
  IFS=_ read -r mode length energy _ <<< "${name}"
  case "${mode}" in
    molecule|atom) (cd "${work}" && run_driver "${input}" "${case_dir}/${name}.out" "${mode}" dataset "${length}" "${energy}") ;;
    *)             (cd "${work}" && run_driver "${input}" "${case_dir}/${name}.out" "${mode}" dataset) ;;
  esac
  record_manifest "${case_dir}" "${case_dir}" "${driver}" "${name}.in" "${work}/dataset" -- "${name}.out"
done
echo "* References of ${target} written to ${case_dir}"

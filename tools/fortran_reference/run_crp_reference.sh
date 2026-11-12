#!/usr/bin/env bash
# Generates the Fortran reference data of the CRP library tests (tests/crplib/reference/).
#
# Builds the CRP reference driver (crp_driver/build.sh, linked against the objects of the
# reference install, reference_config.sh) in a scratch directory and runs it:
#   kernels   Fourier tables and spline kernels: inputs from crp_driver/make_kernel_inputs.py,
#             written with the driver output to tests/crplib/reference/kernels/
#   <case>    a synthetic data set (tests/crplib/reference/<case>/, see tests/crplib/README.md):
#             the data set is regenerated with build/tests/crplib/crplib_make_dataset, then the
#             driver evaluates the point lists *.in of the case (written by the same tool).
#             cu111_srp48: the real data set potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v (point
#             lists only; the data set is not copied into the case)
#
# Usage: tools/fortran_reference/run_crp_reference.sh kernels|<case> [scratch_dir]
#   scratch_dir  where to build and run, a directory that does not exist yet (default: mktemp -d)
# Every driver run must end with CRP_DRIVER_DONE (a Fortran STOP inside the library exits with
# status 0); the script fails otherwise. glibc and gfortran versions are recorded in VERSIONS.
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh).
# Run time: seconds per case. Docs: tests/crplib/README.md, tools/fortran_reference/README.md.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

(($# >= 1 && $# <= 2)) || { echo "usage: run_crp_reference.sh kernels|<case> [scratch_dir]" >&2; exit 2; }
target="$1"
scratch_dir="$(new_scratch_dir "${2:-}")"
reference_dir="${FORTRAN_REPO_DIR}/tests/crplib/reference"

driver="$("${FORTRAN_REPO_DIR}/tools/fortran_reference/crp_driver/build.sh" "${scratch_dir}/driver" | tail -n 1)"

# run_driver <input> <output> <mode> [arguments...]
run_driver() {
  local input="$1" output="$2"
  shift 2
  run_clean "${driver}" "$@" < "${input}" > "${output}"
  tail -n 1 "${output}" | grep -q '^CRP_DRIVER_DONE$' ||
    refuse "the driver did not finish ($*); see ${output}"
}

record_versions() {
  { echo "gfortran $(gfortran -dumpfullversion)"; ldd --version | head -n 1; } > "$1/VERSIONS"
}

if [[ "${target}" == "kernels" ]]; then
  out="${reference_dir}/kernels"
  mkdir -p "${out}"
  python3 -I -B "${FORTRAN_REPO_DIR}/tools/fortran_reference/crp_driver/make_kernel_inputs.py" "${out}"
  run_driver /dev/null "${out}/tables.out" tables
  for input in "${out}"/*.in; do
    name="$(basename "${input}" .in)"
    run_driver "${input}" "${out}/${name}.out" "${name%%_*}"
  done
  record_versions "${out}"
  echo "* Kernel references written to ${out}"
  exit 0
fi

case_dir="${reference_dir}/${target}"
generator="${FORTRAN_REPO_DIR}/build/tests/crplib/crplib_make_dataset"
[[ -x "${generator}" ]] || refuse "${generator} not built; cmake --build build"
# Cases on a real data set of the repository: the generator writes only the point lists; the data
# set stays where it is (not copied into the case).
declare -A real_data_sets=([cu111_srp48]="potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v")
data_set="${case_dir}/dataset"
if [[ -n "${real_data_sets[${target}]:-}" ]]; then
  data_set="${FORTRAN_REPO_DIR}/${real_data_sets[${target}]}"
  "${generator}" "${target}" "${case_dir}" "${data_set}"
else
  "${generator}" "${target}" "${case_dir}"
fi
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
done
record_versions "${case_dir}"
echo "* References of ${target} written to ${case_dir}"

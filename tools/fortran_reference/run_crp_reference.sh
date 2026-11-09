#!/usr/bin/env bash
# Generates the Fortran reference data of the CRP library tests (tests/crplib/reference/).
#
# Builds the CRP reference driver (crp_driver/build.sh, linked against the objects of the Fortran
# reference build) in a scratch directory and runs it:
#   kernels   Fourier tables and spline kernels: inputs from crp_driver/make_kernel_inputs.py,
#             written with the driver output to tests/crplib/reference/kernels/
#   <case>    a synthetic data set (tests/crplib/reference/<case>/, see tests/crplib/README.md):
#             the data set is regenerated with build/tests/crplib/crplib_make_dataset, then the
#             driver evaluates the point lists *.in of the case (written by the same tool)
#
# Usage: tools/fortran_reference/run_crp_reference.sh kernels|<case> [scratch_dir]
# Every driver run must end with CRP_DRIVER_DONE (a Fortran STOP inside the library exits with
# status 0); the script fails otherwise. glibc and gfortran versions are recorded in VERSIONS.
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
target="${1:?usage: run_crp_reference.sh kernels|<case> [scratch_dir]}"
scratch_dir="${2:-$(mktemp -d)}"
reference_dir="${repo_dir}/tests/crplib/reference"

driver="$("${repo_dir}/tools/fortran_reference/crp_driver/build.sh" "${scratch_dir}/driver" | tail -n 1)"

# run_driver <input> <output> <mode> [arguments...]
run_driver() {
  local input="$1" output="$2"
  shift 2
  "${driver}" "$@" < "${input}" > "${output}"
  tail -n 1 "${output}" | grep -q '^CRP_DRIVER_DONE$' || {
    echo "error: the driver did not finish ($*); see ${output}" >&2; exit 1; }
}

record_versions() {
  { echo "gfortran $(gfortran -dumpfullversion)"; ldd --version | head -n 1; } > "$1/VERSIONS"
}

if [[ "${target}" == "kernels" ]]; then
  out="${reference_dir}/kernels"
  mkdir -p "${out}"
  python3 -I "${repo_dir}/tools/fortran_reference/crp_driver/make_kernel_inputs.py" "${out}"
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
generator="${repo_dir}/build/tests/crplib/crplib_make_dataset"
[[ -x "${generator}" ]] || { echo "error: ${generator} not built (cmake --build build)" >&2; exit 1; }
"${generator}" "${target}" "${case_dir}"
# The Fortran opens the data set with a path of at most 50 characters: run with a short relative one.
work="${scratch_dir}/run"
rm -rf "${work}" && mkdir -p "${work}" && cp -r "${case_dir}/dataset" "${work}/dataset"
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

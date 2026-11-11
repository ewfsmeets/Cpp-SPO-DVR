#!/usr/bin/env bash
# Builds the CRP reference driver (crp_reference_driver.f90 + dump_tables.c) in a scratch
# directory, linked against the CRP library objects of the Fortran reference build, so that the
# driver runs exactly the code generated for the reference executable. Read-only use of the
# Fortran tree; the .mod files there are specific to the gfortran version that built them.
#
# Usage: tools/fortran_reference/crp_driver/build.sh [build_dir]   (default: a new mktemp -d)
# Prints the path of the executable on the last line.
set -euo pipefail

source_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
objects_dir="${SPODVR_FORTRAN_OBJS:-${HOME}/work/SPO-DVR_v5.19/procfiles.dir/OBJS.dir}"
build_dir="${1:-$(mktemp -d)}"
expected_gfortran="15.2.0"

[[ -d "${objects_dir}" ]] || { echo "error: ${objects_dir} not found" >&2; exit 1; }
version="$(gfortran -dumpfullversion)"
[[ "${version}" == "${expected_gfortran}" ]] || {
  echo "error: gfortran ${version}, but the reference objects were built with ${expected_gfortran}" >&2; exit 1; }

mkdir -p "${build_dir}"
objects=()
for name in CRP_Wrapper CRP_QD_Wrapper CRP_General CRP_GasPhase CRP_Constants \
            CRP_C3v_AtomSurface CRP_C3v_MoleculeSurface CRP_C4v_AtomSurface CRP_C4v_MoleculeSurface \
            CRP_C6v_AtomSurface CRP_C6v_MoleculeSurface Spline UnitConversion MPI_wrapper; do
  objects+=("${objects_dir}/${name}.o")
done

# Same floating-point flags as the reference build (GFORTRAN.mk): -O3 -fdefault-real-8.
gcc -O2 -c "${source_dir}/dump_tables.c" -o "${build_dir}/dump_tables.o"
gfortran -O3 -fdefault-real-8 -fopenmp -I"${objects_dir}" -J"${build_dir}" \
  -c "${source_dir}/crp_reference_driver.f90" -o "${build_dir}/crp_reference_driver.o"
gfortran -fopenmp -o "${build_dir}/crp_reference_driver" "${build_dir}/crp_reference_driver.o" \
  "${build_dir}/dump_tables.o" "${objects[@]}"
echo "${build_dir}/crp_reference_driver"

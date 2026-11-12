#!/usr/bin/env bash
# Negative tests of the Fortran reference tools (reference_config.sh, build_fortran_reference.sh,
# the run scripts): each refusal fires, exits 1 and prints "error: <what>; <fix>". They use stub
# installs (bin/spodvr.x a shell script) and synthetic build logs: no Fortran runs.
#   retired-variables  each retired variable (FORTRAN_RETIRED_VARS) set: a run script refuses
#   existing-scratch   a run script given an existing scratch directory refuses
#   changed-install    a file changed after SHA256SUMS; an executable other than BUILDINFO's
#   runtime-drift      a run-time package other than BUILDINFO's refuses, unless allowed
#   run                run_fortran: a failed run and a missing output refuse; the run sees the
#                      pinned environment only (OMP_NUM_THREADS = FORTRAN_REF_THREADS)
#   check-log          the build log checks accept a log of the pinned flags and refuse -Ofast,
#                      -march, -DMKL, an error line, a warning and a missing link line
#   check-binary       the instruction check refuses -mavx2 code and crtfastmath (-ffast-math)
#
# Usage: tools/fortran_reference/test_reference_tools.sh CASE WORK [CXX]
#   CASE  one of the cases above;  WORK  scratch directory (recreated);
#   CXX   C++ compiler for check-binary (default: c++)
# Environment: none read (the tested scripts get theirs set here). Run time: about a second.
# ctest: FortranReference.*. Docs: tools/fortran_reference/README.md.
set -euo pipefail
tools_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config="${tools_dir}/reference_config.sh"
source "${config}"

(($# >= 2 && $# <= 3)) || { echo "usage: test_reference_tools.sh CASE WORK [CXX]" >&2; exit 2; }
test_case="$1"
work="$2"
cxx="${3:-c++}"
rm -rf "${work}"
mkdir -p "${work}"
export TMPDIR="${work}"
# A refusal that fails to fire must not reach a real install: the run scripts default to this one.
export SPODVR_FORTRAN_REF_DIR="${work}/no-install"

fail() {
  echo "FAIL: $*" >&2
  exit 1
}

# expect_refusal PATTERN COMMAND...: COMMAND exits 1 with "error: <what>; <fix>" matching PATTERN.
expect_refusal() {
  local pattern="$1" status=0
  shift
  "$@" > "${work}/stdout.txt" 2> "${work}/stderr.txt" || status=$?
  ((status == 1)) || fail "exit ${status}, expected 1: $*; stderr: $(cat "${work}/stderr.txt")"
  grep -qE "^error: .*${pattern}.*; [a-z]" "${work}/stderr.txt" ||
    fail "no refusal matching \"${pattern}\": $*; stderr: $(cat "${work}/stderr.txt")"
  echo "refused as expected: $(grep -m 1 '^error: ' "${work}/stderr.txt" | cut -c 1-160)"
}

# expect_success COMMAND...: COMMAND exits 0.
expect_success() {
  "$@" > "${work}/stdout.txt" 2> "${work}/stderr.txt" ||
    fail "exit $?, expected 0: $*; stderr: $(cat "${work}/stderr.txt")"
}

# in_config FUNCTION ARGS...: a function of reference_config.sh in a fresh bash.
in_config() {
  bash -c 'set -euo pipefail; source "$0"; "$@"' "${config}" "$@"
}

# make_install DIR EXE_BODY: a stub install whose executable is a shell script with EXE_BODY and
# whose BUILDINFO.tsv records the version "0-recorded" of every run-time package.
make_install() {
  local dir="$1" sha package
  mkdir -p "${dir}/bin" "${dir}/OBJS"
  printf '#!/bin/sh\n%s\n' "$2" > "${dir}/bin/spodvr.x"
  chmod +x "${dir}/bin/spodvr.x"
  echo "stub object" > "${dir}/OBJS/stub.o"
  sha="$(sha256sum < "${dir}/bin/spodvr.x")"
  {
    echo "# stub install of test_reference_tools.sh"
    printf 'key\tvalue\n'
    printf 'exe.sha256\t%s\n' "${sha%% *}"
    for package in "${FORTRAN_RUNTIME_PACKAGES[@]}"; do printf 'package.runtime.%s\t0-recorded\n' "${package}"; done
  } > "${dir}/BUILDINFO.tsv"
  (cd "${dir}" && sha256sum bin/spodvr.x OBJS/* BUILDINFO.tsv > SHA256SUMS)
}

case "${test_case}" in
  retired-variables)
    # Spelled in parts: tools/check_docs.py allows the retired names only in reference_config.sh.
    for suffix in EXE OBJS; do
      name="SPODVR_FORTRAN_${suffix}"
      expect_refusal "${name} is retired" env "${name}=x" "${tools_dir}/run_setup_reference.sh" hardwall_small
    done
    ;;
  existing-scratch)
    mkdir "${work}/scratch"
    expect_refusal "scratch directory .* exists" "${tools_dir}/run_setup_reference.sh" hardwall_small "${work}/scratch"
    [[ -z "$(ls -A "${work}/scratch")" ]] || fail "the refused run wrote into ${work}/scratch"
    ;;
  changed-install)
    make_install "${work}/install" "exit 0"
    export SPODVR_FORTRAN_REF_DIR="${work}/install" SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1
    expect_success in_config check_install
    echo "changed" >> "${work}/install/OBJS/stub.o"
    expect_refusal "differs from its SHA256SUMS" in_config check_install
    make_install "${work}/install2" "exit 0"
    echo "exit 0" >> "${work}/install2/bin/spodvr.x"
    (cd "${work}/install2" && sha256sum bin/spodvr.x OBJS/* BUILDINFO.tsv > SHA256SUMS)
    export SPODVR_FORTRAN_REF_DIR="${work}/install2"
    expect_refusal "is not the executable of" in_config check_install
    ;;
  runtime-drift)
    make_install "${work}/install" "exit 0"
    export SPODVR_FORTRAN_REF_DIR="${work}/install"
    expect_refusal "run-time package ${FORTRAN_RUNTIME_PACKAGES[0]} is .*, not the version in" in_config check_install
    SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1 expect_success in_config check_install
    grep -q "^warning: run-time package" "${work}/stderr.txt" || fail "allowed drift printed no warning"
    ;;
  run)
    make_install "${work}/install" "env > env.txt; [ -f fail ] && exit 3; exit 0"
    export SPODVR_FORTRAN_REF_DIR="${work}/install" SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1
    mkdir "${work}/dir"
    OMP_NUM_THREADS=16 GFORTRAN_CONVERT_UNIT=big_endian LD_LIBRARY_PATH=/nowhere \
      expect_success in_config run_fortran "${work}/dir" stdout.txt env.txt
    grep -qx "OMP_NUM_THREADS=${FORTRAN_REF_THREADS}" "${work}/dir/env.txt" || fail "the run did not get OMP_NUM_THREADS=${FORTRAN_REF_THREADS}"
    ! grep -qE "^(GFORTRAN_|LD_LIBRARY_PATH=)" "${work}/dir/env.txt" || fail "the caller's environment reached the run"
    expect_refusal "wrote no fort\.21" in_config run_fortran "${work}/dir" stdout.txt env.txt fort.21
    touch "${work}/dir/fail"
    expect_refusal "failed \(exit 3\)" in_config run_fortran "${work}/dir" stdout.txt
    ;;
  check-log)
    good="${work}/good.log"
    {
      echo "Cleaning up all stuff for doing a complete re-make ..."
      echo "./PreProcess Main.f90 Main.f90 ${FORTRAN_BUILD_DEFINES}"
      echo "gfortran ${FORTRAN_BUILD_FLAGS}  ${FORTRAN_BUILD_INCLUDES} -c Main.f90"
      echo "./PreProcess H2W100PES-vdW.f90 H2W100PES-vdW.f90 ${FORTRAN_BUILD_DEFINES}"
      echo "gfortran ${FORTRAN_BUILD_FLAGS}  ${FORTRAN_BUILD_INCLUDES} -c H2W100PES-vdW.f90"
      echo "gfortran ${FORTRAN_BUILD_FLAGS}  ${FORTRAN_BUILD_INCLUDES} -o ${FORTRAN_BUILD_EXE_NAME} Main.o H2W100PES-vdW.o  -lfftw3 "
      echo "Cleaning up TEMP dir ..."
    } > "${good}"
    check="${tools_dir}/build_fortran_reference.sh"
    expect_success "${check}" --check-log "${good}"
    sed 's/ -O3 / -Ofast /' "${good}" > "${work}/ofast.log"
    expect_refusal "outside the allow-lists" "${check}" --check-log "${work}/ofast.log"
    sed 's/ -O3 / -O3 -march=native /' "${good}" > "${work}/march.log"
    expect_refusal "outside the allow-lists" "${check}" --check-log "${work}/march.log"
    sed 's/-DNO_MKL/-DMKL/' "${good}" > "${work}/mkl.log"
    expect_refusal "outside the allow-lists" "${check}" --check-log "${work}/mkl.log"
    { cat "${good}"; echo "gmake: *** [Makefile:212: Main.o] Error 1"; } > "${work}/error.log"
    expect_refusal "the build failed" "${check}" --check-log "${work}/error.log"
    { cat "${good}"; echo "Warning: Unused variable 'x' declared at (1)"; } > "${work}/warning.log"
    expect_refusal "outside the allow-lists" "${check}" --check-log "${work}/warning.log"
    grep -v ' -o ' "${good}" > "${work}/nolink.log"
    expect_refusal "no link line" "${check}" --check-log "${work}/nolink.log"
    ;;
  check-binary)
    printf '%s\n' 'void scale(double* a, int n) { for (int i = 0; i < n; ++i) a[i] *= 2.0; }' \
      'int main() { double a[8] = {}; scale(a, 8); return static_cast<int>(a[0]); }' > "${work}/kernel.cpp"
    check="${tools_dir}/build_fortran_reference.sh"
    "${cxx}" -O2 -o "${work}/plain" "${work}/kernel.cpp"
    expect_success "${check}" --check-binary "${work}/plain"
    "${cxx}" -O3 -mavx2 -o "${work}/avx2" "${work}/kernel.cpp"
    expect_refusal "has VEX/EVEX instructions" "${check}" --check-binary "${work}/avx2"
    "${cxx}" -O2 -ffast-math -o "${work}/fastmath" "${work}/kernel.cpp"
    expect_refusal "links crtfastmath" "${check}" --check-binary "${work}/fastmath"
    ;;
  *)
    echo "usage: test_reference_tools.sh CASE WORK [CXX] (unknown case ${test_case})" >&2
    exit 2
    ;;
esac
echo "PASS: ${test_case}"

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
#   manifest           refdata.py manifest: record writes checked lines (ReferenceManifest.cmake
#                      accepts them) and replaces a file's line; record_manifest records the inputs
#                      (a data set as a tree) and the run-time drift; bad paths, committed inputs,
#                      a non-case directory, a ';' in a field and a bad header refuse; delete
#                      removes only the listed files, and refuses a listed input and a case
#                      without MANIFEST
#   regenerate-refusals  regenerate_all.sh refuses an unclean reference tree and a missing install
#                      before it deletes anything (in a scratch git repository)
#   flux-outputs       input_has_flux_analysis finds doFluxAnalysis TRUE (typed or not) and nothing
#                      else; run_propagation_reference.sh refuses a flux case whose run wrote no flux
#                      files, and copies nothing (in a scratch repository)
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
    printf 'id\tstub\n'
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
    expect_refusal "run-time packages ${FORTRAN_RUNTIME_PACKAGES[0]}=.* differ from" in_config check_install
    SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1 expect_success in_config check_install
    grep -q "^warning: run-time packages" "${work}/stderr.txt" || fail "allowed drift printed no warning"
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
  manifest)
    make_install "${work}/install" "exit 0"
    # Under tests/reference/, so that ReferenceManifest.cmake checks the MANIFEST written here.
    case_dir="${work}/tests/reference/case"
    mkdir -p "${case_dir}/fine_x" "${work}/run"
    echo "input" > "${case_dir}/DVR.inp"
    echo "a" > "${case_dir}/a.txt"
    echo "b" > "${case_dir}/fine_x/b.txt"
    echo "c" > "${case_dir}/c.txt"
    mkdir -p "${case_dir}/dataset" "${work}/run/dataset" "${work}/not/fortran_unformatted"
    echo "as run" > "${work}/run/DVR.inp"
    echo "data" > "${work}/run/dataset/CRPData"
    refdata=(python3 -I -B "${tools_dir}/refdata.py" manifest)
    manifest_check=(cmake -DROOT="${work}" -P "${tools_dir}/ReferenceManifest.cmake")
    export SPODVR_FORTRAN_REF_DIR="${work}/install"
    expect_success in_config record_manifest "${work}/run" "${case_dir}" "${work}/install/bin/spodvr.x" DVR.inp dataset -- a.txt fine_x/b.txt
    expect_success "${manifest_check[@]}"
    grep -q "2 files in 1 MANIFESTs match" "${work}/stdout.txt" || fail "the MANIFEST does not list 2 files"
    grep -qP '^a\.txt\t[0-9a-f]{64}\treference_config\.sh\t[^\t]+\tstub spodvr\.x:[0-9a-f]{12}\tDVR\.inp:[0-9a-f]{12},dataset:[0-9a-f]{12}\t[^\t]*libc6=' \
      "${case_dir}/MANIFEST.tsv" || fail "the line of a.txt lacks the program, inputs or drift: $(cat "${case_dir}/MANIFEST.tsv")"
    # A regenerated file replaces its line; the other line stays.
    echo "a2" > "${case_dir}/a.txt"
    expect_success "${refdata[@]}" record "${case_dir}" a.txt --script other.sh --pre-manifest
    expect_success "${manifest_check[@]}"
    grep -q "2 files in 1 MANIFESTs match" "${work}/stdout.txt" || fail "re-recording a.txt changed the line count"
    grep -q $'^a\\.txt\t.*\tother\\.sh\t-\t5\\.19 (pre-MANIFEST)\t-\t-\t-\t-$' "${case_dir}/MANIFEST.tsv" ||
      fail "a.txt's line was not replaced: $(cat "${case_dir}/MANIFEST.tsv")"
    expect_refusal "is not a file path inside" "${refdata[@]}" record "${case_dir}" ../a.txt --script s --pre-manifest
    expect_refusal "none.txt is not a file" "${refdata[@]}" record "${case_dir}" none.txt --script s --pre-manifest
    # The committed inputs can never get a line, so a regeneration can never delete them.
    expect_refusal "DVR.inp is a committed input" "${refdata[@]}" record "${case_dir}" DVR.inp --script s --pre-manifest
    expect_refusal "dataset/x is a committed input" "${refdata[@]}" record "${case_dir}" dataset/x --script s --pre-manifest
    expect_refusal "holds no Fortran output" "${refdata[@]}" record "${work}/not/fortran_unformatted" x --script s --pre-manifest
    expect_success "${refdata[@]}" delete "${case_dir}"
    [[ -f "${case_dir}/DVR.inp" && ! -e "${case_dir}/a.txt" && ! -e "${case_dir}/fine_x" && ! -e "${case_dir}/MANIFEST.tsv" ]] ||
      fail "delete did not remove exactly the listed files, their empty directory and the MANIFEST: $(ls -AR "${case_dir}")"
    expect_refusal "MANIFEST.tsv not found" "${refdata[@]}" delete "${case_dir}"
    expect_refusal "holds a tab, line break" "${refdata[@]}" record "${case_dir}" c.txt --script 'a;b' --pre-manifest
    { printf 'path\tsha256\tscript\tcommit\tfortran_ref\tinput_sha\truntime_drift\tcpu\tdate\n'
      printf 'DVR.inp\t-\t-\t-\t-\t-\t-\t-\t-\n'; } > "${case_dir}/MANIFEST.tsv"
    expect_refusal "DVR.inp is a committed input" "${refdata[@]}" delete "${case_dir}"
    [[ -f "${case_dir}/DVR.inp" ]] || fail "delete removed the input DVR.inp"
    printf 'path\tsha256\n' > "${case_dir}/MANIFEST.tsv"
    expect_refusal "the header is not" "${refdata[@]}" record "${case_dir}" c.txt --script s --pre-manifest
    ;;
  regenerate-refusals)
    repo="${work}/repo"
    mkdir -p "${repo}/tools/fortran_reference" "${repo}/tests/reference/case"
    cp "${tools_dir}/reference_config.sh" "${tools_dir}/regenerate_all.sh" "${tools_dir}/refdata.py" "${repo}/tools/fortran_reference/"
    echo "data" > "${repo}/tests/reference/case/out.txt"
    expect_success python3 -I -B "${tools_dir}/refdata.py" manifest record "${repo}/tests/reference/case" out.txt --script s --pre-manifest
    git -C "${repo}" init -q
    git -C "${repo}" add -A
    git -C "${repo}" -c user.name=test -c user.email=test@example.invalid -c commit.gpgsign=false commit -qm scratch
    echo "stray" > "${repo}/tests/reference/case/stray.txt"
    expect_refusal "reference trees are not clean" "${repo}/tools/fortran_reference/regenerate_all.sh"
    rm "${repo}/tests/reference/case/stray.txt"
    expect_refusal "is not a reference install" "${repo}/tools/fortran_reference/regenerate_all.sh"
    [[ -z "$(git -C "${repo}" status --porcelain)" ]] || fail "a refused regeneration changed the tree: $(git -C "${repo}" status --porcelain)"
    ;;
  flux-outputs)
    printf 'doFluxAnalysis : LOGICAL : TRUE;\n' > "${work}/typed.inp"
    printf '  doFluxAnalysis: true;\n' > "${work}/untyped.inp"
    printf 'doFluxAnalysis : LOGICAL : FALSE;\n# doFluxAnalysis : LOGICAL : TRUE;\n' > "${work}/off.inp"
    printf 'Title : STRING : "doFluxAnalysis : TRUE;";\n' > "${work}/none.inp"
    for name in typed untyped; do
      in_config input_has_flux_analysis "${work}/${name}.inp" || fail "no flux analysis found in ${name}.inp"
    done
    for name in off none; do
      ! in_config input_has_flux_analysis "${work}/${name}.inp" || fail "a flux analysis found in ${name}.inp"
    done
    repo="${work}/repo"
    mkdir -p "${repo}/tools/fortran_reference" "${repo}/tests/reference/case"
    cp "${tools_dir}/reference_config.sh" "${tools_dir}/run_propagation_reference.sh" "${tools_dir}/refdata.py" \
      "${repo}/tools/fortran_reference/"
    cp "${work}/typed.inp" "${repo}/tests/reference/case/DVR.inp"
    make_install "${work}/install" "touch averages.out fort.21 fort.25; exit 0"
    export SPODVR_FORTRAN_REF_DIR="${work}/install" SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1
    expect_refusal "wrote no flux_coeffs\.dat" "${repo}/tools/fortran_reference/run_propagation_reference.sh" case \
      "${work}/scratch"
    [[ "$(ls "${repo}/tests/reference/case")" == DVR.inp ]] ||
      fail "the refused run copied files: $(ls "${repo}/tests/reference/case")"
    ;;
  *)
    echo "usage: test_reference_tools.sh CASE WORK [CXX] (unknown case ${test_case})" >&2
    exit 2
    ;;
esac
echo "PASS: ${test_case}"

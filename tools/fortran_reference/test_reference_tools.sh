#!/usr/bin/env bash
# Tests of the Fortran reference tools (reference_config.sh, build_fortran_reference.sh, the run
# scripts): each refusal fires, exits 1 and prints "error: <what>; <fix>"; the run scripts stage a
# case's inputs as the Fortran reads them. They use stub installs (bin/spodvr.x a shell script) and
# synthetic build logs: no Fortran runs.
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
#   stage-scm          run_propagation_reference.sh on an SCM case whose three data files have paths
#                      of over 100 characters (in a scratch git repository): the run finds them as
#                      dataset, params.dat and slab.dat, the DVR.inp as run names them so and is
#                      otherwise the case's, the MANIFEST and staged_inputs (with a prefix) list them
#                      as inputs; a missing data file and a DVR.inp without the field refuse
#   staged-cases       stage_inputs on every case of tests/reference without SCM keys: the DVR.inp
#                      byte for byte, the CRP data set unchanged, nothing else staged, and
#                      staged_inputs lists exactly what was staged
#
# Usage: tools/fortran_reference/test_reference_tools.sh CASE WORK [CXX]
#   CASE  one of the cases above;  WORK  scratch directory (recreated);
#   CXX   C++ compiler for check-binary (default: c++)
# Environment: none read (the tested scripts get theirs set here). Run time: about a second each.
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

# json_value FILE KEY: the value of KEY in the JSON file FILE ("" when absent).
json_value() {
  python3 -I -B -c 'import json, sys; print(json.load(open(sys.argv[1])).get(sys.argv[2], ""))' "$1" "$2"
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
  stage-scm)
    repo="${work}/repo"
    mkdir -p "${repo}/tools/fortran_reference"
    cp "${config}" "${tools_dir}/refdata.py" "${tools_dir}/run_propagation_reference.sh" "${repo}/tools/fortran_reference/"
    data="potentials/$(printf 'a_directory_of_forty_characters_number_%d/' 1 2 3)"
    ((${#data} > 100)) || fail "the data directory ${data} has ${#data} characters, not over 100"
    case_dir="${repo}/tests/reference/scm_case"
    mkdir -p "${repo}/${data}dataset" "${case_dir}"
    echo "CRP data" > "${repo}/${data}dataset/CRPData"
    echo "SCM parameters" > "${repo}/${data}params.dat"
    echo "slab positions" > "${repo}/${data}positions.dat"
    printf '%s\n' '# An SCM case' 'SystemID : STRING : "CRPDistrtSurface";' \
      "CRPlibDataSet : STRING : \"${data}dataset\";" "SCMParFilepath : STRING : \"${data}params.dat\";" \
      "PosFilepath : STRING : \"${data}positions.dat\";" 'PosFileSlabNum : INTEGER : 1;' > "${case_dir}/DVR.inp"
    printf '{"SystemID": "CRPDistrtSurface", "CRPlibDataSet": "%s", "SCMParFilepath": "%s", "PosFilepath": "%s"}\n' \
      "${data}dataset" "${data}params.dat" "${data}positions.dat" > "${case_dir}/input.json"
    git -C "${repo}" init -q
    git -C "${repo}" add -A
    git -C "${repo}" -c user.name=test -c user.email=test@example.invalid -c commit.gpgsign=false commit -qm scratch
    # The stub runs only where it finds the staged files under their short names.
    make_install "${work}/install" "[ -f dataset/CRPData ] && [ -f params.dat ] && [ -f slab.dat ] || exit 3
touch averages.out fort.21 fort.25"
    export SPODVR_FORTRAN_REF_DIR="${work}/install" SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1
    run="${work}/run"
    expect_success "${repo}/tools/fortran_reference/run_propagation_reference.sh" scm_case "${run}"
    cmp -s "${repo}/${data}params.dat" "${run}/params.dat" && cmp -s "${repo}/${data}positions.dat" "${run}/slab.dat" &&
      diff -r "${repo}/${data}dataset" "${run}/dataset" > /dev/null || fail "the staged data files differ from the case's"
    for line in 'CRPlibDataSet : STRING : "dataset";' 'SCMParFilepath : STRING : "params.dat";' \
      'PosFilepath : STRING : "slab.dat";'; do
      grep -qxF "${line}" "${run}/DVR.inp" || fail "the DVR.inp as run has no line ${line}: $(cat "${run}/DVR.inp")"
    done
    paths='^(CRPlibDataSet|SCMParFilepath|PosFilepath) :'
    diff <(grep -vE "${paths}" "${case_dir}/DVR.inp") <(grep -vE "${paths}" "${run}/DVR.inp") > /dev/null ||
      fail "the DVR.inp as run differs from the case's in other lines"
    grep -qP '^fortran_stdout\.txt\t.*\tDVR\.inp:[0-9a-f]{12},dataset:[0-9a-f]{12},params\.dat:[0-9a-f]{12},slab\.dat:[0-9a-f]{12}\t' \
      "${case_dir}/MANIFEST.tsv" || fail "the MANIFEST does not list the staged inputs: $(cat "${case_dir}/MANIFEST.tsv")"
    in_repo() { bash -c 'set -euo pipefail; source "$0"; "$@"' "${repo}/tools/fortran_reference/reference_config.sh" "$@"; }
    # The input list with a prefix, as the restart script records its two runs.
    listed="$(in_repo staged_inputs "${run}" part1/ | paste -sd ' ')"
    [[ "${listed}" == "part1/DVR.inp part1/dataset part1/params.dat part1/slab.dat" ]] ||
      fail "staged_inputs with a prefix lists ${listed}"
    # stage_inputs alone: a data file missing, then a DVR.inp without a field its input.json names.
    mkdir "${work}/stage1" "${work}/stage2"
    mv "${repo}/${data}positions.dat" "${work}/positions.dat"
    expect_refusal "PosFilepath .*positions\.dat .*not found" in_repo stage_inputs "${case_dir}" "${work}/stage1"
    mv "${work}/positions.dat" "${repo}/${data}positions.dat"
    sed -i '/^PosFilepath /d' "${case_dir}/DVR.inp"
    expect_refusal "field PosFilepath not found" in_repo stage_inputs "${case_dir}" "${work}/stage2"
    ;;
  staged-cases)
    count=0
    for case_dir in "${FORTRAN_REPO_DIR}"/tests/reference/*; do
      [[ -f "${case_dir}/DVR.inp" ]] || continue
      [[ -z "$(json_value "${case_dir}/input.json" SCMParFilepath)" ]] || continue
      name="$(basename "${case_dir}")"
      mkdir "${work}/${name}"
      expect_success in_config stage_inputs "${case_dir}" "${work}/${name}"
      cmp -s "${case_dir}/DVR.inp" "${work}/${name}/DVR.inp" || fail "${name}: the staged DVR.inp differs from the case's"
      data_set="$(json_value "${case_dir}/input.json" CRPlibDataSet)"
      staged="DVR.inp"
      if [[ -n "${data_set}" ]]; then
        diff -r "${FORTRAN_REPO_DIR}/${data_set}" "${work}/${name}/dataset" > /dev/null ||
          fail "${name}: the staged data set differs from ${data_set}"
        staged="DVR.inp dataset"
      fi
      listing=("${work}/${name}"/*)
      [[ "${listing[*]##*/}" == "${staged}" ]] || fail "${name}: staged ${listing[*]##*/}, expected ${staged}"
      listed="$(in_config staged_inputs "${work}/${name}" | paste -sd ' ')"
      [[ "${listed}" == "${staged}" ]] || fail "${name}: staged_inputs lists ${listed}, expected ${staged}"
      count=$((count + 1))
    done
    ((count >= 11)) || fail "only ${count} cases staged; tests/reference has at least 11"
    echo "${count} cases staged as before"
    ;;
  *)
    echo "usage: test_reference_tools.sh CASE WORK [CXX] (unknown case ${test_case})" >&2
    exit 2
    ;;
esac
echo "PASS: ${test_case}"

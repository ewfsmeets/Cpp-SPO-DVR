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
#                      -march, -DMKL, an error line, a warning, a missing link line and no log
#   check-binary       the instruction check refuses -mavx2 code, crtfastmath (-ffast-math), an
#                      AVX libmvec import and a file that is not an object
#   config-refusals    the other refusals of reference_config.sh: a BUILDINFO.tsv without header, a
#                      build package other than BUILDINFO's, record_manifest or build_fortran_driver
#                      without "--", a driver source neither .f90 nor .c, a scratch directory
#                      without parent, a case without DVR.inp, a DVR.inp field without type
#   build-refusals     build_fortran_reference.sh: an existing install, a mirror without the pinned
#                      commit; on a scratch mirror pinned in a copy of reference_config.sh: another
#                      tree, no svn revision, another GNU.mk line 19, a failed DoBuild, no executable
#   manifest           refdata.py manifest: record writes checked lines (ReferenceManifest.cmake
#                      accepts them) and replaces a file's line; record_manifest records the inputs
#                      (a data set as a tree) and the run-time drift; bad paths, committed inputs,
#                      a non-case directory, a ';' in a field and a bad header refuse; delete
#                      removes only the listed files, and refuses a listed input and a case
#                      without MANIFEST; a missing case directory, an install without BUILDINFO.tsv
#                      or its id, a failing git, a short row and a MANIFEST without header refuse
#   regenerate-refusals  regenerate_all.sh refuses an unclean reference tree and a missing install
#                      before it deletes anything (in a scratch git repository); with a stub install,
#                      a step that writes nothing and a tree without MANIFESTs
#   run-script-refusals  run_fine_analysis/crp/scm_reference.sh in a scratch repository with stub
#                      driver builds: each missing input, a driver that fails or does not finish, an
#                      output git does not track (CRP, SCM); scm_driver/make_scm_inputs.py: a
#                      params.dat of 35 values; run_pes_cut_reference.sh: no pes_cuts.json, a file
#                      that is not JSON, no cut, a key that is not a coordinate, a value that is
#                      a string, a bool or not finite, a pes_cut file left from more cuts, a DVR.inp
#                      that sets a write2DPES key (each before the scratch directory is made)
#   unmodified-targets require_unmodified in a scratch git repository passes an unchanged, an absent
#                      and a deleted target, and refuses an edited file, a directory with an edited
#                      file, a file added to the index, an untracked file, a directory with one, a
#                      path outside the repository and a repository git cannot read;
#                      run_propagation_reference.sh and run_fine_analysis_reference.sh refuse such a
#                      target before they run; the propagation runs once its targets are deleted
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
# shellcheck source=reference_config.sh
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

# in_config_of CONFIG FUNCTION ARGS...: a function of the copy CONFIG of reference_config.sh in a
# fresh bash; in_config FUNCTION ARGS...: the same with this directory's reference_config.sh.
in_config_of() {
  bash -c 'set -euo pipefail; source "$0"; "$@"' "$@"
}
in_config() {
  in_config_of "${config}" "$@"
}

# make_install DIR EXE_BODY [ROW...]: a stub install whose executable is a shell script with
# EXE_BODY and whose BUILDINFO.tsv records the version "0-recorded" of every run-time package, and
# the ROWs ("key<TAB>value").
make_install() {
  local dir="$1" body="$2" sha package
  shift 2
  mkdir -p "${dir}/bin" "${dir}/OBJS"
  printf '#!/bin/sh\n%s\n' "${body}" > "${dir}/bin/spodvr.x"
  chmod +x "${dir}/bin/spodvr.x"
  echo "stub object" > "${dir}/OBJS/stub.o"
  sha="$(sha256sum < "${dir}/bin/spodvr.x")"
  {
    echo "# stub install of test_reference_tools.sh"
    printf 'key\tvalue\n'
    printf 'exe.sha256\t%s\n' "${sha%% *}"
    printf 'id\tstub\n'
    for package in "${FORTRAN_RUNTIME_PACKAGES[@]}"; do printf 'package.runtime.%s\t0-recorded\n' "${package}"; done
    (($# == 0)) || printf '%s\n' "$@"
  } > "${dir}/BUILDINFO.tsv"
  (cd "${dir}" && sha256sum bin/spodvr.x OBJS/* BUILDINFO.tsv > SHA256SUMS)
}

# scratch_repo DIR FILE...: DIR/tools/fortran_reference/ with copies of reference_config.sh,
# refdata.py and the named FILEs of this directory, so that DIR is their repository: a git
# repository with these files committed (the run scripts check their targets with git).
scratch_repo() {
  local dir="$1" file
  shift
  mkdir -p "${dir}/tools/fortran_reference"
  for file in reference_config.sh refdata.py "$@"; do cp "${tools_dir}/${file}" "${dir}/tools/fortran_reference/"; done
  git init -q "${dir}"
  commit_all "${dir}" "the tools"
}

# stub_driver_build DIR STATUS: DIR/build.sh, a stand-in for a driver build script that prints the
# path of a driver which exits with STATUS without any output (no end marker).
stub_driver_build() {
  mkdir -p "$1"
  # shellcheck disable=SC2016  # $1 belongs to the stub script being written
  printf '#!/bin/sh\nmkdir -p "$1" && printf "#!/bin/sh\\nexit %s\\n" > "$1/driver" && chmod +x "$1/driver" && echo "$1/driver"\n' \
    "$2" > "$1/build.sh"
  chmod +x "$1/build.sh"
}

# good_log: a DoBuild log of the pinned flags that passes the log checks.
good_log() {
  echo "Cleaning up all stuff for doing a complete re-make ..."
  echo "./PreProcess Main.f90 Main.f90 ${FORTRAN_BUILD_DEFINES}"
  echo "gfortran ${FORTRAN_BUILD_FLAGS}  ${FORTRAN_BUILD_INCLUDES} -c Main.f90"
  echo "./PreProcess H2W100PES-vdW.f90 H2W100PES-vdW.f90 ${FORTRAN_BUILD_DEFINES}"
  echo "gfortran ${FORTRAN_BUILD_FLAGS}  ${FORTRAN_BUILD_INCLUDES} -c H2W100PES-vdW.f90"
  echo "gfortran ${FORTRAN_BUILD_FLAGS}  ${FORTRAN_BUILD_INCLUDES} -o ${FORTRAN_BUILD_EXE_NAME} Main.o H2W100PES-vdW.o  -lfftw3 "
  echo "Cleaning up TEMP dir ..."
}

# commit_all REPO MESSAGE: commits every file of the scratch repository REPO.
commit_all() {
  git -C "$1" add -A
  git -C "$1" -c user.name=test -c user.email=test@example.invalid -c commit.gpgsign=false commit -q --allow-empty -m "$2"
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
    # A case that does not exist, so that only the scratch directory can be refused.
    expect_refusal "scratch directory .* exists" "${tools_dir}/run_setup_reference.sh" no_case "${work}/scratch"
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
    good_log > "${good}"
    check="${tools_dir}/build_fortran_reference.sh"
    expect_success "${check}" --check-log "${good}"
    expect_refusal "build log .*none\.log not found" "${check}" --check-log "${work}/none.log"
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
    # An object calling the AVX2 libmvec exp without VEX code of its own (plain -O2).
    printf '%s\n' 'extern "C" double _ZGVdN4v_exp(double);' 'double call(double x) { return _ZGVdN4v_exp(x); }' \
      > "${work}/libmvec.cpp"
    "${cxx}" -O2 -c -o "${work}/libmvec.o" "${work}/libmvec.cpp"
    expect_refusal "imports AVX libmvec functions _ZGVdN4v_exp" "${check}" --check-binary "${work}/libmvec.o"
    expect_refusal "kernel\.cpp is not an object file" "${check}" --check-binary "${work}/kernel.cpp"
    ;;
  config-refusals)
    make_install "${work}/install" "exit 0" $'package.build.gfortran-15\t0-recorded'
    export SPODVR_FORTRAN_REF_DIR="${work}/install" SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1
    expect_refusal "build package gfortran-15 is .*, not the version in" in_config check_build_package gfortran-15
    expect_refusal "record_manifest: no" in_config record_manifest "${work}" "${work}" program DVR.inp
    expect_refusal "cannot create scratch directory" in_config new_scratch_dir "${work}/missing/scratch"
    mkdir "${work}/no_case" "${work}/stage"
    expect_refusal "no_case/DVR\.inp not found" in_config stage_inputs "${work}/no_case" "${work}/stage"
    echo "tEnd : 7000.0;" > "${work}/DVR.inp"
    expect_refusal "field tEnd of .* is not" in_config set_field "${work}/DVR.inp" tEnd 1.0
    make_install "${work}/install2" "exit 0" $'package.build.gfortran-15\t-'
    export SPODVR_FORTRAN_REF_DIR="${work}/install2"
    expect_refusal "driver source .*DVR\.inp is neither" in_config build_fortran_driver "${work}/d1" name "${work}/DVR.inp" -- o
    expect_refusal "before the object list" in_config build_fortran_driver "${work}/d2" name
    sed -i '/^key\tvalue$/d' "${work}/install2/BUILDINFO.tsv"
    expect_refusal "has no header row" in_config check_install
    ;;
  build-refusals)
    check="${tools_dir}/build_fortran_reference.sh"
    mkdir "${work}/existing"
    SPODVR_FORTRAN_REF_DIR="${work}/existing" expect_refusal "existing exists" "${check}"
    mirror="${work}/mirror"
    git init -q "${mirror}"
    SPODVR_FORTRAN_MIRROR="${mirror}" expect_refusal "commit ${FORTRAN_BUILD_COMMIT} is not in" "${check}"
    # A copy of the build script whose pin names commits of the scratch mirror (pin_mirror).
    scratch_repo "${work}/repo" build_fortran_reference.sh
    check="${work}/repo/tools/fortran_reference/build_fortran_reference.sh"
    export SPODVR_FORTRAN_MIRROR="${mirror}"
    pin_mirror() {  # MESSAGE [TREE]: commits the mirror, pins the copy to it (and TREE, default its own)
      commit_all "${mirror}" "$1"
      local commit tree
      commit="$(git -C "${mirror}" rev-parse HEAD)"
      tree="${2:-$(git -C "${mirror}" rev-parse "HEAD^{tree}")}"
      sed -i -E -e "s/^FORTRAN_BUILD_COMMIT=.*/FORTRAN_BUILD_COMMIT=\"${commit}\"/" \
        -e "s/^FORTRAN_BUILD_TREE=.*/FORTRAN_BUILD_TREE=\"${tree}\"/" "${work}/repo/tools/fortran_reference/reference_config.sh"
      grep -qx "FORTRAN_BUILD_COMMIT=\"${commit}\"" "${work}/repo/tools/fortran_reference/reference_config.sh" ||
        fail "the pin of the copied reference_config.sh was not replaced"
    }
    mkdir "${mirror}/sources.dir"
    gnu_mk() { for line in {1..18}; do echo "# line ${line}"; done; echo "$1"; }
    gnu_mk "${FORTRAN_BUILD_GNU_MK_LINE}" > "${mirror}/sources.dir/GNU.mk"
    good_log > "${work}/good.log"
    printf '#!/bin/sh\ncat %s\n' "${work}/good.log" > "${mirror}/sources.dir/DoBuild"
    chmod +x "${mirror}/sources.dir/DoBuild"
    svn_message="git-svn-id: svn://test/trunk@${FORTRAN_BUILD_SVN#r} 00000000-0000-0000-0000-000000000000"
    pin_mirror "${svn_message}" 0000000000000000000000000000000000000000
    expect_refusal "has tree [0-9a-f]+, not 0+" "${check}"
    pin_mirror "a commit without svn revision"
    expect_refusal "is not svn trunk ${FORTRAN_BUILD_SVN}" "${check}"
    gnu_mk "O3FLAGS = -O2" > "${mirror}/sources.dir/GNU.mk"
    pin_mirror "${svn_message}"
    expect_refusal "line 19 of .*GNU\.mk is not" "${check}"
    gnu_mk "${FORTRAN_BUILD_GNU_MK_LINE}" > "${mirror}/sources.dir/GNU.mk"
    pin_mirror "${svn_message}"
    expect_refusal "DoBuild wrote no .*executable\.dir" "${check}"
    printf '#!/bin/sh\nexit 1\n' > "${mirror}/sources.dir/DoBuild"
    pin_mirror "${svn_message}"
    expect_refusal "DoBuild failed" "${check}"
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
    { printf 'path\tsha256\tscript\tcommit\tfortran_ref\tinput_sha\truntime_drift\tcpu\tdate\n'
      printf 'c.txt\t-\n'; } > "${case_dir}/MANIFEST.tsv"
    expect_refusal "MANIFEST\.tsv:[0-9]+: 2 fields, expected 9" "${refdata[@]}" delete "${case_dir}"
    echo "# a writer line only" > "${case_dir}/MANIFEST.tsv"
    expect_refusal "MANIFEST\.tsv has no header row" "${refdata[@]}" delete "${case_dir}"
    expect_refusal "case directory .*no-case not found" "${refdata[@]}" record "${work}/no-case" c.txt --script s --pre-manifest
    run_options=(--script s --install "${work}/other" --binary "${work}/install/bin/spodvr.x")
    mkdir "${work}/other"
    expect_refusal "other/BUILDINFO\.tsv not found" "${refdata[@]}" record "${case_dir}" c.txt "${run_options[@]}"
    printf 'key\tvalue\n' > "${work}/other/BUILDINFO.tsv"
    expect_refusal "other/BUILDINFO\.tsv has no id row" "${refdata[@]}" record "${case_dir}" c.txt "${run_options[@]}"
    GIT_DIR="${work}/no-git" expect_refusal "git .* failed in" "${refdata[@]}" record "${case_dir}" c.txt "${run_options[@]}"
    ;;
  regenerate-refusals)
    repo="${work}/repo"
    scratch_repo "${repo}" regenerate_all.sh
    mkdir -p "${repo}/tests/reference/case"
    echo "data" > "${repo}/tests/reference/case/out.txt"
    expect_success python3 -I -B "${tools_dir}/refdata.py" manifest record "${repo}/tests/reference/case" out.txt --script s --pre-manifest
    commit_all "${repo}" "a reference case"
    echo "stray" > "${repo}/tests/reference/case/stray.txt"
    expect_refusal "reference trees are not clean" "${repo}/tools/fortran_reference/regenerate_all.sh"
    rm "${repo}/tests/reference/case/stray.txt"
    expect_refusal "is not a reference install" "${repo}/tools/fortran_reference/regenerate_all.sh"
    [[ -z "$(git -C "${repo}" status --porcelain)" ]] || fail "a refused regeneration changed the tree: $(git -C "${repo}" status --porcelain)"
    # With an install that passes: a step that writes nothing (after the deletion), then no MANIFEST.
    make_install "${work}/install" "exit 0" $'package.build.gfortran-15\t-'
    export SPODVR_FORTRAN_REF_DIR="${work}/install" SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1
    printf '#!/bin/sh\nexit 0\n' > "${repo}/tools/fortran_reference/run_setup_reference.sh"
    chmod +x "${repo}/tools/fortran_reference/run_setup_reference.sh"
    expect_refusal "the step wrote no .*hardwall_small/fortran_setup_stdout\.txt" "${repo}/tools/fortran_reference/regenerate_all.sh"
    commit_all "${repo}" "the deleted references"
    expect_refusal "no MANIFEST\.tsv in" "${repo}/tools/fortran_reference/regenerate_all.sh"
    ;;
  stage-scm)
    repo="${work}/repo"
    scratch_repo "${repo}" run_propagation_reference.sh
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
    commit_all "${repo}" "an SCM case"
    # The stub runs only where it finds the staged files under their short names.
    make_install "${work}/install" "[ -f dataset/CRPData ] && [ -f params.dat ] && [ -f slab.dat ] || exit 3
touch averages.out fort.21 fort.25"
    export SPODVR_FORTRAN_REF_DIR="${work}/install" SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1
    run="${work}/run"
    expect_success "${repo}/tools/fortran_reference/run_propagation_reference.sh" scm_case "${run}"
    cmp -s "${repo}/${data}params.dat" "${run}/params.dat" || fail "the staged params.dat differs from the case's"
    cmp -s "${repo}/${data}positions.dat" "${run}/slab.dat" || fail "the staged slab.dat differs from the case's positions.dat"
    diff -r "${repo}/${data}dataset" "${run}/dataset" > /dev/null || fail "the staged dataset differs from the case's"
    for line in 'CRPlibDataSet : STRING : "dataset";' 'SCMParFilepath : STRING : "params.dat";' \
      'PosFilepath : STRING : "slab.dat";'; do
      grep -qxF "${line}" "${run}/DVR.inp" || fail "the DVR.inp as run has no line ${line}: $(cat "${run}/DVR.inp")"
    done
    paths='^(CRPlibDataSet|SCMParFilepath|PosFilepath) :'
    diff <(grep -vE "${paths}" "${case_dir}/DVR.inp") <(grep -vE "${paths}" "${run}/DVR.inp") > /dev/null ||
      fail "the DVR.inp as run differs from the case's in other lines"
    grep -qP '^fortran_stdout\.txt\t.*\tDVR\.inp:[0-9a-f]{12},dataset:[0-9a-f]{12},params\.dat:[0-9a-f]{12},slab\.dat:[0-9a-f]{12}\t' \
      "${case_dir}/MANIFEST.tsv" || fail "the MANIFEST does not list the staged inputs: $(cat "${case_dir}/MANIFEST.tsv")"
    repo_config="${repo}/tools/fortran_reference/reference_config.sh"
    # The input list with a prefix, as the restart script records its two runs.
    listed="$(in_config_of "${repo_config}" staged_inputs "${run}" part1/ | paste -sd ' ')"
    [[ "${listed}" == "part1/DVR.inp part1/dataset part1/params.dat part1/slab.dat" ]] ||
      fail "staged_inputs with a prefix lists ${listed}"
    # stage_inputs alone: a data file missing, then a DVR.inp without a field its input.json names.
    mkdir "${work}/stage1" "${work}/stage2"
    mv "${repo}/${data}positions.dat" "${work}/positions.dat"
    expect_refusal "PosFilepath .*positions\.dat .*not found" in_config_of "${repo_config}" stage_inputs "${case_dir}" "${work}/stage1"
    mv "${work}/positions.dat" "${repo}/${data}positions.dat"
    sed -i '/^PosFilepath /d' "${case_dir}/DVR.inp"
    expect_refusal "field PosFilepath not found" in_config_of "${repo_config}" stage_inputs "${case_dir}" "${work}/stage2"
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
  run-script-refusals)
    repo="${work}/repo"
    scratch_repo "${repo}" run_fine_analysis_reference.sh run_crp_reference.sh run_scm_reference.sh \
      run_pes_cut_reference.sh pes_cuts.py
    scripts="${repo}/tools/fortran_reference"
    expect_refusal "no_case/fortran_fort21\.bin not found" "${scripts}/run_fine_analysis_reference.sh" no_case full
    stub_driver_build "${scripts}/crp_driver" 0
    kernels="${repo}/tests/crplib/reference/kernels"
    mkdir -p "${kernels}"
    expect_refusal "no point lists .* in .*kernels" "${scripts}/run_crp_reference.sh" kernels
    touch "${kernels}/spline.in"
    expect_refusal "the driver did not finish \(tables\)" "${scripts}/run_crp_reference.sh" kernels
    # The failed run left an untracked tables.out: a rerun refuses before the driver runs.
    expect_refusal "kernels/tables\.out: not tracked by git" "${scripts}/run_crp_reference.sh" kernels
    mkdir -p "${repo}/tests/crplib/reference/c3v_case"
    touch "${repo}/tests/crplib/reference/c3v_case/molecule.in"
    expect_refusal "crplib_make_dataset not built" "${scripts}/run_crp_reference.sh" --generate c3v_case
    expect_refusal "data set .*c3v_case/dataset not found" "${scripts}/run_crp_reference.sh" c3v_case
    synthetic="${repo}/tests/scm/reference/synthetic"
    mkdir -p "${synthetic}"
    expect_refusal "no point list .*synthetic/points\.in" "${scripts}/run_scm_reference.sh" synthetic
    touch "${synthetic}/points.in"
    expect_refusal "data set .*H2Cu111SRP_C3v not found" "${scripts}/run_scm_reference.sh" synthetic
    mkdir -p "${repo}/potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v"
    touch "${repo}/potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v/CRPData"
    stub_driver_build "${scripts}/scm_driver" 3
    expect_refusal "input .*params\.in or .*slab_v3\.in not found" "${scripts}/run_scm_reference.sh" synthetic
    touch "${synthetic}/params.in" "${synthetic}/slab_v3.in"
    expect_refusal "the driver failed \(exit 3\)" "${scripts}/run_scm_reference.sh" synthetic
    stub_driver_build "${scripts}/scm_driver" 0
    expect_refusal "the driver did not finish in" "${scripts}/run_scm_reference.sh" synthetic
    echo "stray" > "${synthetic}/v2_slab2_rescaling_off.out"
    expect_refusal "synthetic/v2_slab2_rescaling_off\.out: not tracked by git" "${scripts}/run_scm_reference.sh" synthetic
    # The input generator of --generate: an example params.dat of 35 values (one missing).
    mkdir "${work}/example"
    { echo "VERSION 2"; printf '%s\n' {1..35}; } > "${work}/example/params.dat"
    expect_refusal "params\.dat has 35 values, not 36" \
      python3 -I -B "${tools_dir}/scm_driver/make_scm_inputs.py" example "${work}/scm_case" "${work}/example"
    [[ ! -e "${work}/scm_case" ]] || fail "the refused generator created ${work}/scm_case"
    # run_pes_cut_reference.sh: the cut definitions are checked before the scratch directory is made.
    cut_case="${repo}/tests/reference/cut_case"
    cut_scratch="${work}/cut_scratch"
    mkdir -p "${cut_case}"
    # cut_refusal PATTERN: the script refuses the case with PATTERN and makes no scratch directory.
    cut_refusal() {
      expect_refusal "$1" "${scripts}/run_pes_cut_reference.sh" cut_case "${cut_scratch}"
      [[ ! -e "${cut_scratch}" ]] || fail "the refused run made ${cut_scratch}"
    }
    cut_refusal "cut_case/pes_cuts\.json not found"
    echo '{"cuts": ' > "${cut_case}/pes_cuts.json"
    cut_refusal "pes_cuts\.json is not JSON with a list"
    echo '{"cuts": []}' > "${cut_case}/pes_cuts.json"
    cut_refusal "pes_cuts\.json lists no cut"
    echo '{"cuts": [{}, {"z2DPES": 1.0}]}' > "${cut_case}/pes_cuts.json"
    cut_refusal "cut 2: z2DPES is not a write2DPES coordinate"
    echo '{"cuts": [{"x2DPES": "1.0"}]}' > "${cut_case}/pes_cuts.json"
    cut_refusal "cut 1: x2DPES = '1\.0' is not a finite number"
    echo '{"cuts": [{"x2DPES": true}]}' > "${cut_case}/pes_cuts.json"
    cut_refusal "cut 1: x2DPES = True is not a finite number"
    echo '{"cuts": [{"phi2DPES": NaN}]}' > "${cut_case}/pes_cuts.json"
    cut_refusal "cut 1: phi2DPES = nan is not a finite number"
    echo '{"cuts": [{}]}' > "${cut_case}/pes_cuts.json"
    touch "${cut_case}/pes_cut_2.dat"
    cut_refusal "pes_cut_2\.dat is left from more cuts"
    rm "${cut_case}/pes_cut_2.dat"
    printf 'Title : "x";\nx2DPES : REAL : 1.0;\n' > "${cut_case}/DVR.inp"
    cut_refusal "DVR\.inp:2 sets x2DPES itself"
    ;;
  unmodified-targets)
    repo="${work}/repo"
    scratch_repo "${repo}" run_propagation_reference.sh run_fine_analysis_reference.sh
    case_dir="${repo}/tests/reference/case"
    mkdir -p "${case_dir}/fine_full"
    printf '%s\n' 'writeSMatrix : LOGICAL : FALSE;' > "${case_dir}/DVR.inp"
    for file in fortran_stdout.txt fortran_averages.out fortran_fort21.bin fine_full/stdout.txt; do
      echo "committed" > "${case_dir}/${file}"
    done
    commit_all "${repo}" "the references of a case"
    check=(in_config_of "${repo}/tools/fortran_reference/reference_config.sh" require_unmodified)
    # What git can restore passes: unchanged files and directories, absent and deleted files.
    rm "${case_dir}/fortran_averages.out"
    expect_success "${check[@]}" "${case_dir}/fortran_stdout.txt" "${case_dir}/fine_full" "${case_dir}/absent.txt" \
      "${case_dir}/fortran_averages.out"
    echo "edited" >> "${case_dir}/fortran_stdout.txt"
    expect_refusal "case/fortran_stdout\.txt: differs from the last commit" "${check[@]}" "${case_dir}/fortran_stdout.txt"
    echo "edited" >> "${case_dir}/fine_full/stdout.txt"
    expect_refusal "case/fine_full: differs from the last commit" "${check[@]}" "${case_dir}/fine_full"
    echo "new" > "${case_dir}/fortran_fort25.txt"
    expect_refusal "fortran_fort25\.txt: not tracked by git" "${check[@]}" "${case_dir}/fortran_fort25.txt"
    git -C "${repo}" add "${case_dir}/fortran_fort25.txt"
    expect_refusal "fortran_fort25\.txt: differs from the last commit" "${check[@]}" "${case_dir}/fortran_fort25.txt"
    echo "new" > "${case_dir}/fine_full/fort.3"
    git -C "${repo}" checkout -q -- "${case_dir}/fine_full/stdout.txt"
    expect_refusal "fine_full: not tracked by git \(tests/reference/case/fine_full/fort\.3\)" "${check[@]}" "${case_dir}/fine_full"
    expect_refusal "require_unmodified: .*elsewhere\.txt is not inside" "${check[@]}" "${work}/elsewhere.txt"
    GIT_DIR="${work}/no-git" expect_refusal "is not a git checkout" "${check[@]}" "${case_dir}/fortran_stdout.txt"
    # The run scripts check their targets before they run: the edited file stays as it was.
    make_install "${work}/install" "touch averages.out fort.21 fort.25 fort.10"
    export SPODVR_FORTRAN_REF_DIR="${work}/install" SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1
    scripts="${repo}/tools/fortran_reference"
    expect_refusal "fortran_stdout\.txt: differs" "${scripts}/run_propagation_reference.sh" case "${work}/run1"
    [[ ! -e "${work}/run1" && "$(tail -n 1 "${case_dir}/fortran_stdout.txt")" == edited ]] ||
      fail "the refused propagation ran or overwrote the edited fortran_stdout.txt"
    expect_refusal "fine_full: not tracked by git" "${scripts}/run_fine_analysis_reference.sh" case full
    # Deleted (as regenerate_all.sh leaves them) the targets pass, and the runs write them.
    rm "${case_dir}/fortran_stdout.txt" "${case_dir}/fortran_fort25.txt" "${case_dir}/fine_full/fort.3"
    git -C "${repo}" rm -q --cached "${case_dir}/fortran_fort25.txt"
    expect_success "${scripts}/run_propagation_reference.sh" case "${work}/run2"
    [[ -f "${case_dir}/fortran_stdout.txt" && -f "${case_dir}/fortran_fort25.txt" ]] ||
      fail "the propagation run did not write its outputs"
    ;;
  *)
    echo "usage: test_reference_tools.sh CASE WORK [CXX] (unknown case ${test_case})" >&2
    exit 2
    ;;
esac
echo "PASS: ${test_case}"

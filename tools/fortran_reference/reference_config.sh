# shellcheck shell=bash
# Fortran reference configuration: the one home of the pins (the reference install every script
# runs, the upstream commit build_fortran_reference.sh builds), the paths, the thread count, the
# environment of every Fortran run, and the functions the scripts of this directory share.
#
# Sourced, never executed: source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"
# (from a script in this directory; crp_driver/ uses "/../reference_config.sh").
# Environment read:
#   SPODVR_FORTRAN_REF_DIR              the reference install to run, and the one
#                                       build_fortran_reference.sh installs into
#                                       (default: ${HOME}/work/SPO-DVR_ref/<id>)
#   SPODVR_FORTRAN_MIRROR               the upstream git mirror (default ${HOME}/work/SPO-DVR)
#   SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT  1: run although the installed run-time packages differ
#                                       from the install's BUILDINFO.tsv (default: refuse)
# Retired (refused when set): the variables in FORTRAN_RETIRED_VARS below.
# Run time: sourcing, instant; check_install hashes the install, under a second.
# Docs: tools/fortran_reference/README.md.

shopt -s inherit_errexit  # a failure inside $(...) stops the caller too
export LC_ALL=C

# Refusals: "error: <what>; <fix>" on stderr, exit status 1.
refuse() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

# The single override of the reference location is SPODVR_FORTRAN_REF_DIR; the per-file variables
# of the older scripts are refused so that no run silently uses another executable or object set.
FORTRAN_RETIRED_VARS=(SPODVR_FORTRAN_EXE SPODVR_FORTRAN_OBJS)
for retired_var in "${FORTRAN_RETIRED_VARS[@]}"; do
  [[ -z "${!retired_var+set}" ]] ||
    refuse "${retired_var} is retired; unset it and point SPODVR_FORTRAN_REF_DIR at a reference install"
done
unset retired_var

FORTRAN_REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

# The reference every run script uses: an install with the fixed layout bin/spodvr.x, OBJS/ (*.o,
# *.mod), BUILDINFO.tsv and SHA256SUMS (README.md of this directory). v5.19 is the archive build.
FORTRAN_REF_ROOT="${HOME}/work/SPO-DVR_ref"
FORTRAN_REF_ID="v5.19"
FORTRAN_REF_DIR="${SPODVR_FORTRAN_REF_DIR:-${FORTRAN_REF_ROOT}/${FORTRAN_REF_ID}}"
FORTRAN_REF_EXE="${FORTRAN_REF_DIR}/bin/spodvr.x"
FORTRAN_REF_OBJS="${FORTRAN_REF_DIR}/OBJS"
# Threads of every Fortran run, whatever the caller's OMP_NUM_THREADS.
FORTRAN_REF_THREADS=4

# The upstream commit build_fortran_reference.sh builds: svn r1713 (5.24 plus the fine-analysis
# sign fix). The tree hash identifies the content even if the git-svn mirror is regenerated.
FORTRAN_MIRROR="${HOME}/work/SPO-DVR"
FORTRAN_MIRROR="${SPODVR_FORTRAN_MIRROR:-${FORTRAN_MIRROR}}"
FORTRAN_BUILD_COMMIT="dce0ecc44a401098528e9ef94c54e1086d393ebb"
FORTRAN_BUILD_TREE="3784a9530904190af7a804448cec0260f10edcf6"
FORTRAN_BUILD_SVN="r1713"
FORTRAN_BUILD_VERSION="5.24"
FORTRAN_BUILD_ID="${FORTRAN_BUILD_SVN}-${FORTRAN_BUILD_COMMIT:0:7}"
FORTRAN_BUILD_INSTALL_DIR="${SPODVR_FORTRAN_REF_DIR:-${FORTRAN_REF_ROOT}/${FORTRAN_BUILD_ID}}"
# Line 19 of sources.dir/GNU.mk at that commit (with its trailing space), and its replacement:
# -O3 without -march/-mtune, so no FMA, no AVX libmvec and no crtfastmath.
FORTRAN_BUILD_GNU_MK_LINE="O3FLAGS             = -Ofast -mtune=native -march=native "
FORTRAN_BUILD_GNU_MK_REPLACEMENT="O3FLAGS = -O3"
# What the pinned Makefile settings (COMP=GNU, LIB=FFTW, OpenMP, no MPI) give with that line;
# the build log must show exactly these (build_fortran_reference.sh, the allow-lists).
FORTRAN_BUILD_DEFINES='-DX86 -DGNU -DNO_SCSL -DNO_ESSL -DFFTWL -DNO_MKL -DUSE_OMP -DNO_USE_MPI -D__COMP__=\"GNU\" -D__PROGRAM_VERSION__=\"'"${FORTRAN_BUILD_VERSION}"'\"'
FORTRAN_BUILD_FLAGS="-fdefault-real-8 -ffree-line-length-none -fopenmp -O3"
FORTRAN_BUILD_INCLUDES="-I/usr/include"
FORTRAN_BUILD_EXE_NAME="SPO-DVR_${FORTRAN_BUILD_VERSION}_O3_GNU_OMP_FFTW.x"
# Packages recorded in BUILDINFO.tsv: the toolchain, and the libraries the executable loads.
FORTRAN_BUILD_PACKAGES=(cpp-15 gfortran-15 libc6-dev libfftw3-dev)
FORTRAN_RUNTIME_PACKAGES=(libc6 libfftw3-double3 libgfortran5 libgomp1)

# PATH of the fixed environment (run_clean).
FORTRAN_ENV_PATH="/usr/bin:/bin"

# run_clean COMMAND...: COMMAND in a fixed environment (env -i): none of the caller's OMP_*,
# GOMP_*, GFORTRAN_* (GFORTRAN_CONVERT_UNIT would flip the fort.21 byte order), LD_*, CPATH or
# MAKEFLAGS reaches the Fortran build or run.
run_clean() {
  env -i PATH="${FORTRAN_ENV_PATH}" HOME="${HOME}" LC_ALL=C TMPDIR="${TMPDIR:-/tmp}" \
    OMP_NUM_THREADS="${FORTRAN_REF_THREADS}" "$@"
}

# installed_version PACKAGE: the dpkg version of PACKAGE (native architecture), or "not-installed".
installed_version() {
  dpkg-query -W -f='${Version}' "$1:$(dpkg --print-architecture)" 2> /dev/null || echo "not-installed"
}

# buildinfo_has FILE KEY VALUE: FILE has the row "KEY<TAB>VALUE", optionally followed by " (note)".
buildinfo_has() {
  grep -qxF "$2"$'\t'"$3" "$1" || grep -qF "$2"$'\t'"$3 (" "$1"
}

# check_install: the reference install is complete and unchanged (every file against SHA256SUMS,
# the executable also against BUILDINFO.tsv), and the run-time packages are the recorded ones.
check_install() {
  local buildinfo="${FORTRAN_REF_DIR}/BUILDINFO.tsv"
  [[ -f "${buildinfo}" && -f "${FORTRAN_REF_DIR}/SHA256SUMS" && -x "${FORTRAN_REF_EXE}" ]] ||
    refuse "${FORTRAN_REF_DIR} is not a reference install (bin/spodvr.x, BUILDINFO.tsv, SHA256SUMS);" \
      "build one with build_fortran_reference.sh or set SPODVR_FORTRAN_REF_DIR"
  grep -qxF $'key\tvalue' "${buildinfo}" ||
    refuse "${buildinfo} has no header row \"key<TAB>value\"; restore the install's BUILDINFO.tsv"
  # A file edited, replaced or lost after the build (the symlinked v5.19 archive included).
  (cd "${FORTRAN_REF_DIR}" && sha256sum --check --quiet --strict SHA256SUMS) ||
    refuse "${FORTRAN_REF_DIR} differs from its SHA256SUMS; rebuild the install"
  # SHA256SUMS rewritten together with a different executable.
  local exe_sha
  exe_sha="$(sha256sum < "${FORTRAN_REF_EXE}")"
  buildinfo_has "${buildinfo}" exe.sha256 "${exe_sha%% *}" ||
    refuse "${FORTRAN_REF_EXE} is not the executable of ${buildinfo} (exe.sha256); rebuild the install"
  # A libm/libgfortran/libgomp/FFTW upgrade can move results by an ulp.
  local drift
  drift="$(runtime_drift)"
  [[ "${drift}" != - ]] || return 0
  if [[ "${SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT:-}" == 1 ]]; then
    echo "warning: run-time packages ${drift} differ from ${buildinfo} (recorded in MANIFEST.tsv)" >&2
  else
    refuse "run-time packages ${drift} differ from ${buildinfo};" \
      "rebuild the install, or set SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT=1 to run anyway"
  fi
}

# runtime_drift: the run-time packages whose installed version is not the one in the install's
# BUILDINFO.tsv, as "package=installed-version ..."; "-" when there are none.
runtime_drift() {
  local package installed drift=()
  for package in "${FORTRAN_RUNTIME_PACKAGES[@]}"; do
    installed="$(installed_version "${package}")"
    buildinfo_has "${FORTRAN_REF_DIR}/BUILDINFO.tsv" "package.runtime.${package}" "${installed}" ||
      drift+=("${package}=${installed}")
  done
  if ((${#drift[@]} > 0)); then echo "${drift[*]}"; else echo -; fi
}

# record_manifest RUN_DIR CASE_DIR PROGRAM INPUT... -- FILE...: writes the provenance lines of the
# FILEs (relative to CASE_DIR, just generated there) into CASE_DIR/MANIFEST.tsv, replacing their
# old lines (refdata.py manifest record): the calling script, the install, PROGRAM (the program
# that ran), the INPUTs as run (files or data-set directories, relative to RUN_DIR or absolute)
# and the run-time drift.
record_manifest() {
  local run_dir="$1" case_dir="$2" program="$3" inputs=()
  shift 3
  while (($# > 0)) && [[ "$1" != "--" ]]; do
    inputs+=(--input "$1")
    shift
  done
  (($# > 1)) || refuse "record_manifest: no \"--\" and file list after the inputs; fix the call"
  shift
  (cd "${run_dir}" && python3 -I -B "${FORTRAN_REPO_DIR}/tools/fortran_reference/refdata.py" manifest record \
    "${case_dir}" "$@" --script "$(basename "$0")" --install "${FORTRAN_REF_DIR}" --binary "${program}" \
    --runtime-drift "$(runtime_drift)" "${inputs[@]}")
}

# check_build_package PACKAGE: the installed PACKAGE is the version that built the install ("-":
# not recorded, accepted). Driver builds use it: their .mod files are compiler-version specific.
check_build_package() {
  local buildinfo="${FORTRAN_REF_DIR}/BUILDINFO.tsv" installed
  installed="$(installed_version "$1")"
  buildinfo_has "${buildinfo}" "package.build.$1" - && return 0
  buildinfo_has "${buildinfo}" "package.build.$1" "${installed}" ||
    refuse "build package $1 is ${installed}, not the version in ${buildinfo}; install that version"
}

# new_scratch_dir [DIR]: creates DIR (default: a new mktemp -d directory) and prints its absolute
# path. An existing DIR is refused: a stale .pes or pesmorphparameters.dat there would be read.
new_scratch_dir() {
  if [[ -z "${1:-}" ]]; then
    mktemp -d
    return
  fi
  [[ ! -e "$1" ]] || refuse "scratch directory $1 exists; name a new directory or delete it"
  mkdir "$1" || refuse "cannot create scratch directory $1; its parent must exist"
  (cd "$1" && pwd)
}

# stage_inputs CASE_DIR DIR: copies the case's DVR.inp into DIR and, for a CRP library case (key
# CRPlibDataSet of its input.json, relative to the repository), the data set as DIR/dataset: the
# short relative path its DVR.inp names (the Fortran fails on paths above 50 characters).
stage_inputs() {
  local case_dir="$1" dir="$2" data_set=""
  [[ -f "${case_dir}/DVR.inp" ]] || refuse "${case_dir}/DVR.inp not found; name a case of tests/reference"
  cp "${case_dir}/DVR.inp" "${dir}/DVR.inp"
  if [[ -f "${case_dir}/input.json" ]]; then
    data_set="$(python3 -I -B -c 'import json, sys; print(json.load(open(sys.argv[1])).get("CRPlibDataSet", ""))' \
      "${case_dir}/input.json")"
  fi
  [[ -n "${data_set}" ]] || return 0
  [[ -f "${FORTRAN_REPO_DIR}/${data_set}/CRPData" ]] ||
    refuse "CRP data set ${data_set} (CRPlibDataSet of ${case_dir}/input.json) not found; restore it"
  cp -r "${FORTRAN_REPO_DIR}/${data_set}" "${dir}/dataset"
}

# set_field FILE NAME VALUE: sets the value of the line "NAME : TYPE : value;" of a DVR.inp,
# keeping its spacing.
set_field() {
  grep -q "^$2 *:" "$1" || refuse "field $2 not found in $1; add it to the case's DVR.inp"
  sed -i -E "s/^($2 *:[^:]*: *)[^;]*;/\1$3;/" "$1"
  grep -q "^$2 *:[^:]*: *$3;" "$1" || refuse "field $2 of $1 is not \"$2 : TYPE : value;\"; give it a type"
}

# run_fortran DIR STDOUT REQUIRED...: runs the reference executable in DIR with stdout and stderr
# to DIR/STDOUT, after check_install. Refuses a failed run and a missing REQUIRED output (a glob
# relative to DIR): a Fortran STOP on bad input exits with status 0.
run_fortran() {
  local dir="$1" stdout="$2" status=0 pattern
  shift 2
  check_install
  echo "* Running ${FORTRAN_REF_EXE} in ${dir} (${FORTRAN_REF_THREADS} threads)"
  (cd "${dir}" && run_clean "${FORTRAN_REF_EXE}" > "${stdout}" 2>&1) || status=$?
  ((status == 0)) || refuse "the Fortran run in ${dir} failed (exit ${status}); see ${dir}/${stdout}"
  for pattern in "$@"; do
    compgen -G "${dir}/${pattern}" > /dev/null ||
      refuse "the Fortran run wrote no ${pattern} in ${dir}; see ${dir}/${stdout}"
  done
}

# check_instructions FILE...: the executables or objects hold no VEX/EVEX (v-prefixed SIMD)
# instruction, import no AVX libmvec variant (_ZGVc/_ZGVd/_ZGVe) and do not link crtfastmath
# (set_fast_math): the marks of -march/-mtune, FMA and -Ofast/-ffast-math.
check_instructions() {
  local file found
  for file in "$@"; do
    objdump -f "${file}" > /dev/null || refuse "${file} is not an object file or executable; name a built one"
    found="$(objdump -d --no-show-raw-insn "${file}" | grep -m 3 -E '^ +[0-9a-f]+:[[:space:]]+v[a-z0-9]+' || true)"
    [[ -z "${found}" ]] || refuse "${file} has VEX/EVEX instructions, e.g. ${found//$'\n'/ | }; build with -O3, no -march"
    found="$(nm --undefined-only --format=just-symbols "${file}" | grep -E '^_ZGV[cde]' | paste -sd ' ' || true)"
    [[ -z "${found}" ]] || refuse "${file} imports AVX libmvec functions ${found}; build without -march/-ffast-math"
    # grep without -q reads to the end: an early exit would fail nm with SIGPIPE under pipefail.
    ! nm "${file}" 2> /dev/null | grep -w set_fast_math > /dev/null ||
      refuse "${file} links crtfastmath (set_fast_math); build without -Ofast/-ffast-math"
  done
}

# build_fortran_driver DIR NAME SOURCE... -- OBJECT...: builds the driver NAME in the new directory
# DIR from SOURCE (.f90 with the reference build's flags, .c with -O2), linked against the named
# objects of the reference install (OBJS/<OBJECT>.o), and prints the executable's path last.
build_fortran_driver() {
  local dir="$1" name="$2" source objects=() fortran_objects=()
  shift 2
  check_install
  check_build_package gfortran-15
  dir="$(new_scratch_dir "${dir}")"
  while (($# > 0)) && [[ "$1" != "--" ]]; do
    source="$1"
    shift
    case "${source}" in
      *.c) run_clean gcc -O2 -c "${source}" -o "${dir}/$(basename "${source}" .c).o" ;;
      *.f90)
        # shellcheck disable=SC2086  # FORTRAN_BUILD_FLAGS is a list of flags
        run_clean gfortran ${FORTRAN_BUILD_FLAGS} -I"${FORTRAN_REF_OBJS}" -J"${dir}" \
          -c "${source}" -o "${dir}/$(basename "${source}" .f90).o" ;;
      *) refuse "driver source ${source} is neither .f90 nor .c; fix the build_fortran_driver call" ;;
    esac
    objects+=("${dir}/$(basename "${source%.*}").o")
  done
  (($# > 0)) || refuse "build_fortran_driver: no \"--\" before the object list; fix the call"
  shift
  for source in "$@"; do fortran_objects+=("${FORTRAN_REF_OBJS}/${source}.o"); done
  run_clean gfortran -fopenmp -o "${dir}/${name}" "${objects[@]}" "${fortran_objects[@]}"
  check_instructions "${dir}/${name}"
  echo "${dir}/${name}"
}

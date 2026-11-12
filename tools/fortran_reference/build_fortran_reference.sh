#!/usr/bin/env bash
# Builds the pinned upstream Fortran (reference_config.sh: FORTRAN_BUILD_*) into a read-only
# reference install, reproducibly and with provenance:
#   1. checks the pin in the mirror (commit, tree hash, svn revision) and writes `git archive` of
#      sources.dir to a tarball, the only source of the build;
#   2. extracts it into a new mktemp directory next to the install, replaces the exact GNU.mk line
#      of -Ofast -mtune=native -march=native by -O3, runs DoBuild in a fixed environment
#      (env -i, SOURCE_DATE_EPOCH = the commit time);
#   3. refuses a log with an error or with a line outside the allow-lists (PreProcess, compile and
#      link lines with exactly the pinned defines and flags, DoBuild's messages), and an executable
#      or object with VEX/EVEX instructions, AVX libmvec imports or crtfastmath;
#   4. installs bin/spodvr.x, OBJS/ (*.o, *.mod), BUILDINFO.tsv, SHA256SUMS, build.log and
#      sources.tar with one mv, makes them read-only, and copies BUILDINFO.tsv to
#      tools/fortran_reference/buildinfo/<id>.tsv (commit it).
#
# Usage: tools/fortran_reference/build_fortran_reference.sh
#        tools/fortran_reference/build_fortran_reference.sh --check-log LOG
#        tools/fortran_reference/build_fortran_reference.sh --check-binary FILE...
#   (no argument)        build and install into SPODVR_FORTRAN_REF_DIR, by default
#                        ${HOME}/work/SPO-DVR_ref/<FORTRAN_BUILD_ID>; an existing install is refused
#   --check-log LOG      only step 3's log checks, on a DoBuild log of the pinned version
#   --check-binary FILE  only step 3's instruction check
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_MIRROR (reference_config.sh); TMPDIR.
# Run time: about 40 s (DoBuild compiles serially; 38 s measured 2026-11-12); the checks alone, seconds.
# Docs: tools/fortran_reference/README.md.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/reference_config.sh"

usage() {
  echo "usage: build_fortran_reference.sh [--check-log LOG | --check-binary FILE...]" >&2
  exit 2
}

# ERE of a literal string, spaces turned into " +" (the Makefile pads with variable spacing).
ere() {
  printf '%s' "$1" | sed -e 's/[][\.*^$+?(){}|]/\\&/g' -e 's/ / +/g'
}

# check_log LOG: refuses a DoBuild log that shows a failure, or a line outside the allow-lists.
check_log() {
  local log="$1" bad
  [[ -f "${log}" ]] || refuse "build log ${log} not found; name an existing file"
  # A failed compile or make step: DoMake exits with status 0 even then.
  bad="$(grep -E '^Error|: Error|\*\*\*' "${log}" | head -n 5 || true)"
  [[ -z "${bad}" ]] || refuse "the build failed, ${log} shows: ${bad//$'\n'/ | }; fix the build"
  # Any other define, flag, compiler or link input (-Ofast, -march, -DMKL, -D__EXHAUSTIVE_FFTW_PLANS__,
  # a DEBUG or MPI build, a warning): every line must be one of these, matched whole (grep -x).
  local flags source
  flags="$(ere "${FORTRAN_BUILD_FLAGS} ${FORTRAN_BUILD_INCLUDES}")"
  source='[A-Za-z0-9_-]+\.f90'
  local preprocess="\./PreProcess (${source}) \1 $(ere "${FORTRAN_BUILD_DEFINES}")"
  local compile="gfortran +${flags} +-c +${source}"
  local link="gfortran +${flags} +-o +$(ere "${FORTRAN_BUILD_EXE_NAME}")( +[A-Za-z0-9_-]+\.o)+ +-lfftw3 *"
  # DoBuild/DoMake progress messages; "[[: not found" is upstream's check_variables, a bash test
  # that make runs with /bin/sh (dash), where it is skipped.
  local messages='Cleaning up all stuff for doing a complete re-make \.\.\.|Using /usr/bin/gmake for making\.\.\.'
  messages+='|DoMake script running and checking variables\.\.\.|/bin/sh: 1: \[\[: not found'
  messages+='|Creating TEMP directory for the make process\.\.\.|Copying (all source files|all listings'
  messages+='|all object files|make files) to TEMP dir \.\.\.|Starting /usr/bin/gmake in TEMP dir \.\.\.'
  messages+='|Copying back all make generated files \.\.\.|Cleaning up TEMP dir \.\.\.'
  bad="$(grep -vxE -e "${preprocess}" -e "${compile}" -e "${link}" -e "${messages}" "${log}" || true)"
  [[ -z "${bad}" ]] || refuse "$(wc -l <<< "${bad}") line(s) of ${log} are outside the allow-lists, e.g." \
    "$(head -n 3 <<< "${bad}" | cut -c 1-160 | paste -sd '|'); build with the pinned flags only (reference_config.sh)"
  grep -qxE "${link}" "${log}" || refuse "${log} has no link line; the build did not link"
}

# verify_pin: the mirror has the pinned commit, with the pinned tree and svn revision.
verify_pin() {
  local tree
  git -C "${FORTRAN_MIRROR}" rev-parse --verify --quiet "${FORTRAN_BUILD_COMMIT}^{commit}" > /dev/null ||
    refuse "commit ${FORTRAN_BUILD_COMMIT} is not in ${FORTRAN_MIRROR}; git -C ${FORTRAN_MIRROR} pull, or set SPODVR_FORTRAN_MIRROR"
  tree="$(git -C "${FORTRAN_MIRROR}" rev-parse "${FORTRAN_BUILD_COMMIT}^{tree}")"
  [[ "${tree}" == "${FORTRAN_BUILD_TREE}" ]] ||
    refuse "commit ${FORTRAN_BUILD_COMMIT} has tree ${tree}, not ${FORTRAN_BUILD_TREE}; check the mirror"
  git -C "${FORTRAN_MIRROR}" log -1 --format=%B "${FORTRAN_BUILD_COMMIT}" | grep -qF "trunk@${FORTRAN_BUILD_SVN#r} " ||
    refuse "commit ${FORTRAN_BUILD_COMMIT} is not svn trunk ${FORTRAN_BUILD_SVN}; check the pin in reference_config.sh"
}

# set_o3flags DIR: replaces GNU.mk's O3FLAGS line, refusing any other content of that line.
set_o3flags() {
  local makefile="$1/GNU.mk"
  [[ "$(sed -n 19p "${makefile}")" == "${FORTRAN_BUILD_GNU_MK_LINE}" ]] ||
    refuse "line 19 of ${makefile} is not \"${FORTRAN_BUILD_GNU_MK_LINE}\"; update the pin in reference_config.sh"
  sed -i "19s/.*/${FORTRAN_BUILD_GNU_MK_REPLACEMENT}/" "${makefile}"
}

# write_buildinfo INSTALL TARBALL_SHA: BUILDINFO.tsv of the install (format: README.md here).
write_buildinfo() {
  local install="$1" package exe_sha zgvb
  exe_sha="$(sha256sum < "${install}/bin/spodvr.x")"
  zgvb="$(nm --undefined-only --format=just-symbols "${install}/bin/spodvr.x" | grep '^_ZGVb' | sort | paste -sd ' ' || true)"
  {
    echo "# Fortran reference install ${FORTRAN_BUILD_ID}, written by tools/fortran_reference/build_fortran_reference.sh."
    echo "# Keys: tools/fortran_reference/README.md; \"-\" = not recorded."
    printf 'key\tvalue\n'
    {
      printf 'commit\t%s\n' "${FORTRAN_BUILD_COMMIT}"
      printf 'date\t%s\n' "$(date -u +%F)"
      printf 'defines\t%s\n' "${FORTRAN_BUILD_DEFINES}"
      printf 'exe.sha256\t%s\n' "${exe_sha%% *}"
      printf 'flags\t%s\n' "${FORTRAN_BUILD_FLAGS} ${FORTRAN_BUILD_INCLUDES}"
      printf 'id\t%s\n' "${FORTRAN_BUILD_ID}"
      for package in "${FORTRAN_BUILD_PACKAGES[@]}"; do
        printf 'package.build.%s\t%s\n' "${package}" "$(installed_version "${package}")"
      done
      for package in "${FORTRAN_RUNTIME_PACKAGES[@]}"; do
        printf 'package.runtime.%s\t%s\n' "${package}" "$(installed_version "${package}")"
      done
      printf 'source\t%s\n' "upstream SPO-DVR git (git-svn mirror of the Leiden SVN trunk): git archive <commit> sources.dir"
      printf 'svn\t%s\n' "${FORTRAN_BUILD_SVN}"
      printf 'tarball.sha256\t%s\n' "$2"
      printf 'tree\t%s\n' "${FORTRAN_BUILD_TREE}"
      printf 'version\t%s\n' "${FORTRAN_BUILD_VERSION}"
      printf 'zgvb\t%s\n' "${zgvb:--}"
    } | sort
  } > "${install}/BUILDINFO.tsv"
}

# build: steps 1 to 4 of the header.
build() {
  local target="${FORTRAN_BUILD_INSTALL_DIR}" work tarball_sha epoch
  [[ ! -e "${target}" ]] ||
    refuse "${target} exists; to rebuild: chmod -R u+w ${target} && rm -r ${target}"
  verify_pin
  mkdir -p "$(dirname "${target}")"
  # Next to the install, on the same file system: the final mv is atomic.
  work="$(mktemp -d "$(dirname "${target}")/.build.XXXXXX")"
  echo "* Building ${FORTRAN_BUILD_ID} in ${work}"
  git -C "${FORTRAN_MIRROR}" archive --format=tar "${FORTRAN_BUILD_COMMIT}" sources.dir > "${work}/sources.tar"
  tarball_sha="$(sha256sum < "${work}/sources.tar")"
  mkdir "${work}/src"
  tar -xf "${work}/sources.tar" -C "${work}/src"
  set_o3flags "${work}/src/sources.dir"
  epoch="$(git -C "${FORTRAN_MIRROR}" show -s --format=%ct "${FORTRAN_BUILD_COMMIT}")"
  # DoBuild = DoClean + DoMake; cpp honours SOURCE_DATE_EPOCH for __DATE__/__TIME__.
  (cd "${work}/src/sources.dir" && run_clean SOURCE_DATE_EPOCH="${epoch}" ./DoBuild) > "${work}/build.log" 2>&1 ||
    refuse "DoBuild failed; see ${work}/build.log (delete ${work} afterwards)"
  check_log "${work}/build.log"
  # The directory is new, so an executable there was linked by this build.
  local exe="${work}/src/executable.dir/${FORTRAN_BUILD_EXE_NAME}"
  [[ -x "${exe}" ]] || refuse "DoBuild wrote no ${exe}; see ${work}/build.log (delete ${work} afterwards)"
  check_instructions "${exe}" "${work}/src/sources.dir/OBJS.dir/"*.o

  local install="${work}/install"
  mkdir -p "${install}/bin" "${install}/OBJS"
  cp "${exe}" "${install}/bin/spodvr.x"
  cp "${work}/src/sources.dir/OBJS.dir/"*.o "${work}/src/sources.dir/OBJS.dir/"*.mod "${install}/OBJS/"
  mv "${work}/build.log" "${work}/sources.tar" "${install}/"
  write_buildinfo "${install}" "${tarball_sha%% *}"
  (cd "${install}" && sha256sum bin/spodvr.x OBJS/* BUILDINFO.tsv build.log sources.tar > SHA256SUMS)
  mv "${install}" "${target}"
  rm -rf "${work}"
  chmod -R a-w "${target}"
  cp "${target}/BUILDINFO.tsv" "${FORTRAN_REPO_DIR}/tools/fortran_reference/buildinfo/${FORTRAN_BUILD_ID}.tsv"
  chmod u+w "${FORTRAN_REPO_DIR}/tools/fortran_reference/buildinfo/${FORTRAN_BUILD_ID}.tsv"
  echo "* Installed ${target}"
  echo "* Commit tools/fortran_reference/buildinfo/${FORTRAN_BUILD_ID}.tsv"
}

case "${1:-}" in
  "") (($# == 0)) || usage; build ;;
  --check-log) (($# == 2)) || usage; check_log "$2"; echo "* ${2}: passes the log checks" ;;
  --check-binary) (($# >= 2)) || usage; shift; check_instructions "$@"; echo "* passes the instruction check: $*" ;;
  *) usage ;;
esac

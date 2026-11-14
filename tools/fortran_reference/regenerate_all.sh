#!/usr/bin/env bash
# Regenerates every Fortran reference file of the repository from the committed inputs (each
# case's DVR.inp and input.json, the CRP data sets and point lists), in dependency order:
#   1. deletes the files the MANIFEST.tsv of each case lists, and the MANIFESTs;
#   2. setup runs (gates 1, 2), 3. propagation (gates 3, 4), 4. the Fortran restart,
#   5. fine analysis (gate 5, on the fortran_fort21.bin just regenerated),
#   6. the CRP driver on the committed point lists;
# each through its run_*.sh script, which writes the MANIFEST lines. Each step then checks its
# required outputs, and the end runs the ReferenceManifest check. The steps below ARE the list of
# reference cases: a new case adds its step here in the commit that adds the case.
#
# Refuses unless the reference trees are clean (git status: nothing modified or untracked, so the
# result can be compared with git), and checks the install before it deletes anything.
# Scratch directories: one new mktemp -d directory (printed), kept for inspection.
#
# Usage: tools/fortran_reference/regenerate_all.sh   (no arguments)
# Environment: SPODVR_FORTRAN_REF_DIR, SPODVR_FORTRAN_ALLOW_RUNTIME_DRIFT (reference_config.sh);
# TMPDIR (where the scratch directory is made).
# Run time: about 45 min on 4 threads (2381–2538 s with the r1713 install, 2026-11-12, plus 158 s
# for the case cu111_dailight_short_lowz added later).
# Docs: tools/fortran_reference/README.md.
set -euo pipefail
tools_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${tools_dir}/reference_config.sh"

(($# == 0)) || { echo "usage: regenerate_all.sh (no arguments)" >&2; exit 2; }
reference_trees=(tests/reference tests/crplib/reference)

# Untracked or modified files would be deleted or mixed with the regenerated ones.
dirty="$(git -C "${FORTRAN_REPO_DIR}" status --porcelain -- "${reference_trees[@]}")"
[[ -z "${dirty}" ]] ||
  refuse "the reference trees are not clean (git status: ${dirty//$'\n'/, }); commit, restore or remove those files"
# What the runs and the CRP driver build will check, before anything is deleted.
check_install
check_build_package gfortran-15

# require_outputs DIR FILE...: each FILE (relative to DIR) exists again after its deletion. Only
# this catches a run script that silently stops writing a file the gates read: the script checks
# just its own list, and the gate would skip.
require_outputs() {
  local dir="$1" file
  shift
  for file in "$@"; do
    [[ -f "${dir}/${file}" ]] || refuse "the step wrote no ${dir}/${file}; see its output above"
  done
}

setup() {
  "${tools_dir}/run_setup_reference.sh" "$1"
  require_outputs "${FORTRAN_REPO_DIR}/tests/reference/$1" fortran_setup_stdout.txt debug_asyp.dat pes_reference.txt
}

propagation() {
  "${tools_dir}/run_propagation_reference.sh" "$1"
  require_outputs "${FORTRAN_REPO_DIR}/tests/reference/$1" \
    fortran_stdout.txt fortran_averages.out fortran_fort21.bin fortran_fort25.txt
}

restart() {
  "${tools_dir}/run_restart_reference.sh"
  require_outputs "${FORTRAN_REPO_DIR}/tests/reference/hardwall_small" fortran_restart_stdout.txt
}

# fine CASE VARIANT...: variant "empty" has no energy in the confidence interval, so no fort.3 and
# no fort.25.
fine() {
  local case_name="$1" variant
  "${tools_dir}/run_fine_analysis_reference.sh" "$@"
  shift
  for variant in "$@"; do
    require_outputs "${FORTRAN_REPO_DIR}/tests/reference/${case_name}/fine_${variant}" DVR.inp stdout.txt fort.10
    [[ "${variant}" == empty ]] ||
      require_outputs "${FORTRAN_REPO_DIR}/tests/reference/${case_name}/fine_${variant}" fort.3 fort.25
  done
}

# driver CASE: one <name>.out per point list <name>.in (kernels: also tables.out).
driver() {
  local dir="${FORTRAN_REPO_DIR}/tests/crplib/reference/$1" input
  "${tools_dir}/run_crp_reference.sh" "$1"
  [[ "$1" != kernels ]] || require_outputs "${dir}" tables.out
  for input in "${dir}"/*.in; do require_outputs "${dir}" "$(basename "${input}" .in).out"; done
}

export TMPDIR
TMPDIR="$(mktemp -d)"
echo "* Scratch directories in ${TMPDIR}"
start=$(date +%s)

# 1. The generated files of the previous generation, as their MANIFESTs list them.
shopt -s nullglob
manifests=("${FORTRAN_REPO_DIR}"/tests/reference/*/MANIFEST.tsv "${FORTRAN_REPO_DIR}"/tests/crplib/reference/*/MANIFEST.tsv)
shopt -u nullglob
((${#manifests[@]} > 0)) || refuse "no MANIFEST.tsv in ${reference_trees[*]}; restore them from git"
python3 -I -B "${tools_dir}/refdata.py" manifest delete "${manifests[@]%/MANIFEST.tsv}"

# 2. Setup (potential-only) runs.
setup hardwall_small
setup cu111_dailight_short_sym
setup cu111_dailight_small
setup cu111_persson_small
setup crp_c3v_synthetic_small
setup crp_c4v_synthetic_small
setup crp_c6v_synthetic_small
setup crp_cu111_srp48_small

# 3. Propagation runs.
propagation hardwall_short
propagation hardwall_small
propagation cu111_dailight_short
propagation cu111_dailight_short_sym
propagation cu111_dailight_short_lowz
propagation cu111_dailight_small

# 4. The Fortran restart (hardwall_small in two parts).
restart

# 5. Fine analysis of the regenerated BK coefficients.
fine hardwall_small full truncated nen95 small empty
fine cu111_dailight_small full

# 6. The CRP library driver on the committed point lists and data sets.
driver kernels
driver c3v_impulse_15
driver c3v_impulse_2
driver c3v_impulse_21
driver c3v_impulse_27
driver c3v_impulse_6
driver c3v_impulse_9
driver c3v_missing_extra_cut
driver c3v_pairwise
driver c4v_impulse_13
driver c4v_impulse_2
driver c4v_impulse_6
driver c4v_impulse_9
driver c4v_negative_spacing
driver c4v_pairwise
driver c6v_impulse_1
driver c6v_impulse_10
driver c6v_impulse_15
driver c6v_impulse_5
driver c6v_pairwise
driver c6v_pairwise_angstrom
driver cu111_srp48

# Every MANIFEST line names an existing file with the recorded sha256 (ctest ReferenceManifest.FilesMatch).
cmake -DROOT="${FORTRAN_REPO_DIR}" -P "${tools_dir}/ReferenceManifest.cmake"
echo "* All references regenerated in $(($(date +%s) - start)) s; compare them with git (git status, git diff)"

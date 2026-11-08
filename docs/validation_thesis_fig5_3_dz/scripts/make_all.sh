#!/usr/bin/env bash
# Regenerate the validation report docs/validation_thesis_fig5_3_dz/.
#   make_all.sh            from the extracts in data/: results.json, figures, tables, checks
#   make_all.sh --extract  first re-read the run directories and logs (runs/, ~ minutes)
#   make_all.sh --strict   first the strict-build GPU run: its replays, cleaned log, strict_gpu.py
#   make_all.sh --angular  first the angular-matrix test: helpers, Fortran driver, angular_matrices.py
#                          (needs the extract and strict logs and the potential of the heavy pass; ~1 min)
#   make_all.sh --heavy    also rebuild the C++ helpers and redo the potential, replay,
#                          onlyAnalysis, final-psi, thesis-extraction and sha256 passes (~15 min;
#                          needs the build directories, runs/, ~/work/SPO-DVR, uv for pymupdf)
# Python runs with -B (no __pycache__ in the docs folder). Run from anywhere.
set -euo pipefail
cd "$(dirname "$0")"
PY="python3 -B"
REPO="$(cd ../../.. && pwd)"
WORK="$REPO/runs/thesis_fig5_3_dz_validation"
mode="${1:-default}"

if [[ "$mode" == "--heavy" ]]; then
  for b in build-fast build-cuda-fast; do
    $PY build_helper.py "$b" potential_check.cpp
    $PY build_helper.py "$b" replay.cpp
    OMP_NUM_THREADS=8 "$WORK/helpers/$b/potential_check" "$REPO/inputs/thesis_fig5_3_dz/input.json" "$WORK/potential/$b"
  done
  for f in v_stored.bin v_full.bin grids.txt symmetry_maps.txt; do
    cmp "$WORK/potential/build-fast/$f" "$WORK/potential/build-cuda-fast/$f"   # the two builds must agree bitwise
  done
  for f in asymptotic_potential.dat gas_phase_states.dat; do
    cmp "$WORK/potential/build-fast/$f" "$REPO/runs/thesis_fig5_3_dz_cpu_fast/$f"   # helper self-check
  done
  $PY item04_potential.py build-fast
  for b in build-fast build-cuda-fast; do
    OMP_NUM_THREADS=8 "$WORK/helpers/$b/replay" "$REPO/inputs/thesis_fig5_3_dz/input.json" "$REPO/runs/thesis_fig5_3_dz_cpu_fast/bk_coefficients.bin" native "$WORK/replay/$b/cpu"
    OMP_NUM_THREADS=8 "$WORK/helpers/$b/replay" "$REPO/inputs/thesis_fig5_3_dz/input.json" "$REPO/runs/thesis_fig5_3_dz_cuda_fast/bk_coefficients.bin" native "$WORK/replay/$b/gpu"
    OMP_NUM_THREADS=8 "$WORK/helpers/$b/replay" "$REPO/inputs/thesis_fig5_3_dz/input.json" "$REPO/runs/thesis_fig5_3_dz_fortran/fort.21" fortran "$WORK/replay/$b/fortran"
  done
  mkdir -p "$WORK/oa"
  for spec in fortran:thesis_fig5_3_dz_fortran:fortran cpu:thesis_fig5_3_dz_cpu_fast:native gpu:thesis_fig5_3_dz_cuda_fast:native; do
    IFS=: read -r name dir fmt <<< "$spec"
    python3 -I -c 'import json,sys;j=json.load(open(sys.argv[1]));j.update(onlyAnalysis=True,RestartDir=sys.argv[2],RestartFormat=sys.argv[3]);json.dump(j,open(sys.argv[4],"w"),indent=1)' \
      "$REPO/inputs/thesis_fig5_3_dz/input.json" "$REPO/runs/$dir" "$fmt" "$WORK/oa/oa_$name.json"
    (cd "$REPO" && OMP_NUM_THREADS=8 build-fast/spodvr "$WORK/oa/oa_$name.json" --threads 8 --no-color \
      --output-dir "runs/thesis_fig5_3_dz_oa_$name" > "$WORK/oa/oa_$name.log")
  done
  $PY negative_control.py
  uv run --no-project --with pymupdf python -I thesis_extract.py "$HOME/work/mf_somers_phd_thesis.pdf"
  $PY item17_thesis.py --digitise
  $PY item13_propagation.py --psi
  $PY manifest.py
fi

if [[ "$mode" == "--extract" || "$mode" == "--heavy" ]]; then
  mkdir -p "$WORK/logs"
  for k in cpu:cpu_fast gpu:cuda_fast fortran:fortran; do
    IFS=: read -r short run <<< "$k"
    $PY clean_log.py "$REPO/runs/thesis_fig5_3_dz_$run.log" "$WORK/logs/$short.txt"
    gzip -9 -n -c "$REPO/runs/thesis_fig5_3_dz_$run.log" > "../data/logs/thesis_fig5_3_dz_$run.log.gz"
  done
  cp "$WORK/logs/cpu.txt" ../data/logs/cleaned_cpu_fast.txt
  cp "$WORK/logs/gpu.txt" ../data/logs/cleaned_cuda_fast.txt
  gzip -9 -n -c "$WORK/logs/fortran.txt" > ../data/logs/cleaned_fortran.txt.gz
  $PY provenance.py
  $PY gates.py
  $PY item11_psi_k.py
  $PY setup_items.py
  $PY item13_propagation.py
  $PY item13_propagation.py --derive-psi
  $PY item13_norm_bias.py
  $PY item18b_replay.py
  $PY item14_16_analysis.py
  $PY item18a_only_analysis.py
  $PY item19_physics.py
fi

if [[ "$mode" == "--strict" || "$mode" == "--heavy" ]]; then
  # strict-build GPU run (criteria.md, 2026-11-08): replay its BK file with the build-fast helper
  # (as the other runs; via a header copy) and with a helper linked against build-cuda (self-check)
  $PY build_helper.py build-cuda replay.cpp
  $PY strict_bk_copy.py      # header copy for the build-fast helper (1-ulp ZOptFactor, see the script)
  OMP_NUM_THREADS=8 "$WORK/helpers/build-fast/replay" "$REPO/inputs/thesis_fig5_3_dz/input.json" "$WORK/strict_bk/bk_coefficients.bin" native "$WORK/replay/build-fast/gpu_strict"
  OMP_NUM_THREADS=8 "$WORK/helpers/build-cuda/replay" "$REPO/inputs/thesis_fig5_3_dz/input.json" "$REPO/runs/thesis_fig5_3_dz_cuda/bk_coefficients.bin" native "$WORK/replay/build-cuda/gpu_strict"
  mkdir -p "$WORK/logs"
  $PY clean_log.py "$REPO/runs/thesis_fig5_3_dz_cuda.log" "$WORK/logs/gpu_strict.txt"
  gzip -9 -n -c "$REPO/runs/thesis_fig5_3_dz_cuda.log" > ../data/logs/thesis_fig5_3_dz_cuda.log.gz
  cp "$WORK/logs/gpu_strict.txt" ../data/logs/cleaned_cuda_strict.txt
  $PY strict_gpu.py
fi

if [[ "$mode" == "--angular" || "$mode" == "--heavy" ]]; then
  # angular matrices (criteria.md, 2026-11-08): each code's FBR<->DVR round trip on the thesis grid
  for b in build-cuda build-fast build-cuda-fast; do
    $PY build_helper.py "$b" angular_check.cpp
    OMP_NUM_THREADS=8 "$WORK/helpers/$b/angular_check" "$REPO/inputs/thesis_fig5_3_dz/input.json" "$WORK/angular/$b"
  done
  O="$HOME/work/SPO-DVR/procfiles.dir/OBJS.dir"     # the reference executable's objects (read only)
  F="$WORK/angular/fortran"
  mkdir -p "$F"
  gfortran -O3 -fopenmp -fdefault-real-8 -ffree-line-length-none -I"$O" -J"$F" -c angular_fortran.f90 -o "$F/angular_fortran.o"
  gfortran -O3 -fopenmp "$F/angular_fortran.o" "$O"/{transform,angulargrid,legendre,ffts,fftw_2_scs,constants,errortrap,MPI_wrapper,x86_hack}.o \
    -lfftw3 -o "$F/angular_fortran"
  (cd "$F" && OMP_NUM_THREADS=8 ./angular_fortran)
  $PY angular_matrices.py
fi

$PY item17_thesis.py
$PY results.py
$PY fig1_reaction.py
$PY fig4_time.py
$PY fig3_bk_growth.py
$PY fig2_ladder.py
$PY results.py            # fig4 writes the ladder extract
$PY tables.py
$PY fill_tables.py
$PY check_markdown.py
$PY check_numbers.py

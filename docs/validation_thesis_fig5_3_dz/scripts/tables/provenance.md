| Run | Executable (built) | Build | Started | Wall | Output |
|---|---|---|---|---|---|
| C++ CPU | `build-fast/spodvr` (2026-11-06 15:56) | FAST_MATH ON, ARCH native, oneMKL | 2026-11-06 19:29 | 8:14:37 | `runs/thesis_fig5_3_dz_cpu_fast/` |
| C++ GPU | `build-cuda-fast/spodvr` (2026-11-06 15:58) | FAST_MATH ON, ARCH native, CUDA, cuFFT | 2026-11-06 18:12 | 1:16:47 | `runs/thesis_fig5_3_dz_cuda_fast/` |
| Fortran 5.19 | `~/work/SPO-DVR/EXECUTABLE.DIR/SPO-DVR_5.19_O3_GFORTRAN_OMP_FFTW.x` (2026-11-06 11:39) | -O3 -fopenmp -fdefault-real-8 (GFORTRAN.mk), FFTW_ESTIMATE plans | 2026-11-07 03:43 | 15:47:58 | `runs/thesis_fig5_3_dz_fortran/` |

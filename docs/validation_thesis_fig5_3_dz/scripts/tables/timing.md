| Run | Wall | Propagation | s/step | User CPU | CPU busy | Max RSS | Compiler, flags | FFT library | Threads |
|---|---|---|---|---|---|---|---|---|---|
| Fortran | 15:47:58 (56878 s) | 56816 s | 8.12 | 435768 s | 766% | 9.9 GiB | gfortran 15.2 -O3 -fopenmp -fdefault-real-8 | FFTW 3.3.10, ESTIMATE plans | 8 |
| C++ CPU | 8:14:37 (29677 s) | 29666 s | 4.24 | 230485 s | 776% | 7.9 GiB | g++ 15.2, fast-math, -march=native, LTO | oneMKL 2026.2 (FFTW3 interface) | 8 |
| C++ GPU | 1:16:47 (4607 s) | 4593 s | 0.66 | 4672 s | 101% | 6.3 GiB | g++ 15.2 + nvcc 13.3, fast-math (FMA kernels), -march=native | cuFFT 12.3 (CUDA 13.3) | 8 |
| C++ GPU strict | 1:25:47 (5147 s) | 5124 s | 0.73 | 5217 s | 101% | 6.3 GiB | g++ 15.2 + nvcc 13.3, no fast math, no FMA kernels, x86-64-v3 | cuFFT 12.3 (CUDA 13.3) | 8 |

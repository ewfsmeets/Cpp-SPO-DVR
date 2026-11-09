cmake_minimum_required(VERSION 3.24)
project(SpoDvr VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # for clang-tidy
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# =============================================================================================
# Build options (-D<name>=<value>; README.md lists the usual combinations)
# =============================================================================================
option(SPODVR_USE_MKL "Use Intel oneMKL (BLAS/LAPACK, FFTs) when it is found" ON)
set(SPODVR_FFT_LIBRARY "AUTO" CACHE STRING "FFT library: AUTO (oneMKL if found, else FFTW), MKL or FFTW")
set(SPODVR_ARCH "x86-64-v3" CACHE STRING "-march value: x86-64-v3 (AVX2), x86-64-v4 (AVX-512), native, or empty")
option(SPODVR_FAST_MATH "Experimental: -ffast-math (and FMA in the CUDA kernels); changes rounding" OFF)
set(SPODVR_PGO "OFF" CACHE STRING "Profile-guided optimisation: OFF, GENERATE or USE (docs/performance.md)")
set(SPODVR_PGO_DIR "${CMAKE_BINARY_DIR}/pgo-profile" CACHE PATH "Directory of the PGO profile data")
option(SPODVR_USE_CUDA "Build the GPU backend (--device cuda; docs/cuda.md)" OFF)
option(SPODVR_CUDA_SYNC_DEBUG "Debug: wait for and check every CUDA kernel right after its launch" OFF)
option(SPODVR_POISON_UNINITIALISED "Debug: fill new AlignedVector elements with NaN" OFF)
option(SPODVR_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (CI)" OFF)
option(SPODVR_BUILD_TESTS "Build the unit and regression tests" ON)
option(SPODVR_BUILD_TOOLS "Build the developer tools (operator timing, wave function dump/compare)" ON)

# =============================================================================================
# Libraries
# =============================================================================================
find_package(OpenMP REQUIRED COMPONENTS CXX)

# Intel oneMKL, sequential (it is called from inside our own OpenMP threads) and LP64 (32-bit
# integers). When found it provides BLAS/LAPACK and, by default, the FFTs through its FFTW3
# interface. Without it: the system LAPACK and FFTW. See docs/performance.md.
if(SPODVR_USE_MKL)
  set(MKL_THREADING sequential)
  set(MKL_INTERFACE lp64)
  set(MKL_LINK dynamic)
  find_package(MKL CONFIG QUIET HINTS $ENV{MKLROOT} /opt/intel/oneapi/mkl/latest)
endif()
if(NOT MKL_FOUND)
  find_package(LAPACK REQUIRED)
endif()

set(SPODVR_FFT_LIBRARY_USED ${SPODVR_FFT_LIBRARY})
if(SPODVR_FFT_LIBRARY STREQUAL "AUTO")
  if(MKL_FOUND)
    set(SPODVR_FFT_LIBRARY_USED MKL)
  else()
    set(SPODVR_FFT_LIBRARY_USED FFTW)
  endif()
endif()
if(SPODVR_FFT_LIBRARY_USED STREQUAL "MKL" AND NOT MKL_FOUND)
  message(FATAL_ERROR "SPODVR_FFT_LIBRARY=MKL but oneMKL was not found")
endif()
if(SPODVR_FFT_LIBRARY_USED STREQUAL "FFTW")
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
endif()

include(FetchContent)
FetchContent_Declare(nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(nlohmann_json)

# =============================================================================================
# Compiler flags of our own C++ code. They are applied to C++ files only ($<COMPILE_LANGUAGE:CXX>),
# because the GPU backend also contains .cu files, which nvcc compiles with its own flags.
# =============================================================================================
set(warning_flags -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion)
if(SPODVR_WARNINGS_AS_ERRORS)
  list(APPEND warning_flags -Werror)
endif()

# Floating-point semantics. By default no FMA contraction and no fast-math, so that the results
# stay comparable with the Fortran reference and do not depend on the target CPU. The other flags
# change no result, they only allow optimisations (docs/performance.md):
#   -fcx-fortran-rules  complex * and / as gfortran does (no C99 NaN-recovery call); lets GCC
#                       vectorise complex multiplications
#   -fno-math-errno     math functions need not set errno (sqrt is inlined)
#   -fno-trapping-math  no floating-point traps (more vectorisation)
if(SPODVR_FAST_MATH)
  set(fp_flags -ffast-math -fcx-fortran-rules)
else()
  set(fp_flags -ffp-contract=off -fno-fast-math -fcx-fortran-rules -fno-math-errno -fno-trapping-math)
endif()

# Target instruction set. With -ffp-contract=off our results do not depend on it.
if(SPODVR_ARCH)
  list(APPEND fp_flags -march=${SPODVR_ARCH})
endif()

# Profile-guided optimisation: GENERATE builds a program that writes profile data to
# SPODVR_PGO_DIR when it runs; USE rebuilds with that data.
set(pgo_link_flags)
if(SPODVR_PGO STREQUAL "GENERATE")
  list(APPEND fp_flags -fprofile-generate=${SPODVR_PGO_DIR} -fprofile-update=atomic)
  set(pgo_link_flags -fprofile-generate=${SPODVR_PGO_DIR})
elseif(SPODVR_PGO STREQUAL "USE")
  list(APPEND fp_flags -fprofile-use=${SPODVR_PGO_DIR} -fprofile-correction -Wno-missing-profile)
endif()

# Every target of the project links these two (empty) libraries to get the flags above.
add_library(spodvr_warnings INTERFACE)
target_compile_options(spodvr_warnings INTERFACE "$<$<COMPILE_LANGUAGE:CXX>:${warning_flags}>")
add_library(spodvr_fp INTERFACE)
target_compile_options(spodvr_fp INTERFACE "$<$<COMPILE_LANGUAGE:CXX>:${fp_flags}>")
target_link_options(spodvr_fp INTERFACE ${pgo_link_flags})

# Link-time optimisation in Release builds.
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

# =============================================================================================
# spodvr_crplib: the port of the Fortran CRP library (src/crplib/, docs/crp_library.md), a
# self-contained library (standard library only). It is always built with strict floating point,
# also in the fast-math build: it is set-up code, and its tests compare with the Fortran bit for bit.
# =============================================================================================
set(crplib_fp_flags -ffp-contract=off -fno-fast-math -fno-math-errno -fno-trapping-math)
if(SPODVR_ARCH)
  list(APPEND crplib_fp_flags -march=${SPODVR_ARCH})
endif()
add_library(spodvr_crplib_fp INTERFACE)
target_compile_options(spodvr_crplib_fp INTERFACE "$<$<COMPILE_LANGUAGE:CXX>:${crplib_fp_flags}>")

add_library(spodvr_crplib
  src/crplib/dataset.cpp
  src/crplib/fourier_tables.cpp
  src/crplib/list_directed_reader.cpp
  src/crplib/periodic_spline_3d.cpp
  src/crplib/spline.cpp
  src/crplib/units.cpp
)
target_include_directories(spodvr_crplib PUBLIC src)
target_link_libraries(spodvr_crplib PRIVATE spodvr_warnings spodvr_crplib_fp)

# =============================================================================================
# spodvr_core: everything except the program driver and the GPU backend
# =============================================================================================
add_library(spodvr_core
  src/core/loop_schedule.cpp
  src/core/machine_epsilon.cpp
  src/core/numerical_libraries.cpp
  src/core/system_memory.cpp
  src/input/json_input.cpp
  src/input/validate.cpp
  src/input/process.cpp
  src/input/wave_packet_params.cpp
  src/fft/fft_sizes.cpp
  src/fft/fftw_fft.cpp
  src/linalg/lapack_eigensolver.cpp
  src/grids/momentum_grid.cpp
  src/grids/spatial_grids.cpp
  src/grids/associated_legendre.cpp
  src/grids/gauss_legendre.cpp
  src/grids/angular_grid.cpp
  src/grids/rotational_basis.cpp
  src/propagation/angular_transform.cpp
  src/propagation/propagation_tables.cpp
  src/propagation/initial_wave_packet.cpp
  src/propagation/analysis_surface.cpp
  src/propagation/kinetic_operators.cpp
  src/propagation/potential_operator.cpp
  src/propagation/specular_operator.cpp
  src/propagation/split_operator.cpp
  src/propagation/propagation_backend.cpp
  src/propagation/cpu_backend.cpp
  src/propagation/averages.cpp
  src/propagation/propagation.cpp
  src/propagation/memory_estimate.cpp
  src/analysis/analysis_basis.cpp
  src/analysis/balint_kurti.cpp
  src/analysis/fine_analysis.cpp
  src/analysis/s_matrix.cpp
  src/io/checkpoint.cpp
  src/io/fortran_unformatted.cpp
  src/io/fortran_files.cpp
  src/potential/hard_wall_potential.cpp
  src/potential/dai_light_cu111_potential.cpp
  src/potential/persson_cu111_potential.cpp
  src/potential/potential_factory.cpp
  src/potential/potential_grid.cpp
  src/potential/potential_symmetry.cpp
  src/symmetry/surface_lattice.cpp
  src/symmetry/surface_operation.cpp
  src/symmetry/point_group.cpp
  src/symmetry/grid_symmetry.cpp
  src/gasphase/fourier_grid_hamiltonian.cpp
  src/gasphase/gas_phase.cpp
  src/report/format.cpp
  src/report/sections.cpp
  src/report/propagation_report.cpp
  src/report/analysis_report.cpp
  src/report/analysis_files.cpp
  src/report/fine_analysis_report.cpp
)
target_include_directories(spodvr_core PUBLIC src)
target_compile_definitions(spodvr_core PUBLIC SPODVR_VERSION="${PROJECT_VERSION}")
if(SPODVR_POISON_UNINITIALISED)
  target_compile_definitions(spodvr_core PUBLIC SPODVR_POISON_UNINITIALISED)
endif()
target_link_libraries(spodvr_core
  PUBLIC nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX
  PRIVATE spodvr_crplib spodvr_warnings spodvr_fp)

# The FFTs. The real FFTW must come before oneMKL in the link line: oneMKL also defines the
# fftw_* functions, and the first library that defines a symbol wins.
if(SPODVR_FFT_LIBRARY_USED STREQUAL "FFTW")
  target_link_libraries(spodvr_core PUBLIC PkgConfig::FFTW3)
else()
  target_include_directories(spodvr_core PRIVATE ${MKL_INCLUDE}/fftw ${MKL_INCLUDE})
  target_compile_definitions(spodvr_core PUBLIC SPODVR_FFT_FROM_MKL)
endif()

# BLAS/LAPACK.
if(MKL_FOUND)
  target_link_libraries(spodvr_core PUBLIC MKL::MKL)
  target_compile_definitions(spodvr_core PUBLIC SPODVR_HAVE_MKL SPODVR_MKL_VERSION="${MKL_VERSION}")
  message(STATUS "SPO-DVR: BLAS/LAPACK from oneMKL ${MKL_VERSION} (sequential), FFT library ${SPODVR_FFT_LIBRARY_USED}")
else()
  target_link_libraries(spodvr_core PUBLIC LAPACK::LAPACK)
  message(STATUS "SPO-DVR: BLAS/LAPACK from the system, FFT library ${SPODVR_FFT_LIBRARY_USED}")
endif()

# =============================================================================================
# spodvr_cuda (optional): the GPU backend, the wave packet lives in GPU memory (src/cuda/,
# docs/cuda.md). The kernels (.cu) are compiled by nvcc, the host code (.cpp) by the C++ compiler.
# =============================================================================================
if(SPODVR_USE_CUDA)
  # nvcc from PATH, else from where NVIDIA's packages install it.
  find_program(CMAKE_CUDA_COMPILER nvcc HINTS /usr/local/cuda/bin)
  # The kernels are compiled for the GPU of this machine. For a build that runs elsewhere, set
  # e.g. -DCMAKE_CUDA_ARCHITECTURES="80-real;90-real" (A100, H100).
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
  endif()
  # nvcc compiles the host side of the .cu files with the same compiler as the C++ code (one
  # toolchain; mixed GCC versions break e.g. PGO). -DCMAKE_CUDA_HOST_COMPILER or the environment
  # variable CUDAHOSTCXX override this, e.g. when the CUDA release is too old for that GCC
  # ("unsupported GNU version").
  if(NOT CMAKE_CUDA_HOST_COMPILER AND NOT DEFINED ENV{CUDAHOSTCXX})
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)

  # The cuFFT load callback (src/cuda/fft_callbacks.cu) is not compiled into the library: nvcc
  # turns it into LTO-IR (a fatbin), which is embedded as a byte array and handed to cuFFT at
  # run time (cufftXtSetJITCallback; cuFFT links it with nvJitLink). compute_75 IR runs on every
  # GPU from sm_75 on (cuFFT finalises it for the GPU in use).
  set(callback_dir ${CMAKE_CURRENT_BINARY_DIR}/generated)
  add_custom_command(
    OUTPUT ${callback_dir}/fft_callbacks.fatbin
    COMMAND ${CMAKE_COMMAND} -E make_directory ${callback_dir}
    COMMAND ${CMAKE_CUDA_COMPILER} -ccbin ${CMAKE_CUDA_HOST_COMPILER} -std=c++20 -fmad=false
            -I${CMAKE_CURRENT_SOURCE_DIR}/src -gencode arch=compute_75,code=lto_75 -dc -fatbin
            -o ${callback_dir}/fft_callbacks.fatbin ${CMAKE_CURRENT_SOURCE_DIR}/src/cuda/fft_callbacks.cu
    DEPENDS src/cuda/fft_callbacks.cu src/cuda/fft_callbacks.hpp src/cuda/kernels.hpp src/cuda/device_complex.hpp
    COMMENT "LTO-IR of the cuFFT callback")
  add_custom_command(
    OUTPUT ${callback_dir}/fft_callbacks_fatbin.h
    COMMAND ${CMAKE_COMMAND} -DINPUT=${callback_dir}/fft_callbacks.fatbin
            -DOUTPUT=${callback_dir}/fft_callbacks_fatbin.h -DNAME=spodvr_fft_callbacks_fatbin
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_binary.cmake
    DEPENDS ${callback_dir}/fft_callbacks.fatbin cmake/embed_binary.cmake
    COMMENT "Embedding the cuFFT callback")

  add_library(spodvr_cuda
    src/cuda/cuda_backend.cpp
    src/cuda/cuda_error.cpp
    src/cuda/cuda_fft.cpp
    src/cuda/device_buffer.cpp
    src/cuda/device_kinetic_operator.cpp
    src/cuda/device_potential_operator.cpp
    src/cuda/fft_callbacks.cpp
    ${callback_dir}/fft_callbacks_fatbin.h
    src/cuda/kinetic_kernels.cu
    src/cuda/potential_kernels.cu
    src/cuda/specular_kernels.cu
    src/cuda/wave_function_kernels.cu
  )
  # Only kernels and their launchers are .cu files; all host code of the backend is .cpp. No
  # link-time optimisation (nvcc objects cannot take part in GCC's LTO).
  set_target_properties(spodvr_cuda PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON
                                               INTERPROCEDURAL_OPTIMIZATION OFF)
  if(SPODVR_FAST_MATH)
    # The most aggressive device code (build-cuda-fast). Of --use_fast_math only the FMA
    # contraction matters for our FP64 kernels; the rest (flush-to-zero, fast division and
    # square root) acts on FP32.
    set(cuda_flags --use_fast_math --extra-device-vectorization -Xcompiler=-ffast-math)
  else()
    # No FMA contraction, as on the CPU.
    set(cuda_flags -fmad=false -Xcompiler=-ffp-contract=off,-fno-fast-math)
  endif()
  # Source-line information for Nsight Systems/Compute (no effect on the generated code).
  list(APPEND cuda_flags -lineinfo)
  target_compile_options(spodvr_cuda PRIVATE "$<$<COMPILE_LANGUAGE:CUDA>:${cuda_flags}>")
  if(SPODVR_CUDA_SYNC_DEBUG)
    target_compile_definitions(spodvr_cuda PRIVATE SPODVR_CUDA_SYNC_DEBUG)
  endif()
  target_include_directories(spodvr_cuda PRIVATE ${callback_dir})
  target_link_libraries(spodvr_cuda
    PUBLIC spodvr_core CUDA::cudart CUDA::cufft
    PRIVATE spodvr_warnings spodvr_fp)
  message(STATUS "SPO-DVR: CUDA backend ON (nvcc ${CMAKE_CUDA_COMPILER_VERSION}, host compiler "
                 "${CMAKE_CUDA_HOST_COMPILER}, architectures ${CMAKE_CUDA_ARCHITECTURES})")
endif()

# =============================================================================================
# spodvr_app and the program: command line, run modes, choice of the propagation backend. A
# separate library, so that the optional GPU backend can sit between it and spodvr_core.
# =============================================================================================
add_library(spodvr_app
  src/app/backend_factory.cpp
  src/app/cli.cpp
  src/app/fine_analysis_run.cpp
  src/app/run.cpp
)
target_link_libraries(spodvr_app PUBLIC spodvr_core PRIVATE spodvr_warnings spodvr_fp)
if(SPODVR_USE_CUDA)
  target_link_libraries(spodvr_app PUBLIC spodvr_cuda)
  target_compile_definitions(spodvr_app PUBLIC SPODVR_HAVE_CUDA)
endif()

add_executable(spodvr src/app/main.cpp)
target_link_libraries(spodvr PRIVATE spodvr_app spodvr_warnings spodvr_fp)

# =============================================================================================
# Tests and developer tools
# =============================================================================================
if(SPODVR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
if(SPODVR_BUILD_TOOLS)
  add_subdirectory(tools)
endif()

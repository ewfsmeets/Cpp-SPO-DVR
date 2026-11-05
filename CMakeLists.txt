cmake_minimum_required(VERSION 3.24)
project(SpoDvr VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(SPODVR_WARNINGS_AS_ERRORS "Treat compiler warnings as errors (CI)" OFF)
option(SPODVR_BUILD_TESTS "Build the unit and regression tests" ON)
option(SPODVR_POISON_UNINITIALISED "Debug: fill uninitialised AlignedVector elements with NaN" OFF)
option(SPODVR_BUILD_TOOLS "Build the developer tools (operator timing, wave function dump/compare)" ON)

# ---------------------------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------------------------
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(PkgConfig REQUIRED)

# Intel oneMKL (sequential: the program calls it from inside its own OpenMP threads; LP64:
# 32-bit integers). If found it provides BLAS/LAPACK, replacing the system libraries, and by
# default also the FFTs through its FFTW3 interface. See docs/performance.md.
option(SPODVR_USE_MKL "Use Intel oneMKL when it is found" ON)
set(SPODVR_FFT_LIBRARY "AUTO" CACHE STRING "Implementation of the FFTW3 API: AUTO (oneMKL if found), MKL or FFTW")
set_property(CACHE SPODVR_FFT_LIBRARY PROPERTY STRINGS AUTO MKL FFTW)
if(SPODVR_USE_MKL)
  set(MKL_THREADING sequential)
  set(MKL_INTERFACE lp64)
  set(MKL_LINK dynamic)
  find_package(MKL CONFIG QUIET HINTS $ENV{MKLROOT} /opt/intel/oneapi/mkl/latest)
endif()
set(SPODVR_HAVE_MKL ${MKL_FOUND})  # MKL_FOUND is unset when SPODVR_USE_MKL is OFF

if(SPODVR_FFT_LIBRARY STREQUAL "AUTO")
  if(SPODVR_HAVE_MKL)
    set(SPODVR_FFT_LIBRARY_USED MKL)
  else()
    set(SPODVR_FFT_LIBRARY_USED FFTW)
  endif()
else()
  set(SPODVR_FFT_LIBRARY_USED ${SPODVR_FFT_LIBRARY})
endif()
if(SPODVR_FFT_LIBRARY_USED STREQUAL "MKL" AND NOT SPODVR_HAVE_MKL)
  message(FATAL_ERROR "SPODVR_FFT_LIBRARY=MKL but oneMKL was not found")
endif()
if(SPODVR_FFT_LIBRARY_USED STREQUAL "FFTW")
  pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
endif()
if(NOT SPODVR_HAVE_MKL)
  find_package(LAPACK REQUIRED)
endif()

include(FetchContent)
FetchContent_Declare(nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(nlohmann_json)

# ---------------------------------------------------------------------------------------------
# Warnings
# ---------------------------------------------------------------------------------------------
# (C++ only: these options must not reach nvcc when the CUDA backend is built.)
add_library(spodvr_warnings INTERFACE)
target_compile_options(spodvr_warnings INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Wshadow
                                                 -Wconversion -Wno-sign-conversion>)
if(SPODVR_WARNINGS_AS_ERRORS)
  target_compile_options(spodvr_warnings INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Werror>)
endif()

# Floating-point semantics. No contraction (FMA) and no fast-math, so that results stay
# comparable with the Fortran reference and independent of the target CPU. The other flags do
# not change any result, they only enable optimisations (see docs/performance.md):
#   -fcx-fortran-rules  complex * and / without the C99 NaN-recovery call (__muldc3), i.e. the
#                       gfortran semantics; lets GCC vectorise complex multiplication loops
#   -fno-math-errno     math functions need not set errno (sqrt is inlined)
#   -fno-trapping-math  no floating-point traps are assumed (more if-conversion/vectorisation)
#
# SPODVR_FAST_MATH (experimental, default OFF): -ffast-math instead (FMA contraction, reassociation,
# vectorised math functions). Changes rounding; only for speed experiments, checked with the gates.
option(SPODVR_FAST_MATH "Experimental: compile with -ffast-math (changes rounding)" OFF)
# (C++ only, like the warnings; the CUDA kernels get their own flags, see SPODVR_USE_CUDA.)
add_library(spodvr_fp INTERFACE)
if(SPODVR_FAST_MATH)
  target_compile_options(spodvr_fp INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-ffast-math -fcx-fortran-rules>)
else()
  target_compile_options(spodvr_fp INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-ffp-contract=off -fno-fast-math
                                             -fcx-fortran-rules -fno-math-errno -fno-trapping-math>)
endif()

# Target instruction set. The default x86-64-v3 (AVX2) runs on any recent x86-64 CPU; "native"
# uses everything the build machine has (e.g. AVX-512). With -ffp-contract=off the results of
# our own code do not depend on this choice.
set(SPODVR_ARCH "x86-64-v3" CACHE STRING "Value for -march (e.g. x86-64-v3, x86-64-v4, native; empty: compiler default)")
if(SPODVR_ARCH)
  target_compile_options(spodvr_fp INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-march=${SPODVR_ARCH}>)
endif()

# Profile-guided optimisation (two builds, see docs/performance.md): GENERATE builds an
# instrumented program that writes profiles to SPODVR_PGO_DIR when run; USE rebuilds with them.
set(SPODVR_PGO "OFF" CACHE STRING "Profile-guided optimisation: OFF, GENERATE or USE")
set_property(CACHE SPODVR_PGO PROPERTY STRINGS OFF GENERATE USE)
set(SPODVR_PGO_DIR "${CMAKE_BINARY_DIR}/pgo-profile" CACHE PATH "Directory of the PGO profile data")
if(SPODVR_PGO STREQUAL "GENERATE")
  target_compile_options(spodvr_fp INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-fprofile-generate=${SPODVR_PGO_DIR}
                                             -fprofile-update=atomic>)
  target_link_options(spodvr_fp INTERFACE -fprofile-generate=${SPODVR_PGO_DIR})
elseif(SPODVR_PGO STREQUAL "USE")
  target_compile_options(spodvr_fp INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-fprofile-use=${SPODVR_PGO_DIR}
                                             -fprofile-correction -Wno-missing-profile>)
endif()

# Link-time optimisation for optimised builds.
include(CheckIPOSupported)
check_ipo_supported(RESULT spodvr_ipo_supported OUTPUT spodvr_ipo_output)
if(spodvr_ipo_supported)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

# ---------------------------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------------------------
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
  PRIVATE spodvr_warnings spodvr_fp)
# The real FFTW must come before oneMKL in the link line: oneMKL also exports the fftw_*
# functions, and the first library that defines a symbol wins.
if(SPODVR_FFT_LIBRARY_USED STREQUAL "FFTW")
  target_link_libraries(spodvr_core PUBLIC PkgConfig::FFTW3)
else()
  target_include_directories(spodvr_core PRIVATE ${MKL_INCLUDE}/fftw ${MKL_INCLUDE})
  target_compile_definitions(spodvr_core PUBLIC SPODVR_FFT_FROM_MKL)
endif()
if(SPODVR_HAVE_MKL)
  target_link_libraries(spodvr_core PUBLIC MKL::MKL)
  target_compile_definitions(spodvr_core PUBLIC SPODVR_HAVE_MKL SPODVR_MKL_VERSION="${MKL_VERSION}")
else()
  target_link_libraries(spodvr_core PUBLIC LAPACK::LAPACK)
endif()
if(SPODVR_HAVE_MKL)
  message(STATUS "SPO-DVR: BLAS/LAPACK from oneMKL ${MKL_VERSION} (sequential), FFT library ${SPODVR_FFT_LIBRARY_USED}")
else()
  message(STATUS "SPO-DVR: BLAS/LAPACK from the system, FFT library ${SPODVR_FFT_LIBRARY_USED}")
endif()

# ---------------------------------------------------------------------------------------------
# Optional GPU backend (CUDA): the wave packet lives in GPU memory (src/cuda/, docs/cuda.md).
# Kernels (.cu) are compiled by nvcc; the host code of the backend (.cpp) by the C++ compiler.
# ---------------------------------------------------------------------------------------------
option(SPODVR_USE_CUDA "Build the CUDA propagation backend (--device cuda)" OFF)
option(SPODVR_CUDA_SYNC_DEBUG "Debug: wait for and check every CUDA kernel right after its launch" OFF)
if(SPODVR_USE_CUDA)
  # nvcc compiles the host side of the .cu files with a host compiler it supports: CUDA 12.0
  # rejects GCC > 12 (crt/host_config.h), so use g++-12 when it is installed.
  if(NOT CMAKE_CUDA_HOST_COMPILER AND EXISTS /usr/bin/g++-12)
    set(CMAKE_CUDA_HOST_COMPILER /usr/bin/g++-12)
  endif()
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 89-real)  # RTX 40xx; e.g. "80-real;90-real" for A100/H100
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)

  add_library(spodvr_cuda
    src/cuda/cuda_backend.cpp
    src/cuda/cuda_error.cpp
    src/cuda/cuda_fft.cpp
    src/cuda/device_buffer.cpp
    src/cuda/device_kinetic_operator.cpp
    src/cuda/device_potential_operator.cpp
    src/cuda/kinetic_kernels.cu
    src/cuda/potential_kernels.cu
    src/cuda/specular_kernels.cu
    src/cuda/wave_function_kernels.cu
  )
  # The kernels only need C++17 (nvcc 12.0 does not handle every C++20 construct of the host code,
  # which is why the backend's host code is in .cpp files). No FMA contraction, as on the CPU.
  set_target_properties(spodvr_cuda PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON
                                               INTERPROCEDURAL_OPTIMIZATION OFF)
  target_compile_options(spodvr_cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-fmad=false
                                             -Xcompiler=-ffp-contract=off,-fno-fast-math>)
  if(SPODVR_CUDA_SYNC_DEBUG)
    target_compile_definitions(spodvr_cuda PRIVATE SPODVR_CUDA_SYNC_DEBUG)
  endif()
  target_link_libraries(spodvr_cuda PUBLIC spodvr_core CUDA::cudart CUDA::cufft PRIVATE spodvr_warnings spodvr_fp)
  message(STATUS "SPO-DVR: CUDA backend ON (nvcc ${CMAKE_CUDA_COMPILER_VERSION}, host compiler "
                 "${CMAKE_CUDA_HOST_COMPILER}, architectures ${CMAKE_CUDA_ARCHITECTURES})")
endif()

# ---------------------------------------------------------------------------------------------
# Program driver (command line, run modes, choice of the propagation backend). A separate library
# so that an optional GPU backend can sit between it and spodvr_core.
# ---------------------------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------------------------
if(SPODVR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
if(SPODVR_BUILD_TOOLS)
  add_subdirectory(tools)
endif()

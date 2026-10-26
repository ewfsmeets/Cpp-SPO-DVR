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
add_library(spodvr_warnings INTERFACE)
target_compile_options(spodvr_warnings INTERFACE -Wall -Wextra -Wpedantic -Wshadow -Wconversion
                                                 -Wno-sign-conversion)
if(SPODVR_WARNINGS_AS_ERRORS)
  target_compile_options(spodvr_warnings INTERFACE -Werror)
endif()

# Floating-point semantics. No contraction (FMA) and no fast-math, so that results stay
# comparable with the Fortran reference and independent of the target CPU. The other flags do
# not change any result, they only enable optimisations (see docs/performance.md):
#   -fcx-fortran-rules  complex * and / without the C99 NaN-recovery call (__muldc3), i.e. the
#                       gfortran semantics; lets GCC vectorise complex multiplication loops
#   -fno-math-errno     math functions need not set errno (sqrt is inlined)
#   -fno-trapping-math  no floating-point traps are assumed (more if-conversion/vectorisation)
add_library(spodvr_fp INTERFACE)
target_compile_options(spodvr_fp INTERFACE -ffp-contract=off -fno-fast-math -fcx-fortran-rules
                                           -fno-math-errno -fno-trapping-math)

# Target instruction set. The default x86-64-v3 (AVX2) runs on any recent x86-64 CPU; "native"
# uses everything the build machine has (e.g. AVX-512). With -ffp-contract=off the results of
# our own code do not depend on this choice.
set(SPODVR_ARCH "x86-64-v3" CACHE STRING "Value for -march (e.g. x86-64-v3, x86-64-v4, native; empty: compiler default)")
if(SPODVR_ARCH)
  target_compile_options(spodvr_fp INTERFACE -march=${SPODVR_ARCH})
endif()

# Profile-guided optimisation (two builds, see docs/performance.md): GENERATE builds an
# instrumented program that writes profiles to SPODVR_PGO_DIR when run; USE rebuilds with them.
set(SPODVR_PGO "OFF" CACHE STRING "Profile-guided optimisation: OFF, GENERATE or USE")
set_property(CACHE SPODVR_PGO PROPERTY STRINGS OFF GENERATE USE)
set(SPODVR_PGO_DIR "${CMAKE_BINARY_DIR}/pgo-profile" CACHE PATH "Directory of the PGO profile data")
if(SPODVR_PGO STREQUAL "GENERATE")
  target_compile_options(spodvr_fp INTERFACE -fprofile-generate=${SPODVR_PGO_DIR} -fprofile-update=atomic)
  target_link_options(spodvr_fp INTERFACE -fprofile-generate=${SPODVR_PGO_DIR})
elseif(SPODVR_PGO STREQUAL "USE")
  target_compile_options(spodvr_fp INTERFACE -fprofile-use=${SPODVR_PGO_DIR} -fprofile-correction
                                             -Wno-missing-profile)
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
  src/core/machine_epsilon.cpp
  src/core/numerical_libraries.cpp
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
  src/propagation/averages.cpp
  src/propagation/propagation.cpp
  src/analysis/analysis_basis.cpp
  src/analysis/balint_kurti.cpp
  src/analysis/s_matrix.cpp
  src/io/checkpoint.cpp
  src/potential/hard_wall_potential.cpp
  src/potential/potential_factory.cpp
  src/potential/potential_grid.cpp
  src/gasphase/fourier_grid_hamiltonian.cpp
  src/gasphase/gas_phase.cpp
  src/report/format.cpp
  src/report/sections.cpp
  src/report/propagation_report.cpp
  src/report/analysis_report.cpp
  src/report/analysis_files.cpp
  src/app/cli.cpp
  src/app/run.cpp
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

add_executable(spodvr src/app/main.cpp)
target_link_libraries(spodvr PRIVATE spodvr_core spodvr_warnings spodvr_fp)

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

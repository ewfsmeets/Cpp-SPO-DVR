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

# ---------------------------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------------------------
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
find_package(LAPACK REQUIRED)

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

# No floating-point contraction (FMA) so that results stay comparable with the Fortran
# reference and independent of the target CPU.
add_library(spodvr_fp INTERFACE)
target_compile_options(spodvr_fp INTERFACE -ffp-contract=off -fno-fast-math)

# ---------------------------------------------------------------------------------------------
# Core library
# ---------------------------------------------------------------------------------------------
add_library(spodvr_core
  src/core/machine_epsilon.cpp
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
  src/potential/hard_wall_potential.cpp
  src/potential/potential_factory.cpp
  src/potential/potential_grid.cpp
  src/gasphase/fourier_grid_hamiltonian.cpp
  src/gasphase/gas_phase.cpp
  src/report/format.cpp
  src/report/sections.cpp
  src/app/cli.cpp
  src/app/run.cpp
)
target_include_directories(spodvr_core PUBLIC src)
target_compile_definitions(spodvr_core PUBLIC SPODVR_VERSION="${PROJECT_VERSION}")
target_link_libraries(spodvr_core
  PUBLIC nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX
  PRIVATE PkgConfig::FFTW3 LAPACK::LAPACK spodvr_warnings spodvr_fp)

add_executable(spodvr src/app/main.cpp)
target_link_libraries(spodvr PRIVATE spodvr_core spodvr_warnings spodvr_fp)

# ---------------------------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------------------------
if(SPODVR_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

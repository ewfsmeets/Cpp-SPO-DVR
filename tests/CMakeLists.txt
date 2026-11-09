# One test program per tests/test_*.cpp file (plus tests/cuda/test_*.cpp with the GPU backend;
# those skip when no GPU is present). ctest runs every GoogleTest case as its own test.
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
  GIT_SHALLOW TRUE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

file(GLOB test_sources CONFIGURE_DEPENDS test_*.cpp)
if(SPODVR_USE_CUDA)
  file(GLOB cuda_test_sources CONFIGURE_DEPENDS cuda/test_*.cpp)
  list(APPEND test_sources ${cuda_test_sources})
endif()

foreach(source IN LISTS test_sources)
  get_filename_component(name ${source} NAME_WE)
  add_executable(${name} ${source})
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})  # shared test headers
  target_link_libraries(${name} PRIVATE spodvr_app GTest::gtest_main spodvr_warnings spodvr_fp)
  # Fortran reference output of the regression tests.
  target_compile_definitions(${name} PRIVATE SPODVR_REFERENCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/reference")
  gtest_discover_tests(${name} DISCOVERY_TIMEOUT 60)
endforeach()

# The CRP library's own tests (they link only spodvr_crplib).
add_subdirectory(crplib)

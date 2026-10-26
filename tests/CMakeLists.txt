FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG v1.14.0
  GIT_SHALLOW TRUE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

set(SPODVR_REFERENCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/reference")

function(spodvr_add_test name)
  add_executable(${name} ${name}.cpp ${ARGN})
  target_link_libraries(${name} PRIVATE spodvr_core GTest::gtest_main spodvr_warnings spodvr_fp)
  target_compile_definitions(${name} PRIVATE SPODVR_REFERENCE_DIR="${SPODVR_REFERENCE_DIR}")
  gtest_discover_tests(${name} DISCOVERY_TIMEOUT 60)
endfunction()

file(GLOB SPODVR_UNIT_TESTS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp)
foreach(test_source ${SPODVR_UNIT_TESTS})
  get_filename_component(test_name ${test_source} NAME_WE)
  spodvr_add_test(${test_name})
endforeach()

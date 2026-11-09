# Tests of the CRP library (src/crplib/): one program per tests/crplib/test_*.cpp file. They link
# only spodvr_crplib (and GoogleTest); ctest -R CrpLib runs them. Reference data of the Fortran
# CRP library (tools/fortran_reference/run_crp_reference.sh) is in tests/crplib/reference/.
file(GLOB crplib_test_sources CONFIGURE_DEPENDS test_*.cpp)
foreach(source IN LISTS crplib_test_sources)
  get_filename_component(name ${source} NAME_WE)
  set(target crplib_${name})
  add_executable(${target} ${source})
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})  # shared test headers
  target_link_libraries(${target} PRIVATE spodvr_crplib GTest::gtest_main spodvr_warnings spodvr_crplib_fp)
  target_compile_definitions(${target} PRIVATE CRPLIB_REFERENCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/reference")
  gtest_discover_tests(${target} TEST_PREFIX "CrpLib." DISCOVERY_TIMEOUT 60)
endforeach()

# The library includes only the standard library and its own headers.
add_test(NAME CrpLib.IncludeBoundary
         COMMAND ${CMAKE_COMMAND} -DDIR=${PROJECT_SOURCE_DIR}/src/crplib -P ${CMAKE_CURRENT_SOURCE_DIR}/check_include_boundary.cmake)

# The synthetic reference cases (data set + point lists) for tools/fortran_reference/run_crp_reference.sh.
add_executable(crplib_make_dataset make_crp_dataset.cpp)
target_link_libraries(crplib_make_dataset PRIVATE spodvr_crplib spodvr_warnings spodvr_crplib_fp)

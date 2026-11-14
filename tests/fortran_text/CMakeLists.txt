# Tests of the Fortran text reader (src/fortran_text/): one program per tests/fortran_text/test_*.cpp
# file. They link only spodvr_fortran_text (and GoogleTest); ctest -R FortranText runs them. Its
# include rule (own headers and the standard library only) is checked by the CTest
# CrpLib.IncludeBoundary (tests/crplib), which ctest -R FortranText does not run.
file(GLOB fortran_text_test_sources CONFIGURE_DEPENDS test_*.cpp)
foreach(source IN LISTS fortran_text_test_sources)
  get_filename_component(name ${source} NAME_WE)
  set(target fortran_text_${name})
  add_executable(${target} ${source})
  # Shared test headers: tests/temporary_directory.hpp.
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
  target_link_libraries(${target} PRIVATE spodvr_fortran_text GTest::gtest_main spodvr_warnings spodvr_crplib_fp)
  gtest_discover_tests(${target} TEST_PREFIX "FortranText." DISCOVERY_TIMEOUT 60)
endforeach()

# Tests of the SCM distorted-surface correction against the Fortran reference data in
# tests/scm/reference/ (tools/fortran_reference/run_scm_reference.sh): one program per
# tests/scm/test_*.cpp file, linked against spodvr_core; ctest -R Scm runs them (prefix "Scm.").
file(GLOB scm_test_sources CONFIGURE_DEPENDS test_*.cpp)
foreach(source IN LISTS scm_test_sources)
  get_filename_component(name ${source} NAME_WE)
  set(target scm_${name})
  add_executable(${target} ${source})
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
  target_link_libraries(${target} PRIVATE spodvr_core GTest::gtest_main spodvr_warnings spodvr_fp)
  target_compile_definitions(${target} PRIVATE SPODVR_SCM_REFERENCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/reference"
                             SPODVR_SCM_EXAMPLE_DIR="${PROJECT_SOURCE_DIR}/potentials/H2Cu111SRP_C3v.dir")
  gtest_discover_tests(${target} TEST_PREFIX "Scm." DISCOVERY_TIMEOUT 60)
endforeach()

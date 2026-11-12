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

# The CRP adapter test writes synthetic data sets with the CRP library's test helpers.
target_include_directories(test_crp_potential PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/crplib)
target_link_libraries(test_crp_potential PRIVATE spodvr_crplib)

# The flux analysis is strict in every build (spodvr_strict_objects); its tests compare it bit for
# bit with expressions written out in the test, which must then be strict too (also in build-fast).
target_link_libraries(test_flux_analysis PRIVATE spodvr_strict_fp)
# The same for the time loop's flux tests, which recompute the sums in the test (flux_propagation_run.hpp).
target_link_libraries(test_flux_propagation PRIVATE spodvr_strict_fp)
if(SPODVR_USE_CUDA)
  target_link_libraries(test_cuda_flux PRIVATE spodvr_strict_fp)
endif()

# The docs against the code (tools/check_docs.py; README "Testing"), with spodvr --help for the
# command-line options, and negative tests that show the checker fails on a planted fault.
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(check_docs ${PROJECT_SOURCE_DIR}/tools/check_docs.py)
add_test(NAME Docs.Consistency
         COMMAND ${Python3_EXECUTABLE} -I -B ${check_docs} --root ${PROJECT_SOURCE_DIR} --spodvr $<TARGET_FILE:spodvr>)
set(check_docs_negative ${CMAKE_COMMAND} -DPYTHON=${Python3_EXECUTABLE} -DCHECKER=${check_docs}
    -DROOT=${PROJECT_SOURCE_DIR})
set(check_docs_negative_script ${PROJECT_SOURCE_DIR}/tools/check_docs_negative_test.cmake)
add_test(NAME Docs.UndocumentedKeyFails
         COMMAND ${check_docs_negative} -DCASE=undocumented-key
                 -DWORK=${CMAKE_CURRENT_BINARY_DIR}/check_docs_negative/undocumented-key -P ${check_docs_negative_script})
add_test(NAME Docs.StaleAllowlistRowFails
         COMMAND ${check_docs_negative} -DCASE=stale-allowlist
                 -DWORK=${CMAKE_CURRENT_BINARY_DIR}/check_docs_negative/stale-allowlist -P ${check_docs_negative_script})
add_test(NAME Docs.UpstreamPinMismatchFails
         COMMAND ${check_docs_negative} -DCASE=upstream-pin
                 -DWORK=${CMAKE_CURRENT_BINARY_DIR}/check_docs_negative/upstream-pin -P ${check_docs_negative_script})

# The refusals of the Fortran reference tools (tools/fortran_reference/test_reference_tools.sh,
# tools/fortran_reference/README.md): stub installs and synthetic build logs, no Fortran run.
set(reference_tools_test ${PROJECT_SOURCE_DIR}/tools/fortran_reference/test_reference_tools.sh)
function(add_reference_tools_test name test_case)
  add_test(NAME FortranReference.${name}
           COMMAND bash ${reference_tools_test} ${test_case}
                   ${CMAKE_CURRENT_BINARY_DIR}/fortran_reference_tools/${test_case} ${CMAKE_CXX_COMPILER})
endfunction()
add_reference_tools_test(RetiredVariablesRefused retired-variables)
add_reference_tools_test(ExistingScratchDirRefused existing-scratch)
add_reference_tools_test(ChangedInstallRefused changed-install)
add_reference_tools_test(RuntimeDriftRefused runtime-drift)
add_reference_tools_test(RunChecksExitAndOutputs run)
add_reference_tools_test(BuildLogAllowLists check-log)
add_reference_tools_test(InstructionCheck check-binary)
add_reference_tools_test(ManifestRecordAndDelete manifest)
add_reference_tools_test(RegenerateRefusesBeforeDeleting regenerate-refusals)

# refdata.py compare (tools/fortran_reference/test_refdata.py): expected/r1713.tsv is valid and
# matches the committed files, every class it allows passes, each fault and each refusal fails.
set(refdata_test ${PROJECT_SOURCE_DIR}/tools/fortran_reference/test_refdata.py)
function(add_refdata_test name test_case)
  add_test(NAME FortranReference.Compare${name}
           COMMAND ${Python3_EXECUTABLE} -I -B ${refdata_test} ${test_case}
                   ${CMAKE_CURRENT_BINARY_DIR}/refdata_compare/${test_case})
endfunction()
add_refdata_test(ExpectedFileMatchesTree expected-file)
add_refdata_test(AllowedDifferencesPass allowed)
add_refdata_test(EditedNumberFails edited-number)
add_refdata_test(DeletedFileFails deleted-file)
add_refdata_test(UnknownLabelFails unknown-label)
add_refdata_test(SMagnitudeChangeFails s-magnitude)
add_refdata_test(BkBeyondLimitFails bk-limit)
add_refdata_test(BkNotFiniteFails bk-nan)
add_refdata_test(PrintedNotFiniteFails non-finite-number)
add_refdata_test(NumberCountChangeFails number-count)
add_refdata_test(ProvenanceFaultsFail provenance)
add_refdata_test(OtherInstallFails install-mismatch)
add_refdata_test(FineLineWithoutBkInputFails fine-without-fort22)
add_refdata_test(OverlappingRowsRefused multi-match)
add_refdata_test(BadExpectedFileRefused refusals)

# spodvr_compare_wave (tools/compare/test_compare_wave.py): finite files compare as before, a NaN or
# inf in either file always fails.
if(SPODVR_BUILD_TOOLS)
  set(compare_wave_test ${PROJECT_SOURCE_DIR}/tools/compare/test_compare_wave.py)
  add_test(NAME CompareWave.FiniteFiles
           COMMAND ${Python3_EXECUTABLE} -I -B ${compare_wave_test} $<TARGET_FILE:spodvr_compare_wave> finite
                   ${CMAKE_CURRENT_BINARY_DIR}/compare_wave/finite)
  add_test(NAME CompareWave.NonFiniteValueFails
           COMMAND ${Python3_EXECUTABLE} -I -B ${compare_wave_test} $<TARGET_FILE:spodvr_compare_wave> non-finite
                   ${CMAKE_CURRENT_BINARY_DIR}/compare_wave/non-finite)
endif()

# Every Fortran reference file listed in a MANIFEST.tsv exists with its recorded sha256, so a
# missing or edited reference fails instead of skipping its gate (tools/fortran_reference/README.md).
set(reference_manifest_check ${PROJECT_SOURCE_DIR}/tools/fortran_reference/ReferenceManifest.cmake)
add_test(NAME ReferenceManifest.FilesMatch
         COMMAND ${CMAKE_COMMAND} -DROOT=${PROJECT_SOURCE_DIR} -P ${reference_manifest_check})
add_test(NAME ReferenceManifest.ChangedOrMissingFileFails
         COMMAND ${CMAKE_COMMAND} -DROOT=${PROJECT_SOURCE_DIR}
                 -DWORK=${CMAKE_CURRENT_BINARY_DIR}/reference_manifest_negative -P ${reference_manifest_check})

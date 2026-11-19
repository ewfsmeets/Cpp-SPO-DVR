# One test program per tests/test_*.cpp file (plus tests/cuda/test_*.cpp with the GPU backend;
# those skip when no GPU is present). ctest runs every GoogleTest case as its own test, except the
# suites with a shared set-up registered below the loop (one ctest per suite).
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

# Test programs whose ctest tests are registered below the loop.
set(tests_registered_below test_restart test_regression_pes test_plot_fixtures)

foreach(source IN LISTS test_sources)
  get_filename_component(name ${source} NAME_WE)
  add_executable(${name} ${source})
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})  # shared test headers
  target_link_libraries(${name} PRIVATE spodvr_app GTest::gtest_main spodvr_warnings spodvr_fp)
  # Fortran reference output of the regression tests; the repository root (designer_cases.hpp).
  target_compile_definitions(${name} PRIVATE SPODVR_REFERENCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/reference"
                                             SPODVR_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  if(NOT name IN_LIST tests_registered_below)
    gtest_discover_tests(${name} DISCOVERY_TIMEOUT 60)
  endif()
endforeach()

# Suites whose cases share an expensive set-up (a function-local static in the test file: program runs,
# potential grids) run as one ctest in one process, so the set-up is computed once instead of once per case
# (docs/CONTRIBUTING.md "Selecting tests"). Such a test fails when a case skips (ctest's
# SKIP_REGULAR_EXPRESSION would mark it skipped even when another case failed: it ignores the exit code) or
# when its filter selects nothing (gtest passes an empty run). Disabled cases run, as gtest_discover_tests
# runs them.
function(add_one_process_test name target filter)
  add_test(NAME ${name} COMMAND ${target} --gtest_filter=${filter} --gtest_also_run_disabled_tests)
  set_tests_properties(${name} PROPERTIES FAIL_REGULAR_EXPRESSION "\\[  SKIPPED \\];Running 0 tests from")
endfunction()
add_one_process_test(Restart test_restart "Restart.*")
add_one_process_test(FineAnalysisRuns test_restart "FineAnalysisRuns.*")
# The other suites of test_restart stay one ctest per case (RestartFiles skips one case when run as root).
gtest_discover_tests(test_restart DISCOVERY_TIMEOUT 60 TEST_FILTER "-Restart.*:FineAnalysisRuns.*")
# The negative tests on their own: one of them is a known failure of the fast-math build
# (docs/performance.md), which must not hide the gates of the other suites there.
add_one_process_test(RegressionPes test_regression_pes "-RegressionPesScmNegative.*")
add_one_process_test(RegressionPesScmNegative test_regression_pes "RegressionPesScmNegative.*")
# The runs that the plot tool's tests read (tests/test_plot_fixtures.cpp): one process, the fixture of
# every PlotRuns.* pytest test, named PlotRuns so that the quick tier leaves it out (README "Testing").
set(plot_fixtures_dir ${CMAKE_CURRENT_BINARY_DIR}/plot_fixtures)
target_compile_definitions(test_plot_fixtures PRIVATE SPODVR_PLOT_FIXTURES_DIR="${plot_fixtures_dir}")
add_one_process_test(PlotRuns.Fixtures test_plot_fixtures "PlotRuns.*")
set_tests_properties(PlotRuns.Fixtures PROPERTIES FIXTURES_SETUP plot_fixtures)

# The CRP library's own tests (they link only spodvr_crplib, with its spodvr_fortran_text), and
# those of the Fortran text reader (they link only spodvr_fortran_text).
add_subdirectory(crplib)
add_subdirectory(fortran_text)
# The SCM correction against its Fortran reference data (tests/scm/reference/).
add_subdirectory(scm)

# The CRP adapter test writes synthetic data sets with the CRP library's test helpers.
target_include_directories(test_crp_potential PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/crplib)
target_link_libraries(test_crp_potential PRIVATE spodvr_crplib)

# The SCM file readers read the example files of the repository's SRP48 data set.
target_compile_definitions(test_scm_io PRIVATE SPODVR_SCM_EXAMPLE_DIR="${PROJECT_SOURCE_DIR}/potentials/H2Cu111SRP_C3v.dir")
# The SCM potential test computes expected values bit for bit and tests for NaN: strict floating
# point in every build (spodvr_strict_fp after the loop's spodvr_fp overrides -ffast-math).
target_link_libraries(test_scm_potential PRIVATE spodvr_strict_fp)
# So do the SCM diagnostics tests (their stub potentials' values, NaN, bitwise comparisons).
target_link_libraries(test_scm_diagnostics PRIVATE spodvr_strict_fp)

# The docs against the code (tools/check_docs.py; docs/CONTRIBUTING.md "Testing"), with spodvr --help for the
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
add_test(NAME Docs.SubscriptAnchorFails
         COMMAND ${check_docs_negative} -DCASE=subscript-anchor
                 -DWORK=${CMAKE_CURRENT_BINARY_DIR}/check_docs_negative/subscript-anchor -P ${check_docs_negative_script})

# spodvr --version names the git commit the program was built from (cmake/git_describe.cmake): in a git
# checkout a commit, never "unknown" (the tests' run_program calls always get "unknown").
add_test(NAME Program.VersionNamesTheCommit COMMAND $<TARGET_FILE:spodvr> --version)
set_tests_properties(Program.VersionNamesTheCommit PROPERTIES
                     PASS_REGULAR_EXPRESSION "^spodvr [0-9.]+ \\(git [0-9A-Za-z._-]+\\)")
if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
  set_tests_properties(Program.VersionNamesTheCommit PROPERTIES FAIL_REGULAR_EXPRESSION "git unknown")
endif()

# The refusals and the input staging of the Fortran reference tools (test_reference_tools.sh,
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
add_reference_tools_test(ConfigHelperRefusals config-refusals)
add_reference_tools_test(BuildRefusals build-refusals)
add_reference_tools_test(ManifestRecordAndDelete manifest)
add_reference_tools_test(RegenerateRefusesBeforeDeleting regenerate-refusals)
add_reference_tools_test(ScmInputsStagedUnderShortNames stage-scm)
add_reference_tools_test(StagingOfTheCasesUnchanged staged-cases)
add_reference_tools_test(RunScriptRefusals run-script-refusals)
add_reference_tools_test(ModifiedTargetsRefused unmodified-targets)

# The Python tests (pytest, in the repository's Python environment: Python3_EXECUTABLE, the root
# CMakeLists.txt). One ctest per selection, a test module or one test function, with the tool's own
# pytest.ini and its directory as the root; label python (ctest -L python). pytest recreates the
# basetemp of each ctest. Every selection is also listed in pytest_selections.txt, which
# Python.EveryTestRegistered compares with what pytest collects (docs/CONTRIBUTING.md "Testing").
set(pytest_selections ${CMAKE_BINARY_DIR}/pytest_selections.txt)
file(WRITE ${pytest_selections} "# tool<TAB>selection<TAB>built|not-built, written by tests/CMakeLists.txt\n")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/pytest)
# add_pytest_test(NAME <ctest name> TOOL <tool directory> SELECT <module or module::test> [ARGS …])
function(add_pytest_test)
  cmake_parse_arguments(PARSE_ARGV 0 arg "" "NAME;TOOL;SELECT" "ARGS")
  set(tool ${PROJECT_SOURCE_DIR}/${arg_TOOL})
  add_test(NAME ${arg_NAME}
           COMMAND ${Python3_EXECUTABLE} -I -B -m pytest -c ${tool}/pytest.ini --rootdir=${tool}
                   -p no:cacheprovider -q --basetemp=${CMAKE_CURRENT_BINARY_DIR}/pytest/${arg_NAME}
                   ${tool}/${arg_SELECT} ${arg_ARGS})
  set_tests_properties(${arg_NAME} PROPERTIES
                       LABELS python
                       ENVIRONMENT "MPLBACKEND=Agg;MPLCONFIGDIR=${CMAKE_CURRENT_BINARY_DIR}/pytest/matplotlib")
  file(APPEND ${pytest_selections} "${arg_TOOL}\t${arg_SELECT}\tbuilt\n")
endfunction()
# A selection whose ctest exists only in some configurations, listed when it is not registered.
function(note_pytest_not_built tool select)
  file(APPEND ${pytest_selections} "${tool}\t${select}\tnot-built\n")
endfunction()

# extract_pes.py on synthetic .pes files (tools/fortran_reference/test_extract_pes.py): files without
# SCM records give the output of before, the six SCM records are detected, wrong files refused.
add_pytest_test(NAME FortranReference.ExtractPesWithoutScmUnchanged TOOL tools/fortran_reference
                SELECT test_extract_pes.py::test_without_scm)
add_pytest_test(NAME FortranReference.ExtractPesScmRecords TOOL tools/fortran_reference
                SELECT test_extract_pes.py::test_scm)
add_pytest_test(NAME FortranReference.ExtractPesRefusals TOOL tools/fortran_reference
                SELECT test_extract_pes.py::test_refusals)

# refdata.py compare (tools/fortran_reference/test_refdata.py): expected/r1713.tsv is valid and
# matches the committed files, every class it allows passes, each fault and each refusal fails. The
# copies of refdata.py's MANIFEST header, reference trees and NOT_A_CASE (ReferenceManifest.cmake,
# tools/check_docs.py, reference_config.sh, .pre-commit-config.yaml) agree; planted differences fail.
add_pytest_test(NAME FortranReference.CompareExpectedFileMatchesTree TOOL tools/fortran_reference
                SELECT test_refdata.py::test_expected_file)
add_pytest_test(NAME FortranReference.CompareAllowedDifferencesPass TOOL tools/fortran_reference
                SELECT test_refdata.py::test_allowed)
add_pytest_test(NAME FortranReference.CompareEditedNumberFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_edited_number)
add_pytest_test(NAME FortranReference.CompareDeletedFileFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_deleted_file)
add_pytest_test(NAME FortranReference.CompareUnknownLabelFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_unknown_label)
add_pytest_test(NAME FortranReference.CompareSMagnitudeChangeFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_s_magnitude)
add_pytest_test(NAME FortranReference.CompareBkBeyondLimitFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_bk_limit)
add_pytest_test(NAME FortranReference.CompareBkNotFiniteFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_bk_nan)
add_pytest_test(NAME FortranReference.ComparePrintedNotFiniteFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_non_finite_number)
add_pytest_test(NAME FortranReference.CompareNumberCountChangeFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_number_count)
add_pytest_test(NAME FortranReference.CompareProvenanceFaultsFail TOOL tools/fortran_reference
                SELECT test_refdata.py::test_provenance)
add_pytest_test(NAME FortranReference.CompareOtherInstallFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_install_mismatch)
add_pytest_test(NAME FortranReference.CompareFineLineWithoutBkInputFails TOOL tools/fortran_reference
                SELECT test_refdata.py::test_fine_without_fort22)
add_pytest_test(NAME FortranReference.CompareOverlappingRowsRefused TOOL tools/fortran_reference
                SELECT test_refdata.py::test_multi_match)
add_pytest_test(NAME FortranReference.CompareBadExpectedFileRefused TOOL tools/fortran_reference
                SELECT test_refdata.py::test_refusals)
add_pytest_test(NAME FortranReference.ReferenceListsAgree TOOL tools/fortran_reference
                SELECT test_refdata.py::test_reference_lists)

# spodvr_compare_wave (tools/compare/test_compare_wave.py): finite files compare as before, a NaN or
# inf in either file always fails, bad arguments are usage errors and unreadable files are refused.
if(SPODVR_BUILD_TOOLS)
  add_pytest_test(NAME CompareWave.FiniteFiles TOOL tools/compare SELECT test_compare_wave.py::test_finite
                  ARGS --compare-wave-tool=$<TARGET_FILE:spodvr_compare_wave>)
  add_pytest_test(NAME CompareWave.NonFiniteValueFails TOOL tools/compare SELECT test_compare_wave.py::test_non_finite
                  ARGS --compare-wave-tool=$<TARGET_FILE:spodvr_compare_wave>)
  add_pytest_test(NAME CompareWave.BadArgumentsAndFilesRefused TOOL tools/compare
                  SELECT test_compare_wave.py::test_arguments
                  ARGS --compare-wave-tool=$<TARGET_FILE:spodvr_compare_wave>)
else()
  note_pytest_not_built(tools/compare test_compare_wave.py::test_finite)
  note_pytest_not_built(tools/compare test_compare_wave.py::test_non_finite)
  note_pytest_not_built(tools/compare test_compare_wave.py::test_arguments)
endif()

# The averages comparison of tools/compare_propagation.sh (tools/compare/test_compare_averages.py):
# equal files pass, a NaN or inf in either file always fails, a file without rows is refused, a bad
# tolerance is a usage error.
add_pytest_test(NAME ComparePropagation.EqualFilesPass TOOL tools/compare SELECT test_compare_averages.py::test_equal)
add_pytest_test(NAME ComparePropagation.NonFiniteValueFails TOOL tools/compare
                SELECT test_compare_averages.py::test_non_finite)
add_pytest_test(NAME ComparePropagation.FileWithoutRowsRefused TOOL tools/compare
                SELECT test_compare_averages.py::test_rows)
add_pytest_test(NAME ComparePropagation.BadToleranceIsUsageError TOOL tools/compare
                SELECT test_compare_averages.py::test_arguments)

# The wave-packet designer (tools/wave_packet_designer, numpy and matplotlib): one ctest per test module.
# They compare the designer's copies of the program's rules with the C++ (key table, check list, grid
# indices, packet values) and include negative cases.
add_pytest_test(NAME WavePacketDesigner.FftSizes TOOL tools/wave_packet_designer SELECT tests/test_fft_sizes.py)
add_pytest_test(NAME WavePacketDesigner.RealCompare TOOL tools/wave_packet_designer SELECT tests/test_real_compare.py)
add_pytest_test(NAME WavePacketDesigner.InputReading TOOL tools/wave_packet_designer SELECT tests/test_input_reading.py)
add_pytest_test(NAME WavePacketDesigner.Validate TOOL tools/wave_packet_designer SELECT tests/test_validate.py)
add_pytest_test(NAME WavePacketDesigner.Absorbers TOOL tools/wave_packet_designer SELECT tests/test_absorbers.py)
add_pytest_test(NAME WavePacketDesigner.Grids TOOL tools/wave_packet_designer SELECT tests/test_grids.py)
add_pytest_test(NAME WavePacketDesigner.Packet TOOL tools/wave_packet_designer SELECT tests/test_packet.py)
add_pytest_test(NAME WavePacketDesigner.Memory TOOL tools/wave_packet_designer SELECT tests/test_memory.py)
add_pytest_test(NAME WavePacketDesigner.ExpectedCases TOOL tools/wave_packet_designer
                SELECT tests/test_expected_cases.py)
add_pytest_test(NAME WavePacketDesigner.GasPhase TOOL tools/wave_packet_designer SELECT tests/test_gas_phase.py)
add_pytest_test(NAME WavePacketDesigner.Design TOOL tools/wave_packet_designer SELECT tests/test_design.py)
add_pytest_test(NAME WavePacketDesigner.Plots TOOL tools/wave_packet_designer SELECT tests/test_plots.py)
add_pytest_test(NAME WavePacketDesigner.InputFile TOOL tools/wave_packet_designer SELECT tests/test_input_file.py)
add_pytest_test(NAME WavePacketDesigner.Cli TOOL tools/wave_packet_designer SELECT tests/test_cli.py)
add_pytest_test(NAME WavePacketDesigner.Structure TOOL tools/wave_packet_designer SELECT tests/test_structure.py)

# The plot tool (tools/plots): one ctest per test module, found by the glob so that a new module is never
# left out. Modules named test_runs_* read the fixture runs: they are PlotRuns.<name> (out of the quick tier)
# and require the fixture writer; the others are Plots.<name> and run without it.
file(GLOB plot_test_modules CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/tools/plots/tests/test_*.py)
foreach(module IN LISTS plot_test_modules)
  get_filename_component(stem ${module} NAME_WE)
  if(stem MATCHES "^test_runs_(.+)$")
    add_pytest_test(NAME PlotRuns.${CMAKE_MATCH_1} TOOL tools/plots SELECT tests/${stem}.py
                    ARGS --fixture-runs=${plot_fixtures_dir})
    set_tests_properties(PlotRuns.${CMAKE_MATCH_1} PROPERTIES FIXTURES_REQUIRED plot_fixtures)
  else()
    string(REGEX REPLACE "^test_" "" name ${stem})
    add_pytest_test(NAME Plots.${name} TOOL tools/plots SELECT tests/${stem}.py)
  endif()
endforeach()

# Every pytest test is selected by exactly one ctest (tools/check_python_tests.py), and the check finds a
# test without a ctest and a test with two (--plant-faults).
set(check_python_tests ${PROJECT_SOURCE_DIR}/tools/check_python_tests.py)
add_test(NAME Python.EveryTestRegistered
         COMMAND ${Python3_EXECUTABLE} -I -B ${check_python_tests} ${pytest_selections} --root ${PROJECT_SOURCE_DIR})
add_test(NAME Python.UnregisteredTestFails
         COMMAND ${Python3_EXECUTABLE} -I -B ${check_python_tests} ${pytest_selections} --root ${PROJECT_SOURCE_DIR}
                 --plant-faults)
set_tests_properties(Python.EveryTestRegistered Python.UnregisteredTestFails PROPERTIES
                     LABELS python
                     ENVIRONMENT "MPLBACKEND=Agg;MPLCONFIGDIR=${CMAKE_CURRENT_BINARY_DIR}/pytest/matplotlib")

# Every Fortran reference file listed in a MANIFEST.tsv exists with its recorded sha256, so a
# missing or edited reference fails instead of skipping its gate (tools/fortran_reference/README.md).
set(reference_manifest_check ${PROJECT_SOURCE_DIR}/tools/fortran_reference/ReferenceManifest.cmake)
add_test(NAME ReferenceManifest.FilesMatch
         COMMAND ${CMAKE_COMMAND} -DROOT=${PROJECT_SOURCE_DIR} -P ${reference_manifest_check})
add_test(NAME ReferenceManifest.ChangedOrMissingFileFails
         COMMAND ${CMAKE_COMMAND} -DROOT=${PROJECT_SOURCE_DIR}
                 -DWORK=${CMAKE_CURRENT_BINARY_DIR}/reference_manifest_negative -P ${reference_manifest_check})

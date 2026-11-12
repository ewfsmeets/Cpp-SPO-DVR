# Negative tests of tools/check_docs.py (ctest Docs.UndocumentedKeyFails, Docs.StaleAllowlistRowFails;
# README.md "Testing"): the checker must exit 1 and name the fault when
#   undocumented-key  src/input/json_input.cpp reads a key that docs/input_schema.md does not name
#                     (on a minimal copy of those two files in WORK, which is first shown to be clean);
#   stale-allowlist   the allow-list has one more row, which suppresses nothing (on the source tree).
# Usage: cmake -DPYTHON=<python3> -DCHECKER=<tools/check_docs.py> -DROOT=<source tree>
#              -DWORK=<scratch directory, recreated> -DCASE=undocumented-key|stale-allowlist
#              -P tools/check_docs_negative_test.cmake
# Environment: none read. Run time: under a second per case.

# Runs the checker with the given arguments; fails unless it exits with `expected_code` and its
# output matches `expected_regex`.
function(run_checker expected_code expected_regex)
  execute_process(COMMAND ${PYTHON} -I -B ${CHECKER} ${ARGN}
                  OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE code)
  if(NOT code EQUAL expected_code OR NOT out MATCHES "${expected_regex}")
    list(JOIN ARGN " " arguments)
    message(FATAL_ERROR "check_docs.py ${arguments}: exit ${code}, expected ${expected_code} and output "
                        "matching '${expected_regex}'; output:\n${out}${err}")
  endif()
endfunction()

foreach(variable PYTHON CHECKER ROOT WORK CASE)
  if(NOT DEFINED ${variable})
    message(FATAL_ERROR "error: -D${variable}=... is missing; see the usage at the top of this script")
  endif()
endforeach()
file(REMOVE_RECURSE ${WORK})

if(CASE STREQUAL "undocumented-key")
  file(COPY ${ROOT}/src/input/json_input.cpp DESTINATION ${WORK}/src/input)
  file(COPY ${ROOT}/docs/input_schema.md DESTINATION ${WORK}/docs)
  file(WRITE ${WORK}/allowlist.tsv "check\tfile\tmatch\treason\n")
  set(args --root ${WORK} --checks input-keys --allowlist ${WORK}/allowlist.tsv)
  run_checker(0 "check_docs: 0 finding" ${args})

  file(READ ${WORK}/src/input/json_input.cpp source)
  set(anchor "in.read(\"Title\", p.title, required);")
  string(FIND "${source}" "${anchor}" position)
  if(position EQUAL -1)
    message(FATAL_ERROR "error: '${anchor}' not found in json_input.cpp; update this test's anchor")
  endif()
  string(REPLACE "${anchor}" "${anchor}\n  in.read(\"UndocumentedTestKey\", p.title, required);" source "${source}")
  file(WRITE ${WORK}/src/input/json_input.cpp "${source}")
  run_checker(1 "src/input/json_input.cpp:[0-9]+: input-keys: key \"UndocumentedTestKey\" is read but not documented"
              ${args})
elseif(CASE STREQUAL "stale-allowlist")
  # The extra row sorts last (check and file names), so only its staleness is wrong.
  file(READ ${ROOT}/tools/check_docs_allowlist.tsv rows)
  file(WRITE ${WORK}/allowlist.tsv "${rows}output-files\tzz/no_such_file.md\tno finding has this text\tnegative test\n")
  run_checker(1 "allowlist: row \"output-files  zz/no_such_file.md  no finding has this text\" suppresses nothing"
              --root ${ROOT} --allowlist ${WORK}/allowlist.tsv --no-report)
else()
  message(FATAL_ERROR "error: unknown CASE '${CASE}'; use undocumented-key or stale-allowlist")
endif()

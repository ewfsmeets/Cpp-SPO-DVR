# Checks the Fortran reference files against their provenance manifests (ctest ReferenceManifest.*;
# tools/fortran_reference/README.md "MANIFEST.tsv"): every case directory of the reference trees has
# a MANIFEST.tsv, and every line of every tests/**/MANIFEST.tsv names an existing file inside its
# case directory whose sha256 is the recorded one. A missing or edited reference file, or a case
# without MANIFEST, fails here; the regression tests that read them would skip or compare with them.
#
# Usage: cmake -DROOT=<source tree> -P tools/fortran_reference/ReferenceManifest.cmake
#          the check (ctest ReferenceManifest.FilesMatch; regenerate_all.sh runs it at the end)
#        cmake -DROOT=<source tree> -DWORK=<scratch directory, recreated> -P ...
#          the negative test (ctest ReferenceManifest.ChangedOrMissingFileFails): on a copy of one
#          case, the check passes, then fails for one changed byte, one deleted file and a
#          deleted MANIFEST
# Environment: none read. Run time: under a second (it hashes every listed file: 0.1 s for the 172
# files listed on 2026-11-12).

cmake_minimum_required(VERSION 3.24)
set(manifest_header "path\tsha256\tscript\tcommit\tfortran_ref\tinput_sha\truntime_drift\tcpu\tdate")
set(manifest_columns 9)
# Every subdirectory of these trees is a reference case and needs a MANIFEST.tsv, except the ones
# that hold no Fortran output (refdata.py NOT_A_CASE).
set(reference_trees tests/reference tests/crplib/reference)
set(not_a_case fortran_unformatted)

# report_error(<text>): prints "error: <text>" unformatted (one line, so it can be searched) and
# counts it in the caller's `errors`.
macro(report_error text)
  message(NOTICE "error: ${text}")
  math(EXPR errors "${errors} + 1")
endmacro()

# check_manifest(<manifest> <errors_var> <files_var>): reports every bad line of <manifest> and adds
# the number of errors and of listed files to <errors_var> and <files_var>.
function(check_manifest manifest errors_var files_var)
  get_filename_component(case_dir "${manifest}" DIRECTORY)
  set(errors ${${errors_var}})
  set(files ${${files_var}})
  # One CMake list element per line; refdata.py never writes ';', '[' or ']', which would split a
  # line or join it with the next.
  file(STRINGS "${manifest}" lines ENCODING UTF-8)
  set(header_seen FALSE)
  foreach(line IN LISTS lines)
    if(line MATCHES "^#")
      continue()
    endif()
    string(FIND "${line}" "[" open)
    string(FIND "${line}" "]" close)
    if(line STREQUAL "" OR open GREATER -1 OR close GREATER -1)
      report_error("${manifest}: an empty line, or one with '[' or ']' (CMake joins it with the next); restore the file from git")
      continue()
    endif()
    if(NOT header_seen)
      if(NOT line STREQUAL manifest_header)
        break()
      endif()
      set(header_seen TRUE)
      continue()
    endif()
    string(REPLACE "\t" ";" fields "${line}")
    list(LENGTH fields count)
    list(GET fields 0 path)
    math(EXPR files "${files} + 1")
    if(NOT count EQUAL manifest_columns)
      report_error("${manifest}: the line of ${path} has ${count} fields, not ${manifest_columns}; restore the file from git")
      continue()
    endif()
    list(GET fields 1 recorded)
    if(path MATCHES "^/|(^|/)\\.\\.(/|$)")
      report_error("${manifest}: ${path} leaves the case directory; restore the file from git")
    elseif(NOT EXISTS "${case_dir}/${path}" OR IS_DIRECTORY "${case_dir}/${path}")
      report_error("${case_dir}/${path}: missing; restore it (git checkout) or regenerate it")
    else()
      file(SHA256 "${case_dir}/${path}" actual)
      if(NOT actual STREQUAL recorded)
        report_error("${case_dir}/${path}: sha256 differs from the MANIFEST; restore it (git checkout) or regenerate it")
      endif()
    endif()
  endforeach()
  if(NOT header_seen)
    report_error("${manifest}: no header row with the ${manifest_columns} MANIFEST columns; restore the file from git")
  endif()
  set(${errors_var} ${errors} PARENT_SCOPE)
  set(${files_var} ${files} PARENT_SCOPE)
endfunction()

# check_cases(<root> <errors_var> <cases_var>): reports every case directory without MANIFEST.tsv:
# after an interrupted regeneration, or for a case without a regeneration step, its gates would skip.
function(check_cases root errors_var cases_var)
  set(errors ${${errors_var}})
  set(cases 0)
  foreach(tree IN LISTS reference_trees)
    file(GLOB case_dirs LIST_DIRECTORIES true "${root}/${tree}/*")
    foreach(case_dir IN LISTS case_dirs)
      get_filename_component(name "${case_dir}" NAME)
      list(FIND not_a_case "${name}" excluded)
      if(NOT IS_DIRECTORY "${case_dir}" OR excluded GREATER -1)
        continue()
      endif()
      math(EXPR cases "${cases} + 1")
      if(NOT EXISTS "${case_dir}/MANIFEST.tsv")
        report_error("${case_dir}: no MANIFEST.tsv; restore it from git, or regenerate the case (tools/fortran_reference/README.md)")
      endif()
    endforeach()
  endforeach()
  set(${errors_var} ${errors} PARENT_SCOPE)
  set(${cases_var} ${cases} PARENT_SCOPE)
endfunction()

# check_tree(<root>): fails (FATAL_ERROR, exit 1) unless every case under <root> has a MANIFEST and
# every MANIFEST under <root>/tests matches.
function(check_tree root)
  set(errors 0)
  check_cases("${root}" errors cases)
  file(GLOB_RECURSE manifests LIST_DIRECTORIES false "${root}/tests/MANIFEST.tsv")
  list(LENGTH manifests manifest_count)
  if(cases EQUAL 0 AND manifest_count EQUAL 0)
    message(FATAL_ERROR "error: no reference case under ${root}/tests; pass the source tree as -DROOT")
  endif()
  set(files 0)
  foreach(manifest IN LISTS manifests)
    check_manifest("${manifest}" errors files)
  endforeach()
  if(errors GREATER 0)
    message(FATAL_ERROR "${errors} errors in ${cases} cases, ${manifest_count} MANIFESTs listing ${files} files (above)")
  endif()
  message(STATUS "ReferenceManifest: ${cases} cases, ${files} files in ${manifest_count} MANIFESTs match")
endfunction()

# expect_check(<root> <pass|fail> <regex>): runs the check on <root> in a new cmake; fails unless it
# passes or fails as expected with output matching <regex>.
function(expect_check root expected regex)
  execute_process(COMMAND "${CMAKE_COMMAND}" -DROOT=${root} -P "${CMAKE_CURRENT_LIST_FILE}"
                  OUTPUT_VARIABLE out ERROR_VARIABLE err RESULT_VARIABLE code)
  set(result fail)
  if(code EQUAL 0)
    set(result pass)
  endif()
  if(NOT result STREQUAL expected OR NOT "${out}${err}" MATCHES "${regex}")
    message(FATAL_ERROR "the check on ${root} should ${expected} with output matching '${regex}' "
                        "(exit ${code}); output:\n${out}${err}")
  endif()
  message(STATUS "as expected (${expected}): ${regex}")
endfunction()

if(NOT DEFINED ROOT)
  message(FATAL_ERROR "error: -DROOT=<source tree> is missing; see the usage at the top of this script")
endif()
if(NOT DEFINED WORK)
  check_tree("${ROOT}")
  return()
endif()

# The negative test, on a copy of a small case: the copy passes, then a changed byte, a deleted
# file and a deleted MANIFEST each fail, naming the file.
set(case tests/reference/cu111_persson_small)
set(changed fortran_setup_stdout.txt)
set(deleted debug_asyp.dat)
file(REMOVE_RECURSE "${WORK}")
file(COPY "${ROOT}/${case}/" DESTINATION "${WORK}/${case}")
expect_check("${WORK}" pass "1 cases, 3 files in 1 MANIFESTs match")

# The first byte replaced by another one: same size, so only the hash can tell.
file(READ "${WORK}/${case}/${changed}" first LIMIT 1)
file(READ "${WORK}/${case}/${changed}" rest OFFSET 1)
if(first STREQUAL "X")
  set(first "Y")
else()
  set(first "X")
endif()
file(WRITE "${WORK}/${case}/${changed}" "${first}${rest}")
expect_check("${WORK}" fail "${changed}: sha256 differs from the MANIFEST")

file(COPY "${ROOT}/${case}/${changed}" DESTINATION "${WORK}/${case}")
file(REMOVE "${WORK}/${case}/${deleted}")
expect_check("${WORK}" fail "${deleted}: missing")

file(COPY "${ROOT}/${case}/${deleted}" DESTINATION "${WORK}/${case}")
file(REMOVE "${WORK}/${case}/MANIFEST.tsv")
expect_check("${WORK}" fail "cu111_persson_small: no MANIFEST.tsv")

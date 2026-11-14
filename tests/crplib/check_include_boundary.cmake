# Checks the include boundaries of the CRP library and of the Fortran text reader (SRC = the src/
# directory):
# - a file of src/fortran_text includes only its own headers ("fortran_text/...") and the standard
#   library;
# - a file of the CRP library (src/crplib) includes only its own headers ("crplib/..."), those of
#   src/fortran_text and the standard library;
# - outside the library only the adapter src/potential/crp_potential.cpp includes crplib/ headers
#   (written "crplib/..." or <crplib/...>).
# A standard-library header is written <name>, the name made of letters, digits and underscores.
# C headers such as <math.h> are rejected on purpose: the house style is the C++ headers (<cmath>).
set(violations "")
file(GLOB_RECURSE files "${SRC}/*.hpp" "${SRC}/*.cpp" "${SRC}/*.cu" "${SRC}/*.cuh")
foreach(file IN LISTS files)
  file(STRINGS "${file}" includes REGEX "^[ \t]*#[ \t]*include[ \t]*[\"<]")
  string(FIND "${file}" "${SRC}/crplib/" in_library)
  string(FIND "${file}" "${SRC}/fortran_text/" in_reader)
  foreach(line IN LISTS includes)
    set(includes_library FALSE)
    set(includes_reader FALSE)
    set(includes_standard FALSE)
    if(line MATCHES "#[ \t]*include[ \t]*[\"<]crplib/")
      set(includes_library TRUE)
    elseif(line MATCHES "#[ \t]*include[ \t]*[\"<]fortran_text/")
      set(includes_reader TRUE)
    elseif(line MATCHES "#[ \t]*include[ \t]*<[A-Za-z0-9_]+>")
      set(includes_standard TRUE)
    endif()
    if(in_reader EQUAL 0)
      if(NOT includes_reader AND NOT includes_standard)
        string(APPEND violations
               "\n  ${file}: ${line} (src/fortran_text includes only its own headers and the standard library)")
      endif()
    elseif(in_library EQUAL 0)
      if(NOT includes_library AND NOT includes_reader AND NOT includes_standard)
        string(APPEND violations
               "\n  ${file}: ${line} (the library includes only crplib/ and fortran_text/ headers and the standard library)")
      endif()
    elseif(includes_library AND NOT file STREQUAL "${SRC}/potential/crp_potential.cpp")
      string(APPEND violations "\n  ${file}: ${line} (only potential/crp_potential.cpp uses the library)")
    endif()
  endforeach()
endforeach()
if(violations)
  message(FATAL_ERROR "Include boundary of the CRP library or the Fortran text reader violated:${violations}")
endif()
list(LENGTH files n)
message(STATUS "Include boundaries of the CRP library and the Fortran text reader OK (${n} files)")

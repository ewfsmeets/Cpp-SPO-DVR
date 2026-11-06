# Writes the bytes of a binary file as a C++ array into a header (used for the cuFFT callback's
# LTO-IR, src/cuda/fft_callbacks.cu). Run as a script:
#   cmake -DINPUT=<file> -DOUTPUT=<header> -DNAME=<array name> -P embed_binary.cmake
file(READ "${INPUT}" hex HEX)
string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," bytes "${hex}")
file(WRITE "${OUTPUT}"
     "// Generated from ${INPUT} by cmake/embed_binary.cmake. Do not edit.\n"
     "#pragma once\n"
     "alignas(8) inline constexpr unsigned char ${NAME}[] = {${bytes}};\n")

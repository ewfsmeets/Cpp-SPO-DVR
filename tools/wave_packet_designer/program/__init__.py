"""What spodvr computes or rejects, copied from the C++ so the designer needs no build.

Every module docstring has a `Shadows:` line naming the C++ files and functions it copies. Most modules
shadow one C++ file and have its name (validate.py <-> src/input/validate.cpp); some also copy small pieces
of other files (pi from constants.hpp in unit_conversion.py). To find the copy of a C++ file, search the
`Shadows:` lines: `grep -rn "spatial_grids.cpp" tools/wave_packet_designer`.

The test case `structure` checks that every cited file exists and contains the cited names, and that these
modules import only the standard library, numpy and each other (README, "Architecture").
"""

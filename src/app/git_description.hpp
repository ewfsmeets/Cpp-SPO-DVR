#pragma once

#include <string_view>

namespace spodvr {

// The git commit the spodvr executable was built from: `git describe --always --dirty` at build
// time, "unknown" outside a git checkout. Defined in the generated git_description.cpp
// (cmake/git_describe.cmake), which only the spodvr executable compiles; the program passes it to
// run_program as ProgramInfo (--version, run_info.json), never into the printed output.
[[nodiscard]] std::string_view git_description();

}  // namespace spodvr

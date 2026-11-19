#pragma once

#include <ostream>

#include "app/cli.hpp"
#include "core/program_info.hpp"

namespace spodvr {

// Runs the program for a parsed command line; returns the process exit code. Errors are
// reported on `err`. `program` (version and git commit) goes into run_info.json only; the tests
// pass the default ("unknown" commit).
int run_program(const CommandLine& command_line, std::ostream& out, std::ostream& err, bool color,
                const ProgramInfo& program = {});

}  // namespace spodvr

#pragma once

#include <ostream>

#include "app/cli.hpp"

namespace spodvr {

// Runs the program for a parsed command line; returns the process exit code. Errors are
// reported on `err`.
int run_program(const CommandLine& command_line, std::ostream& out, std::ostream& err, bool color);

}  // namespace spodvr

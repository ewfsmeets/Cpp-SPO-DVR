#pragma once

#include <optional>
#include <ostream>

#include "app/cli.hpp"
#include "core/program_info.hpp"
#include "propagation/propagation_backend.hpp"

namespace spodvr {

// Runs the program for a parsed command line; returns the process exit code. Errors are
// reported on `err`. `program` (version and git commit) goes into run_info.json only; the tests
// pass the default ("unknown" commit). `backend` replaces the device's backend factory: for tests only
// (tests/counting_backend.hpp counts the backend calls of the time loop).
int run_program(const CommandLine& command_line, std::ostream& out, std::ostream& err, bool color,
                const ProgramInfo& program = {}, const std::optional<BackendFactory>& backend = std::nullopt);

}  // namespace spodvr

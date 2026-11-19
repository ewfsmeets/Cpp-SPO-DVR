#pragma once

#include <string>

namespace spodvr {

// The program that wrote a run: its version and the git commit it was built from
// (`git describe --always --dirty`, compiled into the spodvr executable only; "unknown" in the
// tests and outside a git checkout). Recorded in run_info.json and printed by --version.
struct ProgramInfo {
  std::string version = SPODVR_VERSION;
  std::string git = "unknown";
};

}  // namespace spodvr

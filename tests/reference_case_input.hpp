#pragma once
// The C++ input of a reference case (tests/reference/<case>/input.json): its data-file paths are
// relative to the repository root there, as the Fortran reference scripts read them
// (tools/fortran_reference/reference_config.sh stage_inputs), while the tests run in the build tree.

#include <filesystem>
#include <memory>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

#include "input/config.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "potential/potential.hpp"
#include "potential/potential_factory.hpp"

namespace spodvr::testing {

// The repository root (SPODVR_REFERENCE_DIR is <root>/tests/reference).
inline std::filesystem::path repository_root() {
  return std::filesystem::path(SPODVR_REFERENCE_DIR).parent_path().parent_path();
}

// A path of input.json relative to the repository root, made absolute; empty (key absent) and
// absolute paths stay as they are.
inline void resolve_repository_path(std::string& path) {
  if (!path.empty()) path = (repository_root() / path).string();
}

// The data-file paths of a reference case's input made absolute: the CRP data set (CRPlibDataSet),
// the SCM parameter file (SCMParFilepath) and the slab position file (PosFilepath).
inline void resolve_repository_paths(InputParameters& input) {
  resolve_repository_path(input.crp_dataset);
  resolve_repository_path(input.scm.parameter_file);
  resolve_repository_path(input.scm.position_file);
}

// A reference case's PES and configuration, the PES built once (as the program does) from the resolved input.
struct CaseSetup {
  std::unique_ptr<Potential> pes;
  Config config;
};
inline CaseSetup case_setup(const nlohmann::json& input) {
  InputParameters parameters = parse_input(input);
  resolve_repository_paths(parameters);
  std::unique_ptr<Potential> pes = make_potential(parameters);
  Config config = build_config(parameters, *pes);
  return {std::move(pes), std::move(config)};
}

// The configuration of a reference case's JSON input, with its data-file paths resolved (case_setup's).
inline Config case_config(const nlohmann::json& input) { return case_setup(input).config; }

}  // namespace spodvr::testing

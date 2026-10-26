#pragma once

#include <filesystem>

#include <nlohmann/json.hpp>

#include "input/config.hpp"

namespace spodvr {

// Reads the JSON input. Keys are the Fortran DVR.inp field names (case-sensitive, see
// docs/input_schema.md). Unknown keys, missing required keys and wrongly typed values are
// InputErrors; keys of Fortran features that are not supported give a NotSupportedError.
// Fortran: InputNew.f90 ReadInput (plus the ERROR calls made while reading).
[[nodiscard]] InputParameters parse_input(const nlohmann::json& json);
[[nodiscard]] InputParameters read_input_file(const std::filesystem::path& path);

}  // namespace spodvr

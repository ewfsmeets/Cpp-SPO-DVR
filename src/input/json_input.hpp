#pragma once

#include <filesystem>
#include <istream>
#include <string>

#include <nlohmann/json.hpp>

#include "input/config.hpp"

namespace spodvr {

// Reads the JSON input. Keys are the Fortran DVR.inp field names (case-sensitive, see
// docs/input_schema.md). Unknown keys, missing required keys and wrongly typed values are
// InputErrors; keys of Fortran features that are not supported give a NotSupportedError.
// Fortran: InputNew.f90 ReadInput (plus the ERROR calls made while reading).
[[nodiscard]] InputParameters parse_input(const nlohmann::json& json);
[[nodiscard]] InputParameters read_input_file(const std::filesystem::path& path);

// The JSON object of an input text (`source` names it in messages). Comments are allowed; a key given
// twice in one object is an InputError (the Fortran reader stops on it too: inputfields.f90
// inp_Handle_Line, error code 4), and so is invalid JSON.
[[nodiscard]] nlohmann::json parse_input_text(std::istream& text, const std::string& source);

}  // namespace spodvr

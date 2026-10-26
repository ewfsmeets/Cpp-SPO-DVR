#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace spodvr {

struct CommandLine {
  std::filesystem::path input_file;
  std::optional<int> threads;
  std::filesystem::path output_dir = ".";
  bool no_color = false;
  bool help = false;
  bool version = false;
};

// Parses `spodvr <input.json> [--threads N] [--output-dir DIR] [--no-color] [--help] [--version]`.
// Throws std::invalid_argument with a message on errors.
[[nodiscard]] CommandLine parse_command_line(const std::vector<std::string>& args);
[[nodiscard]] std::string usage();

}  // namespace spodvr

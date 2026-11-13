#pragma once

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace spodvr {

// A tab-separated table of the result files (docs/output_files.md): '#' comment lines, one header
// line, then one row per record. Throws std::runtime_error if the file cannot be opened.
class TsvFile {
 public:
  TsvFile(const std::filesystem::path& path, const std::vector<std::string>& comments,
          const std::vector<std::string>& columns)
      : file_(path) {
    if (!file_) throw std::runtime_error("cannot open output file '" + path.string() + "'");
    for (const std::string& line : comments) file_ << "# " << line << '\n';
    write_row(columns);
  }

  // One line: the cells separated by tabs.
  void write_row(const std::vector<std::string>& cells) {
    for (std::size_t i = 0; i < cells.size(); ++i) file_ << (i == 0 ? "" : "\t") << cells[i];
    file_ << '\n';
  }

 private:
  std::ofstream file_;
};

}  // namespace spodvr

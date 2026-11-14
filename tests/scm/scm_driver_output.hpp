#pragma once
// Reading the output of the Fortran SCM reference driver (tools/fortran_reference/scm_driver/,
// tests/scm/reference/): lines "<tag> <index> <16 hex digits> <decimal>" (one double each, the
// index 1-based and consecutive per tag) and "<tag> <value>..." (text), ending with
// SCM_DRIVER_DONE. The tags and their order: scm_reference_driver.f90.

#include <bit>
#include <charconv>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace scm_test {

// The parsed output of one driver run.
struct ScmDriverOutput {
  std::map<std::string, std::vector<double>> values;  // per tag, its doubles in index order
  std::map<std::string, std::string> text;            // per tag, the rest of its text line
  std::vector<int> atom_layers;                       // the "atom_layer <i> <layer>" lines, in order
  bool done = false;                                  // the output ended with SCM_DRIVER_DONE

  // The one double of a tag such as "cutoff".
  double value(const std::string& tag) const { return values.at(tag).at(0); }
};

// True for a token of 16 hexadecimal digits (a double's bit pattern).
inline bool is_hex_word(const std::string& token) {
  return token.size() == 16 && token.find_first_not_of("0123456789ABCDEFabcdef") == std::string::npos;
}

// The double whose bit pattern a 16-digit hex token gives.
inline double from_hex_word(const std::string& token) {
  std::uint64_t bits = 0;
  std::from_chars(token.data(), token.data() + token.size(), bits, 16);
  return std::bit_cast<double>(bits);
}

// Reads a driver output. Throws std::runtime_error, naming the file and line, for a double line
// whose index is not the next one of its tag, or for a line after SCM_DRIVER_DONE.
inline ScmDriverOutput read_scm_driver_output(const std::filesystem::path& file) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("cannot open " + file.string());
  ScmDriverOutput output;
  std::string line;
  for (int number = 1; std::getline(in, line); ++number) {
    const std::string where = file.string() + ":" + std::to_string(number);
    if (output.done) throw std::runtime_error(where + ": a line after SCM_DRIVER_DONE");
    if (line == "SCM_DRIVER_DONE") {
      output.done = true;
      continue;
    }
    std::istringstream fields(line);
    std::string tag, index, word;
    fields >> tag >> index >> word;
    if (is_hex_word(word)) {
      std::vector<double>& values = output.values[tag];
      if (std::stoi(index) != std::ssize(values) + 1)
        throw std::runtime_error(where + ": index " + index + " of " + tag + " out of order");
      values.push_back(from_hex_word(word));
    } else if (tag == "atom_layer") {
      if (std::stoi(index) != std::ssize(output.atom_layers) + 1)
        throw std::runtime_error(where + ": atom " + index + " out of order");
      output.atom_layers.push_back(std::stoi(word));
    } else {
      output.text[tag] = line.size() > tag.size() ? line.substr(tag.size() + 1) : "";
    }
  }
  return output;
}

}  // namespace scm_test

#pragma once
// Reading the input and output files of the Fortran CRP reference driver
// (tools/fortran_reference/crp_driver/, tests/crplib/reference/). Inputs are whitespace-separated
// tokens: integers and 16-digit hexadecimal IEEE doubles. Output lines are
// "<label> <index> <hex> <value>" (results), or text (e.g. "Q outside", "CRP_DRIVER_DONE").

#include <gtest/gtest.h>

#include <bit>
#include <charconv>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace crplib_test {

inline std::filesystem::path reference_file(const std::string& name) {
  return std::filesystem::path(CRPLIB_REFERENCE_DIR) / name;
}

// The data set of a reference case: <case>/dataset for the synthetic cases; the case cu111_srp48
// uses the real data set of the repository (potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v).
inline std::filesystem::path reference_data_set(const std::string& case_name) {
  if (case_name == "cu111_srp48") return CRPLIB_REAL_DATA_SET;
  return reference_file(case_name + "/dataset");
}

#define REQUIRE_CRP_REFERENCE(path)   \
  if (!std::filesystem::exists(path)) \
  GTEST_SKIP() << "reference file " << (path) << " missing; run tools/fortran_reference/run_crp_reference.sh"

inline std::uint64_t parse_hex(const std::string& text) {
  std::uint64_t bits = 0;
  std::from_chars(text.data(), text.data() + text.size(), bits, 16);
  return bits;
}
inline double hex_to_double(const std::string& text) { return std::bit_cast<double>(parse_hex(text)); }
inline std::uint64_t bits_of(double x) { return std::bit_cast<std::uint64_t>(x); }

// The tokens of a driver input file, read in order.
class DriverInput {
 public:
  explicit DriverInput(const std::filesystem::path& file) {
    std::ifstream in(file);
    std::string token;
    while (in >> token) tokens_.push_back(token);
  }
  int integer() { return std::stoi(tokens_.at(next_++)); }
  double real() { return hex_to_double(tokens_.at(next_++)); }
  std::vector<double> reals(int n) {
    std::vector<double> values(n);
    for (double& v : values) v = real();
    return values;
  }

 private:
  std::vector<std::string> tokens_;
  std::size_t next_ = 0;
};

// One line of driver output: a result (label, 1-based index, bits) or text.
struct DriverRecord {
  std::string label;
  int index = 0;
  bool has_value = false;
  std::uint64_t bits = 0;
  std::string text;
};

// The output lines in order; fails the test if the run did not end with CRP_DRIVER_DONE.
inline std::vector<DriverRecord> read_driver_output(const std::filesystem::path& file) {
  std::ifstream in(file);
  std::vector<DriverRecord> records;
  std::string line;
  bool done = false;
  while (std::getline(in, line)) {
    if (line == "CRP_DRIVER_DONE") {
      done = true;
      continue;
    }
    std::istringstream fields(line);
    DriverRecord r;
    std::string index, hex;
    fields >> r.label >> index >> hex;
    r.text = line;
    if (hex.size() == 16 && !index.empty() && index.find_first_not_of("0123456789") == std::string::npos) {
      r.index = std::stoi(index);
      r.bits = parse_hex(hex);
      r.has_value = true;
    }
    records.push_back(r);
  }
  EXPECT_TRUE(done) << file << " does not end with CRP_DRIVER_DONE";
  return records;
}

// Sequential access to the result records of a driver output.
class DriverOutput {
 public:
  explicit DriverOutput(const std::filesystem::path& file) : records_(read_driver_output(file)) {}
  // The next record; it must carry `label`.
  const DriverRecord& next(const std::string& label) {
    const DriverRecord& r = records_.at(next_++);
    EXPECT_EQ(r.label, label) << r.text;
    return r;
  }
  // The values of the next n records with `label` (an array printed element by element).
  std::vector<std::uint64_t> next_bits(const std::string& label, int n) {
    std::vector<std::uint64_t> bits(n);
    for (auto& b : bits) b = next(label).bits;
    return bits;
  }
  [[nodiscard]] bool at_end() const { return next_ >= records_.size(); }

 private:
  std::vector<DriverRecord> records_;
  std::size_t next_ = 0;
};

// Expects bitwise equality of a C++ result with the Fortran bits (prints both values on failure).
#define EXPECT_BITS(cpp_value, fortran_bits)                                                                        \
  EXPECT_EQ(crplib_test::bits_of(cpp_value), (fortran_bits))                                                        \
      << #cpp_value << " = " << std::hexfloat << (cpp_value) << ", Fortran " << std::bit_cast<double>(fortran_bits) \
      << std::defaultfloat

}  // namespace crplib_test

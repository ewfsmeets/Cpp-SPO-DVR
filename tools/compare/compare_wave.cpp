// Compares two wave packets written by spodvr_dump_wave: bitwise identical, and the largest
// difference relative to the largest amplitude; a NaN or inf in either file always fails, naming
// the first such element. Docs: docs/performance.md, section "Benchmarking checklist".
// Usage: spodvr_compare_wave <a> <b> [tolerance]   (tolerance: a finite number >= 0, default 0)
// Exit 0 within the tolerance, 1 different or a file refused ("error: <what>; <fix>"), 2 usage.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <format>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include "io/binary_file.hpp"

namespace {

// Printed on stderr for every bad command line.
constexpr std::string_view usage_line =
    "usage: spodvr_compare_wave <a> <b> [tolerance]   (tolerance: a finite number >= 0, default 0)\n";

// The tolerance argument: a finite number >= 0, written in full; nullopt otherwise.
std::optional<double> parse_tolerance(const std::string& text) {
  std::size_t used = 0;
  double value = 0.0;
  try {
    value = std::stod(text, &used);
  } catch (const std::exception&) {
    return std::nullopt;
  }
  if (used != text.size() || !std::isfinite(value) || value < 0.0) return std::nullopt;
  return value;
}

// The doubles (re, im pairs) of a file written by spodvr_dump_wave; throws CheckpointError for a
// file that cannot be read, is empty or does not hold a whole number of complex values.
std::vector<double> read_values(const std::filesystem::path& path) {
  spodvr::BinaryReader in(path);
  std::error_code error;
  const std::uintmax_t bytes = std::filesystem::file_size(path, error);
  if (error || bytes == 0 || bytes % sizeof(std::complex<double>) != 0)
    throw spodvr::CheckpointError("'" + path.string() + "' is empty, unreadable or not a list of complex values");
  std::vector<double> values(bytes / sizeof(double));
  in.read_array<double>(values);
  return values;
}

// The first element index at which a or b is not finite, or a.size() if there is none.
std::size_t first_non_finite(const std::vector<double>& a, const std::vector<double>& b) {
  std::size_t i = 0;
  while (i < a.size() && std::isfinite(a[i]) && std::isfinite(b[i])) ++i;
  return i;
}

// Prints the comparison of two wave packets of equal size (named name_a, name_b); returns the
// exit status.
int compare(const std::vector<double>& a, const std::vector<double>& b, const std::string& name_a,
            const std::string& name_b, double tolerance) {
  // A NaN or inf fails whatever the tolerance (the std::max reductions below would drop a NaN).
  if (const std::size_t i = first_non_finite(a, b); i < a.size()) {
    const bool in_a = !std::isfinite(a[i]);
    std::cout << std::format("DIFFERENT: non-finite value {:g} in {} at element {} (complex value {}, {} part)\n",
                             in_a ? a[i] : b[i], in_a ? name_a : name_b, i, i / 2, i % 2 == 0 ? "real" : "imaginary");
    return 1;
  }
  bool identical = true;
  double max_difference = 0.0;
  double max_amplitude = 0.0;
  for (std::size_t i = 0; i + 1 < a.size(); i += 2) {
    const double dre = a[i] - b[i];
    const double dim = a[i + 1] - b[i + 1];
    max_difference = std::max(max_difference, std::hypot(dre, dim));
    max_amplitude = std::max(max_amplitude, std::hypot(a[i], a[i + 1]));
    if (a[i] != b[i] || a[i + 1] != b[i + 1]) identical = false;
  }
  const double relative = max_amplitude > 0.0 ? max_difference / max_amplitude : max_difference;
  std::cout << std::format("{}: max |a - b| / max |a| = {:.3e}\n", identical ? "bitwise identical" : "different",
                           relative);
  return relative <= tolerance ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  const std::vector<std::string> args(argv + 1, argv + argc);
  const std::optional<double> tolerance = args.size() == 3 ? parse_tolerance(args[2]) : std::optional(0.0);
  if (args.size() < 2 || args.size() > 3 || !tolerance) {
    std::cerr << usage_line;
    return 2;
  }
  std::vector<double> a;
  std::vector<double> b;
  try {
    a = read_values(args[0]);
    b = read_values(args[1]);
  } catch (const spodvr::CheckpointError& e) {
    std::cerr << "error: " << e.what() << "; give two files written by spodvr_dump_wave\n";
    return 1;
  }
  if (a.size() != b.size()) {
    std::cout << std::format("DIFFERENT: sizes {} and {}\n", a.size(), b.size());
    return 1;
  }
  return compare(a, b, args[0], args[1], *tolerance);
}

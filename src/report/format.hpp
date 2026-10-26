#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <string_view>
#include <vector>

namespace spodvr {

// Output styling: ANSI colours on or off.
struct Style {
  bool color = false;

  [[nodiscard]] std::string heading(std::string_view text) const;
  [[nodiscard]] std::string emphasis(std::string_view text) const;
  [[nodiscard]] std::string warning(std::string_view text) const;
};

// Full-precision scientific notation (%.16e), used for every regression-relevant number.
[[nodiscard]] std::string full_precision(double x);
// Fixed or scientific formatting with the given number of significant/decimal digits.
[[nodiscard]] std::string fixed(double x, int decimals);
[[nodiscard]] std::string scientific(double x, int decimals);
// Memory size in binary units: "512 B", "3.2 MiB", "4.0 GiB".
[[nodiscard]] std::string format_bytes(std::size_t bytes);

// Section header with a rule above and below.
void print_section(std::ostream& out, const Style& style, std::string_view title);

// Aligned "label ..... value [unit]" lines.
class KeyValueList {
 public:
  void add(std::string label, std::string value, std::string unit = "");
  void print(std::ostream& out, const Style& style) const;

 private:
  struct Row {
    std::string label;
    std::string value;
    std::string unit;
  };
  std::vector<Row> rows_;
};

// Simple right-aligned text table with column headers (units belong in the headers).
class Table {
 public:
  explicit Table(std::vector<std::string> headers) : headers_(std::move(headers)) {}
  void add_row(std::vector<std::string> cells) { rows_.push_back(std::move(cells)); }
  void print(std::ostream& out, const Style& style) const;

 private:
  std::vector<std::string> headers_;
  std::vector<std::vector<std::string>> rows_;
};

}  // namespace spodvr

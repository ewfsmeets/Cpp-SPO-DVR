// Text formatting: number formats, ANSI styles, aligned lists and tables.
#include "report/format.hpp"

#include <algorithm>
#include <array>
#include <cstdio>

namespace spodvr {

namespace {

std::string wrap(bool color, std::string_view code, std::string_view text) {
  if (!color) return std::string(text);
  return std::string("\033[") + std::string(code) + "m" + std::string(text) + "\033[0m";
}

std::string format(const char* fmt, int decimals, double x) {
  char buffer[64];
  std::snprintf(buffer, sizeof buffer, fmt, decimals, x);
  return buffer;
}

}  // namespace

std::string Style::heading(std::string_view text) const { return wrap(color, "1;36", text); }
std::string Style::emphasis(std::string_view text) const { return wrap(color, "1", text); }
std::string Style::warning(std::string_view text) const { return wrap(color, "1;33", text); }

std::string full_precision(double x) { return format("%.*e", 16, x); }
std::string fixed(double x, int decimals) { return format("%.*f", decimals, x); }

std::string format_bytes(std::size_t bytes) {
  if (bytes < 1024) return std::to_string(bytes) + " B";
  const std::array<std::string, 4> units = {"KiB", "MiB", "GiB", "TiB"};
  double size = static_cast<double>(bytes) / 1024.0;
  int unit = 0;
  while (size >= 1024.0 && unit < 3) {
    size = size / 1024.0;
    unit = unit + 1;
  }
  return fixed(size, 1) + " " + units[unit];
}
std::string scientific(double x, int decimals) { return format("%.*e", decimals, x); }

void print_section(std::ostream& out, const Style& style, std::string_view title) {
  const std::string rule(80, '=');
  out << '\n'
      << style.heading(rule) << '\n'
      << style.heading(std::string("  ") + std::string(title)) << '\n'
      << style.heading(rule) << "\n\n";
}

void KeyValueList::add(std::string label, std::string value, std::string unit) {
  rows_.push_back({std::move(label), std::move(value), std::move(unit)});
}

void KeyValueList::print(std::ostream& out, const Style& /*style*/) const {
  std::size_t width = 0;
  for (const auto& row : rows_) width = std::max(width, row.label.size());
  for (const auto& row : rows_) {
    out << "  " << row.label << ' ' << std::string(width - row.label.size() + 2, '.') << ' ' << row.value;
    if (!row.unit.empty()) out << ' ' << row.unit;
    out << '\n';
  }
}

void Table::print(std::ostream& out, const Style& style) const {
  std::vector<std::size_t> widths(headers_.size());
  for (std::size_t c = 0; c < headers_.size(); ++c) widths[c] = headers_[c].size();
  for (const auto& row : rows_)
    for (std::size_t c = 0; c < row.size() && c < widths.size(); ++c) widths[c] = std::max(widths[c], row[c].size());

  auto line = [&](const std::vector<std::string>& cells) {
    std::string text;
    for (std::size_t c = 0; c < widths.size(); ++c) {
      const std::string cell = c < cells.size() ? cells[c] : "";
      text += "  " + std::string(widths[c] - cell.size(), ' ') + cell;
    }
    return text;
  };
  out << style.emphasis(line(headers_)) << '\n';
  std::size_t total = 0;
  for (const auto w : widths) total += w + 2;
  out << std::string(total, '-') << '\n';
  for (const auto& row : rows_) out << line(row) << '\n';
}

}  // namespace spodvr

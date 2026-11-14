// List-directed reading of text files as gfortran's READ(unit,*) does (the subset that the CRP
// data-set files use). Fortran: the READ statements of CRP_Wrapper.f90 and CRP_C*v_*Surface.f90.
#include "fortran_text/list_directed_reader.hpp"

#include <cctype>
#include <charconv>
#include <fstream>
#include <system_error>

namespace fortran_text {

namespace {

bool is_blank(char c) { return c == ' ' || c == '\t' || c == '\r'; }
bool is_digit(char c) { return std::isdigit(static_cast<unsigned char>(c)) != 0; }

// The number of digits starting at text[pos].
std::size_t count_digits(std::string_view text, std::size_t pos) {
  std::size_t n = 0;
  while (pos + n < text.size() && is_digit(text[pos + n])) ++n;
  return n;
}

// std::from_chars over a whole string; false unless every character is consumed. (The one place
// where the reader hands character pointers to the standard library.)
template <typename T>
bool from_chars_whole(const std::string& text, T& value) {
  const char* first = text.data();
  const char* last = first + text.size();
  const auto [end, error] = std::from_chars(first, last, value);
  return error == std::errc() && end == last;
}

}  // namespace

// Fortran real constant: [sign] (digits [. [digits]] | . digits) [(E|e|D|d) [sign] digits].
// gfortran reads it correctly rounded, as std::from_chars does.
std::optional<double> parse_fortran_real(std::string_view text) {
  std::size_t pos = 0;
  std::string normalised;
  if (pos < text.size() && (text[pos] == '+' || text[pos] == '-')) {
    if (text[pos] == '-') normalised += '-';
    ++pos;
  }
  const std::size_t integer_digits = count_digits(text, pos);
  normalised += text.substr(pos, integer_digits);
  pos += integer_digits;
  std::size_t fraction_digits = 0;
  if (pos < text.size() && text[pos] == '.') {
    ++pos;
    fraction_digits = count_digits(text, pos);
    normalised += '.';
    normalised += text.substr(pos, fraction_digits);
    pos += fraction_digits;
  }
  if (integer_digits + fraction_digits == 0) return std::nullopt;
  if (pos < text.size()) {
    const char letter = text[pos];
    if (letter != 'E' && letter != 'e' && letter != 'D' && letter != 'd') return std::nullopt;
    ++pos;
    normalised += 'e';
    if (pos < text.size() && (text[pos] == '+' || text[pos] == '-')) {
      if (text[pos] == '-') normalised += '-';
      ++pos;
    }
    const std::size_t exponent_digits = count_digits(text, pos);
    if (exponent_digits == 0) return std::nullopt;
    normalised += text.substr(pos, exponent_digits);
    pos += exponent_digits;
  }
  if (pos != text.size()) return std::nullopt;
  double value = 0.0;
  if (!from_chars_whole(normalised, value)) return std::nullopt;
  return value;
}

void ListStatement::fail(int i, std::string_view expected) const {
  const Item& item = items_[i];
  throw ReadError(file_ + ":" + std::to_string(item.line) + ": '" + item.text + "' is not " + std::string(expected));
}

double ListStatement::real(int i) const {
  const std::optional<double> value = parse_fortran_real(items_[i].text);
  if (!value) fail(i, "a real number");
  return *value;
}

int ListStatement::integer(int i) const {
  std::string text = items_[i].text;
  if (!text.empty() && text[0] == '+') text.erase(0, 1);
  int value = 0;
  if (text.empty() || !from_chars_whole(text, value)) fail(i, "an integer");
  return value;
}

// Fortran list-directed LOGICAL: an optional period, then T or F (any case), then anything.
bool ListStatement::logical(int i) const {
  std::string_view text = items_[i].text;
  if (!text.empty() && text[0] == '.') text.remove_prefix(1);
  if (!text.empty() && (text[0] == 'T' || text[0] == 't')) return true;
  if (!text.empty() && (text[0] == 'F' || text[0] == 'f')) return false;
  fail(i, "a logical value (T or F)");
}

const std::string& ListStatement::text(int i) const { return items_[i].text; }

ListDirectedReader::ListDirectedReader(const std::filesystem::path& file) : file_(file.string()) {
  std::ifstream in(file);
  if (!in) throw ReadError(file_ + ": cannot open the file");
  std::string line;
  while (std::getline(in, line)) lines_.push_back(line);
}

void ListDirectedReader::fail(std::size_t line_index, const std::string& message) const {
  throw ReadError(file_ + ":" + std::to_string(line_index + 1) + ": " + message);
}

// Appends the values of record `line_index` until `n_values` are collected; the rest of the
// record is not looked at (Fortran skips it). `expecting_value` carries across records: a comma
// while a value is expected would be a null value.
void ListDirectedReader::split_record(std::size_t line_index, std::size_t n_values,
                                      std::vector<ListStatement::Item>& items, bool& expecting_value) const {
  const std::string& line = lines_[line_index];
  const int line_number = static_cast<int>(line_index) + 1;
  std::size_t pos = 0;
  while (pos < line.size() && items.size() < n_values) {
    const char c = line[pos];
    if (is_blank(c)) {
      ++pos;
    } else if (c == ',') {
      if (expecting_value) fail(line_index, "null values (two commas in a row) are not supported");
      expecting_value = true;
      ++pos;
    } else if (c == '/') {
      fail(line_index, "'/' (end of the input values) is not supported");
    } else if (c == '\'' || c == '"') {
      std::string text;
      bool closed = false;
      for (++pos; pos < line.size(); ++pos) {
        if (line[pos] != c) {
          text += line[pos];
        } else if (pos + 1 < line.size() && line[pos + 1] == c) {
          text += c;  // doubled quote
          ++pos;
        } else {
          closed = true;
          ++pos;
          break;
        }
      }
      if (!closed) fail(line_index, "unterminated character string");
      items.push_back({text, line_number});
      expecting_value = false;
    } else {
      const std::size_t start = pos;
      while (pos < line.size() && !is_blank(line[pos]) && line[pos] != ',' && line[pos] != '/') ++pos;
      std::string text = line.substr(start, pos - start);
      if (text.find('*') != std::string::npos) fail(line_index, "repeat counts (r*c) are not supported: " + text);
      items.push_back({text, line_number});
      expecting_value = false;
    }
  }
}

std::optional<ListStatement> ListDirectedReader::try_read(int n_values) {
  if (next_line_ >= lines_.size()) return std::nullopt;
  std::vector<ListStatement::Item> items;
  bool expecting_value = true;
  if (n_values == 0) ++next_line_;  // an empty statement still skips a record
  while (static_cast<int>(items.size()) < n_values) {
    if (next_line_ >= lines_.size()) return std::nullopt;
    split_record(next_line_, n_values, items, expecting_value);
    ++next_line_;
  }
  return ListStatement(file_, std::move(items));
}

ListStatement ListDirectedReader::read(int n_values) {
  std::optional<ListStatement> statement = try_read(n_values);
  if (!statement) {
    throw ReadError(file_ + ": unexpected end of the file (reading " + std::to_string(n_values) +
                    " values after line " + std::to_string(next_line_) + ")");
  }
  return std::move(*statement);
}

std::optional<std::string> ListDirectedReader::read_line() {
  if (next_line_ >= lines_.size()) return std::nullopt;
  std::string line = lines_[next_line_++];
  if (!line.empty() && line.back() == '\r') line.pop_back();
  return line;
}

}  // namespace fortran_text

#pragma once
// Reading text files the way gfortran's list-directed READ(unit,*) does (the subset that the CRP
// data-set files need). Fortran: the READ statements of CRP_Wrapper.f90 and CRP_C*v_*Surface.f90.

#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace fortran_text {

// A file cannot be opened, or its text is not what a READ statement expects. The message names
// the file (and the line, for reading errors). Users translate it into their own error type.
class ReadError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// The values of one READ statement, converted on access. Errors name the file and line.
class ListStatement {
 public:
  struct Item {
    std::string text;  // the value as written (quotes removed from character values)
    int line;          // 1-based line number in the file
  };
  ListStatement(std::string file, std::vector<Item> items) : file_(std::move(file)), items_(std::move(items)) {}

  [[nodiscard]] int size() const { return static_cast<int>(items_.size()); }
  // Item i (0-based) as a REAL (D or E exponent), INTEGER, LOGICAL (T/F/.TRUE./.FALSE., any
  // case) or CHARACTER value. Throws ReadError if the item is not of that type.
  [[nodiscard]] double real(int i) const;
  [[nodiscard]] int integer(int i) const;
  [[nodiscard]] bool logical(int i) const;
  [[nodiscard]] const std::string& text(int i) const;

 private:
  [[noreturn]] void fail(int i, std::string_view expected) const;

  std::string file_;
  std::vector<Item> items_;
};

// A text file read statement by statement, as with READ(unit,*). Only the subset of list-directed
// input that the CRP files need is accepted:
// - values are separated by blanks, tabs or a comma; one statement reads from the next record
//   (line) on, continues on the following records until it has all its values, and skips the
//   rest of its last record; a statement with no values skips one record;
// - character values may be quoted ('...' or "...", doubled quotes inside); unquoted ones end at
//   a blank or comma.
// Rejected with a ReadError (instead of the Fortran's special meaning): repeat counts (r*c),
// null values (two commas in a row), the slash that ends the input of a statement.
class ListDirectedReader {
 public:
  // Reads the whole file; throws ReadError if it cannot be read.
  explicit ListDirectedReader(const std::filesystem::path& file);

  // READ(unit,*) with n values. Throws ReadError at the end of the file.
  ListStatement read(int n_values);
  // READ(unit,*,END=...) with n values: nullopt at the end of the file.
  std::optional<ListStatement> try_read(int n_values);
  // READ(unit,'(A)',END=...): the next record as it is; nullopt at the end of the file.
  std::optional<std::string> read_line();

  [[nodiscard]] const std::string& file_name() const { return file_; }

 private:
  void split_record(std::size_t line_index, std::size_t n_values, std::vector<ListStatement::Item>& items,
                    bool& expecting_value) const;
  [[noreturn]] void fail(std::size_t line_index, const std::string& message) const;

  std::string file_;
  std::vector<std::string> lines_;
  std::size_t next_line_ = 0;
};

// A Fortran real constant (optional sign, digits with an optional decimal point, optional E or D
// exponent), correctly rounded like gfortran's reading; nullopt if the text is not one.
std::optional<double> parse_fortran_real(std::string_view text);

}  // namespace fortran_text

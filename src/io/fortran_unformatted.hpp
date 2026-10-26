#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "io/binary_file.hpp"

namespace spodvr {

// gfortran sequential unformatted files (FORM='unformatted'): every WRITE is one record,
// framed by 4-byte length markers. Records longer than the maximum subrecord length
// (2^31 - 9 bytes) are split into subrecords; the head marker of a subrecord is negative if
// another subrecord follows, the tail marker is negative if a subrecord precedes it.
// INTEGER and LOGICAL are 4 bytes (.TRUE. = 1); REAL and COMPLEX are 8 and 16 bytes in the
// reference build (-fdefault-real-8).

// The bytes of one record, assembled value by value.
class FortranRecord {
 public:
  template <typename T>
  FortranRecord& add(const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    return add_bytes(std::as_bytes(std::span(&value, 1)));
  }
  template <typename T>
  FortranRecord& add_array(std::span<const T> values) {
    static_assert(std::is_trivially_copyable_v<T>);
    return add_bytes(std::as_bytes(values));
  }
  [[nodiscard]] std::span<const std::byte> bytes() const { return bytes_; }

 private:
  FortranRecord& add_bytes(std::span<const std::byte> b) {
    bytes_.insert(bytes_.end(), b.begin(), b.end());
    return *this;
  }
  std::vector<std::byte> bytes_;
};

class FortranRecordWriter {
 public:
  static constexpr std::int64_t default_max_subrecord = 2147483639;  // gfortran default

  explicit FortranRecordWriter(const std::filesystem::path& path, std::int64_t max_subrecord = default_max_subrecord)
      : file_(path), max_subrecord_(max_subrecord) {}

  void write(const FortranRecord& record) { write_bytes(record.bytes()); }
  // A record holding one array (avoids a copy for large arrays).
  template <typename T>
  void write_array(std::span<const T> values) {
    write_bytes(std::as_bytes(values));
  }
  void flush() { file_.flush(); }
  void close() { file_.close(); }

 private:
  void write_bytes(std::span<const std::byte> data);

  BinaryWriter file_;
  std::int64_t max_subrecord_;
};

// Values of one record, read in order; reading past the end or leaving data unread is an error.
class FortranRecordParser {
 public:
  FortranRecordParser(std::vector<std::byte> bytes, std::string description)
      : bytes_(std::move(bytes)), description_(std::move(description)) {}

  template <typename T>
  [[nodiscard]] T read() {
    T value;
    read_bytes(std::as_writable_bytes(std::span(&value, 1)));
    return value;
  }
  template <typename T>
  void read_array(std::span<T> values) {
    read_bytes(std::as_writable_bytes(values));
  }
  [[nodiscard]] std::size_t size() const { return bytes_.size(); }
  void expect_end() const;

 private:
  void read_bytes(std::span<std::byte> out);

  std::vector<std::byte> bytes_;
  std::string description_;
  std::size_t position_ = 0;
};

class FortranRecordReader {
 public:
  explicit FortranRecordReader(const std::filesystem::path& path) : file_(path) {}

  [[nodiscard]] bool at_end() { return file_.at_end(); }
  // The next record; `what` names it in error messages.
  [[nodiscard]] FortranRecordParser read(const std::string& what);
  [[nodiscard]] const std::filesystem::path& path() const { return file_.path(); }

 private:
  BinaryReader file_;
};

}  // namespace spodvr

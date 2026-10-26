#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace spodvr {

// Error in a checkpoint or BK coefficient file (missing, corrupt, or not matching the input).
class CheckpointError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Raw binary output of plain values and arrays in host byte order (little-endian on x86-64),
// for the native files of this program.
class BinaryWriter {
 public:
  BinaryWriter(const std::filesystem::path& path, std::ios::openmode mode = std::ios::trunc)
      : path_(path), file_(path, std::ios::binary | std::ios::out | mode) {
    if (!file_) throw CheckpointError("cannot open '" + path.string() + "' for writing");
  }

  template <typename T>
  void write(const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    write_bytes(std::as_bytes(std::span(&value, 1)));
  }
  template <typename T>
  void write_array(std::span<const T> values) {
    static_assert(std::is_trivially_copyable_v<T>);
    write_bytes(std::as_bytes(values));
  }
  void write_string(const std::string& s) {
    write(static_cast<std::uint64_t>(s.size()));
    write_bytes(std::as_bytes(std::span(s)));
  }
  void flush() {
    file_.flush();
    check();
  }
  void close() {
    file_.close();
    check();
  }

 private:
  // The only place where bytes meet the iostream char* interface.
  void write_bytes(std::span<const std::byte> bytes) {
    file_.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    check();
  }
  void check() const {
    if (!file_) throw CheckpointError("error writing '" + path_.string() + "'");
  }

  std::filesystem::path path_;
  std::ofstream file_;
};

class BinaryReader {
 public:
  explicit BinaryReader(const std::filesystem::path& path) : path_(path), file_(path, std::ios::binary) {
    if (!file_) throw CheckpointError("cannot open '" + path.string() + "'");
  }

  template <typename T>
  [[nodiscard]] T read() {
    static_assert(std::is_trivially_copyable_v<T>);
    T value;
    read_bytes(std::as_writable_bytes(std::span(&value, 1)));
    return value;
  }
  template <typename T>
  void read_array(std::span<T> values) {
    static_assert(std::is_trivially_copyable_v<T>);
    read_bytes(std::as_writable_bytes(values));
  }
  [[nodiscard]] std::string read_string() {
    const auto size = read<std::uint64_t>();
    if (size > (1u << 20)) throw CheckpointError("corrupt file '" + path_.string() + "' (string length)");
    std::string s(size, '\0');
    read_bytes(std::as_writable_bytes(std::span(s)));
    return s;
  }
  // True if no data is left.
  [[nodiscard]] bool at_end() { return file_.peek() == std::ifstream::traits_type::eof(); }
  [[nodiscard]] const std::filesystem::path& path() const { return path_; }

 private:
  void read_bytes(std::span<std::byte> bytes) {
    file_.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (!file_) throw CheckpointError("file '" + path_.string() + "' is truncated or unreadable");
  }

  std::filesystem::path path_;
  std::ifstream file_;
};

}  // namespace spodvr

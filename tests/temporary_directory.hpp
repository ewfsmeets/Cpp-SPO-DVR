#pragma once
// A directory for files written by a test, removed again when the test ends.

#include <atomic>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>

#include <unistd.h>

namespace test_support {

class TemporaryDirectory {
 public:
  explicit TemporaryDirectory(std::string_view tag) {
    static std::atomic<int> counter{0};
    path_ = std::filesystem::temp_directory_path() /
            ("spodvr_test_" + std::string(tag) + "_" + std::to_string(::getpid()) + "_" + std::to_string(counter++));
    std::filesystem::remove_all(path_);
    std::filesystem::create_directories(path_);
  }
  ~TemporaryDirectory() {
    std::error_code ignored;
    std::filesystem::remove_all(path_, ignored);
  }
  TemporaryDirectory(const TemporaryDirectory&) = delete;
  TemporaryDirectory& operator=(const TemporaryDirectory&) = delete;

  [[nodiscard]] const std::filesystem::path& path() const { return path_; }

  // Writes `text` to the file `name` (relative to the directory, sub-directories created).
  std::filesystem::path write(const std::string& name, std::string_view text) const {
    const std::filesystem::path file = path_ / name;
    std::filesystem::create_directories(file.parent_path());
    std::ofstream(file) << text;
    return file;
  }

 private:
  std::filesystem::path path_;
};

}  // namespace test_support

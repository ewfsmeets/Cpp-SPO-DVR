// Tests: series directories (io/series_writer.*): the manifest, records appended data first and the time
// last, the sizes checked.
#include <gtest/gtest.h>

#include <cstdint>
#include <fstream>
#include <span>
#include <vector>

#include <nlohmann/json.hpp>

#include "io/series_writer.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using test_support::TemporaryDirectory;
namespace fs = std::filesystem;

namespace {

std::vector<DatasetSpec> two_datasets() {
  return {{"values.bin", "<f8", {2, 3}, {"a", "b"}, "probability", "six numbers"},
          {"index.bin", "<i4", {}, {}, "", "one index"}};
}

std::vector<char> bytes_of_file(const fs::path& path) {
  std::ifstream file(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

}  // namespace

TEST(SeriesWriter, WritesTheManifestAndAppendsRecords) {
  TemporaryDirectory dir("series_writer");
  const fs::path series = dir.path() / "series";
  SeriesWriter writer(series, "test", two_datasets(), {{"record_interval_au", 100.0}});
  EXPECT_EQ(writer.record_bytes(), 8u + 6u * 8u + 4u);
  for (int k = 0; k < 2; ++k) {
    std::vector<double> values(6, 0.5 * k);
    const std::int32_t index = 10 + k;
    writer.append(100.0 * k, {std::as_bytes(std::span(values)), std::as_bytes(std::span(&index, 1))});
  }
  EXPECT_EQ(bytes_of_file(series / "values.bin").size(), 2u * 48u);
  EXPECT_EQ(bytes_of_file(series / "index.bin").size(), 2u * 4u);
  EXPECT_EQ(bytes_of_file(series / "times.bin").size(), 2u * 8u);
  std::ifstream manifest_file(series / "manifest.json");
  const nlohmann::json manifest = nlohmann::json::parse(manifest_file);
  EXPECT_EQ(manifest["kind"], "test");
  EXPECT_EQ(manifest["datasets"]["values"]["shape"], nlohmann::json({2, 3}));
  EXPECT_EQ(manifest["datasets"]["index"]["dtype"], "<i4");
  EXPECT_EQ(manifest["record_interval_au"], 100.0);
}

TEST(SeriesWriter, RefusesAValueOfTheWrongSize) {
  TemporaryDirectory dir("series_writer");
  SeriesWriter writer(dir.path() / "series", "test", two_datasets(), nlohmann::ordered_json::object());
  std::vector<double> five(5, 0.0);
  const std::int32_t index = 0;
  EXPECT_THROW(writer.append(0.0, {std::as_bytes(std::span(five)), std::as_bytes(std::span(&index, 1))}),
               std::logic_error);
  EXPECT_THROW(writer.append(0.0, {std::as_bytes(std::span(five))}), std::logic_error);
  EXPECT_THROW(static_cast<void>(dtype_bytes("<f4")), std::invalid_argument);
}

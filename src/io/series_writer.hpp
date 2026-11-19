#pragma once

#include <cstddef>
#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "io/binary_file.hpp"

namespace spodvr {

// One dataset of a series: its file (a literal like "prob_scat_zr.bin"; the dataset's name is the file
// name without ".bin"), its element type ("<f8": double, "<i4": 32-bit int; little-endian), the shape of one
// record (empty for a scalar), the names of the axes, the units and a description.
struct DatasetSpec {
  std::string file;
  std::string dtype;
  std::vector<int> shape;
  std::vector<std::string> axes;
  std::string units;
  std::string description;
};

// A directory of time series (docs/output_files.md, "Series directories"): manifest.json describes the
// datasets; times.bin and one file per dataset hold the records, appended one record at a time. Each
// append writes and flushes the data first and the time last, so a reader that takes
// min(len(times), size / record bytes) records never sees a time without its data.
class SeriesWriter {
 public:
  // Creates `directory` (removing an older series there) and writes manifest.json; `extra` holds
  // further manifest fields (grids, intervals) and `kind` names the series.
  SeriesWriter(const std::filesystem::path& directory, const std::string& kind, std::vector<DatasetSpec> datasets,
               const nlohmann::ordered_json& extra);

  // One record: the bytes of every dataset, in the order of the specs (each must have the spec's size).
  void append(double time, const std::vector<std::span<const std::byte>>& values);

  // Bytes of one record of all datasets and its time.
  [[nodiscard]] std::size_t record_bytes() const;

 private:
  std::vector<DatasetSpec> datasets_;
  std::vector<std::unique_ptr<BinaryWriter>> files_;
  std::unique_ptr<BinaryWriter> times_;
};

// A directory of arrays written once (docs/output_files.md, "Series directories": the same manifest with
// time_dependent false and no times.bin), one file per dataset holding its one record. The caller writes each
// file through file(k), in the dataset's memory order, and then calls finish().
class ArraysWriter {
 public:
  // Creates `directory` (removing an older one there) and writes manifest.json.
  ArraysWriter(const std::filesystem::path& directory, const std::string& kind, std::vector<DatasetSpec> datasets,
               const nlohmann::ordered_json& extra);

  // The file of dataset k (in the order of the specs).
  [[nodiscard]] BinaryWriter& file(std::size_t k);

  // Closes the files; throws std::logic_error when a file does not hold exactly its dataset's bytes.
  void finish();

 private:
  std::filesystem::path directory_;
  std::vector<DatasetSpec> datasets_;
  std::vector<std::unique_ptr<BinaryWriter>> files_;
};

// Bytes of one value of a dtype ("<f8" 8, "<i4" 4); throws std::invalid_argument for others.
[[nodiscard]] std::size_t dtype_bytes(const std::string& dtype);

// Bytes of one record of a dataset: dtype_bytes times the product of its shape.
[[nodiscard]] std::size_t dataset_bytes(const DatasetSpec& spec);

}  // namespace spodvr

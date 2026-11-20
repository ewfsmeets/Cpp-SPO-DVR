// Series directories: a manifest and record-major binary files appended during a run
// (docs/output_files.md, "Series directories"). C++ only (no Fortran counterpart).
#include "io/series_writer.hpp"

#include <fstream>
#include <stdexcept>
#include <utility>

namespace spodvr {

std::size_t dtype_bytes(const std::string& dtype) {
  if (dtype == "<f8") return 8;
  if (dtype == "<i4") return 4;
  if (dtype == "<c16") return 16;
  throw std::invalid_argument("unknown dtype '" + dtype + "'");
}

std::size_t dataset_bytes(const DatasetSpec& spec) {
  std::size_t elements = 1;
  for (const int extent : spec.shape) elements *= static_cast<std::size_t>(extent);
  return elements * dtype_bytes(spec.dtype);
}

namespace {

// The name of a dataset: its file name without ".bin".
std::string dataset_name(const DatasetSpec& spec) { return std::filesystem::path(spec.file).stem().string(); }

// The manifest of a series (time_dependent: times.bin and records) or of arrays written once.
nlohmann::ordered_json manifest(const std::string& kind, const std::vector<DatasetSpec>& datasets,
                                const nlohmann::ordered_json& extra, bool time_dependent) {
  nlohmann::ordered_json m;
  m["format_version"] = 1;
  m["kind"] = kind;
  m["byte_order"] = "little";
  m["memory_order"] = "C (the last axis fastest)";
  m["time_dependent"] = time_dependent;
  if (time_dependent) m["times"] = {{"file", "times.bin"}, {"dtype", "<f8"}, {"units", "a.u. of time"}};
  nlohmann::ordered_json described = nlohmann::ordered_json::object();
  for (const DatasetSpec& d : datasets)
    described[dataset_name(d)] = {{"file", d.file}, {"dtype", d.dtype}, {"shape", d.shape},
                                  {"axes", d.axes}, {"units", d.units}, {"description", d.description}};
  m["datasets"] = described;
  for (const auto& [key, value] : extra.items()) m[key] = value;
  return m;
}

// Creates `directory` afresh (removing an older one) with its manifest.json.
void start_directory(const std::filesystem::path& directory, const nlohmann::ordered_json& manifest) {
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);
  std::ofstream manifest_file(directory / "manifest.json");
  manifest_file << manifest.dump(1) << '\n';
  manifest_file.close();
  if (!manifest_file) throw CheckpointError("error writing '" + (directory / "manifest.json").string() + "'");
}

}  // namespace

SeriesWriter::SeriesWriter(const std::filesystem::path& directory, const std::string& kind,
                           std::vector<DatasetSpec> datasets, const nlohmann::ordered_json& extra)
    : datasets_(std::move(datasets)) {
  start_directory(directory, manifest(kind, datasets_, extra, true));
  for (const DatasetSpec& d : datasets_) files_.push_back(std::make_unique<BinaryWriter>(directory / d.file));
  times_ = std::make_unique<BinaryWriter>(directory / "times.bin");
}

void SeriesWriter::append(double time, const std::vector<std::span<const std::byte>>& values) {
  if (values.size() != datasets_.size()) throw std::logic_error("SeriesWriter::append: one value per dataset");
  for (std::size_t k = 0; k < values.size(); ++k) {
    if (values[k].size() != dataset_bytes(datasets_[k]))
      throw std::logic_error("SeriesWriter::append: " + datasets_[k].file + " has the wrong size");
    files_[k]->write_array(values[k]);
    files_[k]->flush();
  }
  times_->write(time);
  times_->flush();
}

std::size_t SeriesWriter::record_bytes() const {
  std::size_t bytes = sizeof(double);
  for (const DatasetSpec& d : datasets_) bytes += dataset_bytes(d);
  return bytes;
}

ArraysWriter::ArraysWriter(const std::filesystem::path& directory, const std::string& kind,
                           std::vector<DatasetSpec> datasets, const nlohmann::ordered_json& extra)
    : directory_(directory), datasets_(std::move(datasets)) {
  start_directory(directory, manifest(kind, datasets_, extra, false));
  for (const DatasetSpec& d : datasets_) files_.push_back(std::make_unique<BinaryWriter>(directory / d.file));
}

BinaryWriter& ArraysWriter::file(std::size_t k) { return *files_.at(k); }

void ArraysWriter::finish() {
  for (std::size_t k = 0; k < files_.size(); ++k) {
    files_[k]->close();
    const std::filesystem::path path = directory_ / datasets_[k].file;
    if (std::filesystem::file_size(path) != dataset_bytes(datasets_[k]))
      throw std::logic_error("ArraysWriter: " + path.string() + " has not the size of its dataset");
  }
}

}  // namespace spodvr

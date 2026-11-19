// The plot data of a run (densities/, potential/): their projected sizes, the free-space check before the set-up
// and the "PLOT DATA" lines. C++ only.
#include "report/plot_data.hpp"

#include <cstddef>
#include <fstream>

#include <nlohmann/json.hpp>

#include "core/error.hpp"
#include "report/densities_files.hpp"
#include "report/potential_files.hpp"

namespace spodvr {

std::vector<PlotDataSize> plot_data_sizes(const Config& config, const SpatialGrids& grids, const AngularGrid& angular) {
  std::vector<PlotDataSize> sizes;
  if (const auto densities = densities_data_size(config, grids)) sizes.push_back(*densities);
  if (const auto potential = potential_data_size(config, grids, angular)) sizes.push_back(*potential);
  return sizes;
}

std::uintmax_t free_space(const std::filesystem::path& output_dir) {
  std::filesystem::path existing = std::filesystem::absolute(output_dir);
  while (!std::filesystem::exists(existing) && existing.has_parent_path() && existing.parent_path() != existing)
    existing = existing.parent_path();
  return std::filesystem::space(existing).available;
}

std::uintmax_t replaced_plot_data_bytes(const std::filesystem::path& output_dir,
                                        const std::vector<PlotDataSize>& sizes) {
  std::uintmax_t bytes = 0;
  for (const PlotDataSize& s : sizes) {
    const std::filesystem::path replaced = output_dir / s.directory;
    if (!std::filesystem::is_directory(replaced)) continue;
    for (const auto& entry : std::filesystem::recursive_directory_iterator(replaced))
      if (entry.is_regular_file()) bytes += entry.file_size();
  }
  return bytes;
}

namespace {

// Whether `directory` holds a manifest.json of the given kind (plot data this program wrote).
bool is_plot_data(const std::filesystem::path& directory, const std::string& kind) {
  std::ifstream file(directory / "manifest.json");
  if (!file) return false;
  const nlohmann::json manifest = nlohmann::json::parse(file, nullptr, false);
  return manifest.is_object() && manifest.value("kind", std::string()) == kind;
}

}  // namespace

void check_replaceable(const std::filesystem::path& output_dir, const std::vector<PlotDataSize>& sizes) {
  for (const PlotDataSize& s : sizes) {
    const std::filesystem::path directory = output_dir / s.directory;
    const std::string kind = std::filesystem::path(s.directory).parent_path().string();  // "potential/" -> "potential"
    if (std::filesystem::exists(directory) && !is_plot_data(directory, kind))
      throw InputError(s.key, (output_dir / kind).string() +
                                  " exists and is not plot data written by spodvr (the run "
                                  "would remove it); move it or choose another --output-dir");
  }
}

void check_free_space(const std::vector<PlotDataSize>& sizes, std::uintmax_t available) {
  std::uintmax_t total = 0;
  for (const PlotDataSize& s : sizes) total += s.bytes;
  if (sizes.empty() || total <= available) return;
  std::size_t largest = 0;
  for (std::size_t k = 1; k < sizes.size(); ++k)
    if (sizes[k].bytes > sizes[largest].bytes) largest = k;
  throw InputError(sizes[largest].key,
                   "the plot data (" + format_bytes(total) + ") exceed the free space of the output directory (" +
                       format_bytes(available) + "); " + sizes[largest].directory + " is the largest part");
}

void print_plot_data(std::ostream& out, const Style& style, const std::vector<PlotDataSize>& sizes) {
  if (sizes.empty()) return;
  print_section(out, style, "PLOT DATA");
  Table table({"directory", "key", "size", "cost"});
  for (const PlotDataSize& s : sizes) table.add_row({s.directory, s.key, format_bytes(s.bytes), s.cost});
  table.print(out, style);
}

}  // namespace spodvr

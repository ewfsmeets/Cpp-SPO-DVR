// densities/ (writeDensities): the probability per (Z, r) cell at every averages time, a series directory
// (io/series_writer.hpp); the projected size of the plot data and the free-space check. C++ only.
#include "report/densities_files.hpp"

#include <cstddef>
#include <span>

#include "core/error.hpp"

namespace spodvr {

namespace {

// The datasets of densities/, in the order DensitiesFiles::append writes them.
std::vector<DatasetSpec> density_datasets(const Config& config, const SpatialGrids& grids) {
  const int nz = grids.z_scattering.size();
  const int nz_specular = grids.z_specular.size();
  const int nr = config.input.nr;
  return {
      {"prob_scat_zr.bin",
       "<f8",
       {nz, nr},
       {"Z", "r"},
       "probability",
       "probability per (Z, r) cell of the scattering wave (summed over x, y and the rotational states): "
       "|psi|^2 dZ dr dx dy sin(skew)"},
      {"prob_spec_zr.bin",
       "<f8",
       {nz_specular, nr},
       {"Z", "r"},
       "probability",
       "probability per (Z, r) cell of the specular wave (the incoming channel): |chi|^2 dZ dr"},
      {"surface_index.bin", "<i4", {}, {}, "", "0-based Z index of the specular analysis surface"},
      {"norm.bin", "<f8", {}, {}, "probability", "the total norm of averages.dat at this time"},
  };
}

nlohmann::ordered_json density_grids(const Config& config, const SpatialGrids& grids) {
  const InputParameters& p = config.input;
  return {{"record_interval_au", p.dt_averages},
          {"grids",
           {{"z_start", grids.z_scattering.start()},
            {"dz", grids.z_scattering.spacing()},
            {"nz", grids.z_scattering.size()},
            {"nz_specular", grids.z_specular.size()},
            {"r_start", grids.r.start()},
            {"dr", grids.r.spacing()},
            {"nr", p.nr}}}};
}

}  // namespace

DensitiesFiles::DensitiesFiles(const std::filesystem::path& output_dir, const Config& config, const SpatialGrids& grids)
    : series_(output_dir / "densities", "densities", density_datasets(config, grids), density_grids(config, grids)) {}

void DensitiesFiles::append(const Averages& averages, const CellProbabilities& cells, int surface_index) {
  const std::int32_t index = surface_index;
  series_.append(averages.time, {std::as_bytes(std::span(cells.scattering)), std::as_bytes(std::span(cells.specular)),
                                 std::as_bytes(std::span(&index, 1)), std::as_bytes(std::span(&averages.norm, 1))});
}

std::vector<PlotDataSize> plot_data_sizes(const Config& config, const SpatialGrids& grids) {
  std::vector<PlotDataSize> sizes;
  if (config.input.write_densities) {
    std::size_t record = sizeof(double);
    for (const DatasetSpec& d : density_datasets(config, grids)) record += dataset_bytes(d);
    const int analysis_steps = config.input.n_steps / config.derived.steps_per_analysis;
    const int records = 1 + analysis_steps / config.derived.analyses_per_average;  // tStart and every dtAve
    sizes.push_back({"densities/", "writeDensities", static_cast<std::uintmax_t>(records) * record});
  }
  return sizes;
}

void check_free_space(const std::vector<PlotDataSize>& sizes, std::uintmax_t available) {
  std::uintmax_t total = 0;
  for (const PlotDataSize& s : sizes) total += s.bytes;
  if (sizes.empty() || total <= available) return;
  std::size_t largest = 0;
  for (std::size_t k = 1; k < sizes.size(); ++k)
    if (sizes[k].bytes > sizes[largest].bytes) largest = k;
  throw InputError(sizes[largest].key, "the plot data (" + format_bytes(total) +
                                           ") exceed the free space of the "
                                           "output directory (" +
                                           format_bytes(available) + "); " + sizes[largest].directory +
                                           " is the largest part");
}

void print_plot_data(std::ostream& out, const Style& style, const std::vector<PlotDataSize>& sizes) {
  if (sizes.empty()) return;
  print_section(out, style, "PLOT DATA");
  Table table({"directory", "key", "size", "cost in the time loop"});
  for (const PlotDataSize& s : sizes)
    table.add_row({s.directory, s.key, format_bytes(s.bytes), "a file write per averages row (no extra reduction)"});
  table.print(out, style);
}

}  // namespace spodvr

// densities/ (writeDensities): the probability per (Z, r) cell at every averages time, a series directory
// (io/series_writer.hpp), and its projected size. C++ only.
#include "report/densities_files.hpp"

#include <array>
#include <cstddef>
#include <span>

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
       "probability per (Z, r) cell of the specular wave (the initial channel, incoming and specularly "
       "reflected, every v): |chi|^2 dZ dr"},
      {"surface_index.bin", "<i4", {}, {}, "", "0-based Z index of the specular analysis surface"},
      {"norm.bin", "<f8", {}, {}, "probability", "the total norm of averages.dat at this time"},
  };
}

nlohmann::ordered_json density_grids(const Config& config, const SpatialGrids& grids) {
  const InputParameters& p = config.input;
  // The averaging interval the program uses: its time steps per average (Fortran nave * nanspo, truncated counts)
  // times dtStep. Not always dtAve: dtAve just below m dtAna gives m - 1 analyses per average (analysis D4), and the
  // product need not be bitwise the key (3 * 0.3 is 0.8999999999999999).
  const double steps = static_cast<double>(config.derived.analyses_per_average) * config.derived.steps_per_analysis;
  const double record_interval = steps * config.input.time_step;
  return {{"record_interval_au", record_interval},
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
  const std::array<std::int32_t, 1> index{surface_index};
  const std::array<double, 1> norm{averages.norm};
  series_.append(averages.time, {std::as_bytes(std::span(cells.scattering)), std::as_bytes(std::span(cells.specular)),
                                 std::as_bytes(std::span(index)), std::as_bytes(std::span(norm))});
}

std::optional<PlotDataSize> densities_data_size(const Config& config, const SpatialGrids& grids) {
  if (!config.input.write_densities) return std::nullopt;
  std::size_t record = sizeof(double);
  for (const DatasetSpec& d : density_datasets(config, grids)) record += dataset_bytes(d);
  const int analysis_steps = config.input.n_steps / config.derived.steps_per_analysis;
  const int records = 1 + analysis_steps / config.derived.analyses_per_average;  // tStart and every dtAve
  return PlotDataSize{"densities/", "writeDensities", static_cast<std::uintmax_t>(records) * record,
                      "a file write per averages row (no extra reduction)"};
}

}  // namespace spodvr

// frames/ (dtFrames): one record of the wave packet's reductions every dtFrames, a series directory, its projected
// size and cost, and its measured cost. C++ only.
#include "report/frame_files.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <format>
#include <span>
#include <string>

#include "frames/frame_source.hpp"

namespace spodvr {

namespace {

// Measured on the thesis DZ grid with 8 threads (the plot plan's measurements): seconds per complex value of
// psi (nZ nR nX nY nrot) for a frame (the time step's: report/plot_data.hpp).
constexpr double frame_seconds_per_value_cpu = 0.52e-9;
constexpr double frame_seconds_per_value_gpu = 2.45e-9;

// The datasets of frames/, in the order FrameFiles::append writes them.
std::vector<DatasetSpec> frame_datasets(const FrameLayout& layout, const Config& config) {
  const InputParameters& p = config.input;
  return {
      {"prob_zrxy.bin",
       "<f8",
       {layout.nz, layout.nr, p.nx, p.ny},
       {"Z", "r", "x", "y"},
       "probability",
       "probability per (Z, r, x, y) cell of psi_full = psi_scat + s2s chi (chi in the initial channel), summed "
       "over the rotational states: |psi_full|^2 dZ dr dx dy sin(skew)"},
      {"prob_tail_zr.bin",
       "<f8",
       {layout.nz_specular - layout.nz, layout.nr},
       {"Z", "r"},
       "probability",
       "probability per (Z, r) cell of the specular wave beyond the scattering grid (Z index nZ on): |chi|^2 dZ dr"},
      {"pop_jmj_z.bin",
       "<f8",
       {layout.nz, layout.n_rot},
       {"Z", "rotational state"},
       "probability",
       "probability per (Z, rotational state (J, mJ)) of psi_full, summed over r, x, y"},
      {"surface_index.bin", "<i4", {}, {}, "", "0-based Z index of the specular analysis surface"},
  };
}

nlohmann::ordered_json frame_manifest_fields(const Config& config, const SpatialGrids& grids,
                                             const RotationalBasis& basis) {
  const InputParameters& p = config.input;
  nlohmann::ordered_json states = nlohmann::ordered_json::array();
  for (int irot = 0; irot < basis.size(); ++irot) states.push_back({basis.j(irot), basis.mj(irot)});
  return {{"record_interval_au", *p.dt_frames},
          {"grids",
           {{"z_start", grids.z_scattering.start()},
            {"dz", grids.z_scattering.spacing()},
            {"nz", grids.z_scattering.size()},
            {"nz_specular", grids.z_specular.size()},
            {"r_start", grids.r.start()},
            {"dr", grids.r.spacing()},
            {"nr", p.nr},
            {"nx", p.nx},
            {"ny", p.ny}}},
          {"rotational_states", states},
          {"initial_index", basis.initial_index()},
          {"sites", plot_sites_json(p)}};
}

}  // namespace

FrameFiles::FrameFiles(const std::filesystem::path& output_dir, const Config& config, const SpatialGrids& grids,
                       const RotationalBasis& basis, const FrameLayout& layout)
    : series_(output_dir / "frames", "frames", frame_datasets(layout, config),
              frame_manifest_fields(config, grids, basis)) {}

void FrameFiles::append(double time, const FrameRecord& record, int surface_index) {
  const std::array<std::int32_t, 1> index{surface_index};
  series_.append(time, {std::as_bytes(std::span(record.prob_zrxy)), std::as_bytes(std::span(record.prob_tail_zr)),
                        std::as_bytes(std::span(record.pop_jmj_z)), std::as_bytes(std::span(index))});
}

std::optional<PlotDataSize> frames_data_size(const Config& config, const SpatialGrids& grids,
                                             const RotationalBasis& basis) {
  if (!config.input.dt_frames) return std::nullopt;
  const FrameLayout layout = frame_layout(config, grids, basis, {});
  std::size_t record = sizeof(double);
  for (const DatasetSpec& d : frame_datasets(layout, config)) record += dataset_bytes(d);
  const int analysis_steps = config.input.n_steps / config.derived.steps_per_analysis;
  const int records = 1 + analysis_steps / config.derived.analyses_per_frame;  // tStart and every dtFrames
  const double values = double(layout.nz) * layout.nr * layout.nxy * layout.n_rot;
  const double steps_between = double(config.derived.analyses_per_frame) * config.derived.steps_per_analysis;
  const double cpu = frame_seconds_per_value_cpu * values;
  const double gpu = frame_seconds_per_value_gpu * values;
  const double cpu_share = cpu / (steps_between * step_seconds_per_value_cpu * values);
  const double gpu_share = gpu / (steps_between * step_seconds_per_value_gpu * values);
  std::string cost = "about " + fixed(cpu, 2) + " s (CPU) or " + fixed(gpu, 2) + " s (GPU) per frame, " +
                     percent(cpu_share) + " or " + percent(gpu_share) + " of the propagation between frames";
  if (std::max(cpu_share, gpu_share) > cost_warning_fraction) cost += "; above 2 %: a larger dtFrames is cheaper";
  return PlotDataSize{"frames/", "dtFrames", static_cast<std::uintmax_t>(records) * record, cost};
}

void print_frame_cost(std::ostream& out, int frames_written, double seconds, double propagation_seconds) {
  if (frames_written == 0) return;
  out << "  Frames (dtFrames): " << frames_written << " written, " << fixed(seconds, 2) << " s in all, "
      << fixed(seconds / frames_written, 3) << " s per frame (" << percent(seconds / propagation_seconds)
      << " of the propagation time)\n";
}

}  // namespace spodvr

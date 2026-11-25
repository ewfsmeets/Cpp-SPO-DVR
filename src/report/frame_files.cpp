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
#include "input/validate_outputs.hpp"

namespace spodvr {

namespace {

// Measured on the thesis DZ grid with 8 threads (the plot plan's measurements): seconds per complex value of
// psi (nZ nR nX nY nrot) for a frame (the time step's: report/plot_data.hpp).
constexpr double frame_seconds_per_value_cpu = 0.52e-9;
constexpr double frame_seconds_per_value_gpu = 2.45e-9;

// FrameExtras, measured on the same grid (200 steps, 6 frames, every extra on; 2026-11-20): host seconds per complex
// value of psi per frame (the extras run on the host, so both backends pay them); the vibration extra's per value and
// r point (each value enters nR overlaps). The angular extra measured again after its code changed (2 CPU runs,
// 2026-11-25: 0.84 s per frame instead of 1.50 s); the others then within 5 % of these values, their code unchanged.
constexpr double angular_seconds_per_value = 2.11e-9;
constexpr double diffraction_seconds_per_value = 1.15e-9;
constexpr double vibration_seconds_per_value_and_r = 0.132e-9;
constexpr double rotation_moments_seconds_per_value = 0.88e-9;
constexpr double momentum_seconds_per_value = 4.89e-9;

// The projected host seconds per frame of the extras that are on.
double extras_seconds(const FrameExtraSwitches& on, double values, int nr) {
  double seconds = 0.0;
  if (on.angular) seconds += angular_seconds_per_value * values;
  if (on.diffraction) seconds += diffraction_seconds_per_value * values;
  if (on.vibration) seconds += vibration_seconds_per_value_and_r * values * nr;
  if (on.rotation_moments) seconds += rotation_moments_seconds_per_value * values;
  if (on.momentum) seconds += momentum_seconds_per_value * values;
  return seconds;
}

// The datasets of FrameExtras that are on, in the order FrameFiles::append writes them after the core ones.
std::vector<DatasetSpec> extra_datasets(const FrameLayout& layout, const Config& config, const AngularGrid& angular,
                                        const RotationalBasis& basis) {
  const InputParameters& p = config.input;
  const FrameExtraSwitches on = FrameExtraSwitches::of(p);
  std::vector<DatasetSpec> specs;
  if (on.angular)
    specs.push_back({"prob_angle_z.bin",
                     "<f8",
                     {layout.nz, angular.n_theta(), angular.n_phi()},
                     {"Z", "theta node", "phi node"},
                     "probability",
                     "probability per (Z, DVR node (theta_t, phi_k)) of psi_full, summed over r, x, y; divided by "
                     "theta_weights[t] 2 pi / nphi it is a density per steradian"});
  if (on.diffraction)
    specs.push_back(
        {"pop_nm_z.bin",
         "<f8",
         {layout.nz, p.nx, p.ny},
         {"Z", "n", "m"},
         "probability",
         "probability per (Z, parallel momentum k_par + n b1 + m b2) of psi_full (indices in FFT order: "
         "n for n <= nX/2, else n - nX, as the analysis labels them), summed over r and the rotational states"});
  if (on.vibration)
    specs.push_back({"pop_vj_z.bin",
                     "<f8",
                     {layout.nz, layout.nr, static_cast<int>(basis_j_values(basis).size())},
                     {"Z", "v", "J"},
                     "probability",
                     "probability per (Z, gas-phase state phi_vJ(r)) of psi_full, all nR states of each J of the basis "
                     "(vibration_j_values), summed over x, y and mJ"});
  if (on.rotation_moments) {
    specs.push_back({"p2_numerator_z.bin",
                     "<f8",
                     {layout.nz},
                     {"Z"},
                     "probability",
                     "<P2(cos theta)> times the probability at each Z (psi_full; analytic J, J +- 2 elements)"});
    specs.push_back({"rot_numerator_z.bin",
                     "<f8",
                     {layout.nz},
                     {"Z"},
                     "hartree",
                     "<J(J+1) / (2 mu r^2)> times the probability at each Z (psi_full)"});
  }
  if (on.momentum) {
    specs.push_back(
        {"prob_kr.bin",
         "<f8",
         {layout.nr},
         {"k_r index"},
         "probability",
         "probability per k_r point (FFT order, spacing 2 pi / (nR dr)) of psi_full on the scattering grid"});
    specs.push_back({"prob_kz_spec.bin",
                     "<f8",
                     {layout.nz_specular},
                     {"k_Z index"},
                     "probability",
                     "probability per k_Z point (FFT order, spacing 2 pi / (nZSp dZ)) of the specular wave"});
  }
  return specs;
}

// The datasets of frames/, in the order FrameFiles::append writes them.
std::vector<DatasetSpec> frame_datasets(const FrameLayout& layout, const Config& config, const AngularGrid& angular,
                                        const RotationalBasis& basis) {
  const InputParameters& p = config.input;
  std::vector<DatasetSpec> specs = {
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
  for (DatasetSpec& spec : extra_datasets(layout, config, angular, basis)) specs.push_back(std::move(spec));
  return specs;
}

// The manifest fields of the extras: their names, the DVR grid (angular), the J values (vibration).
nlohmann::ordered_json extras_manifest(const Config& config, const AngularGrid& angular, const RotationalBasis& basis) {
  const FrameExtraSwitches on = FrameExtraSwitches::of(config.input);
  nlohmann::ordered_json extras = {{"names", config.input.frame_extras}};
  if (on.angular)
    extras["angular"] = {
        {"theta_rad", std::vector<double>(angular.thetas().begin(), angular.thetas().end())},
        {"theta_weights", std::vector<double>(angular.theta_weights().begin(), angular.theta_weights().end())},
        {"phi_rad", std::vector<double>(angular.phis().begin(), angular.phis().end())},
        {"homonuclear", angular.homonuclear()}};
  if (on.vibration) extras["vibration_j_values"] = basis_j_values(basis);
  return extras;
}

nlohmann::ordered_json frame_manifest_fields(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, const RotationalBasis& basis) {
  const InputParameters& p = config.input;
  nlohmann::ordered_json states = nlohmann::ordered_json::array();
  for (int irot = 0; irot < basis.size(); ++irot) states.push_back({basis.j(irot), basis.mj(irot)});
  // The frame interval the program uses: its time steps per frame (analyses per frame * nanspo) times dtStep; dtFrames
  // to the tolerances of the dtAna and dtFrames checks.
  const double steps = static_cast<double>(config.derived.analyses_per_frame) * config.derived.steps_per_analysis;
  const double record_interval = steps * config.input.time_step;
  nlohmann::ordered_json fields = {{"record_interval_au", record_interval},
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
  // Only with FrameExtras: frames/ of the other runs keep their manifest.
  if (!p.frame_extras.empty()) fields["extras"] = extras_manifest(config, angular, basis);
  return fields;
}

}  // namespace

FrameFiles::FrameFiles(const std::filesystem::path& output_dir, const Config& config, const SpatialGrids& grids,
                       const AngularGrid& angular, const RotationalBasis& basis, const FrameLayout& layout)
    : series_(output_dir / "frames", "frames", frame_datasets(layout, config, angular, basis),
              frame_manifest_fields(config, grids, angular, basis)) {}

void FrameFiles::append(double time, const FrameRecord& record, int surface_index) {
  const std::array<std::int32_t, 1> index{surface_index};
  std::vector<std::span<const std::byte>> values = {
      std::as_bytes(std::span(record.prob_zrxy)), std::as_bytes(std::span(record.prob_tail_zr)),
      std::as_bytes(std::span(record.pop_jmj_z)), std::as_bytes(std::span(index))};
  // The extras that are on, in the order of extra_datasets (an extra that is off has empty vectors).
  const std::vector<std::span<const double>> extras = {record.prob_angle_z,   record.pop_nm_z,        record.pop_vj_z,
                                                       record.p2_numerator_z, record.rot_numerator_z, record.prob_kr,
                                                       record.prob_kz_spec};
  for (const std::span<const double> extra : extras)
    if (!extra.empty()) values.push_back(std::as_bytes(extra));
  series_.append(time, values);
}

std::optional<PlotDataSize> frames_data_size(const Config& config, const SpatialGrids& grids,
                                             const AngularGrid& angular, const RotationalBasis& basis) {
  if (!config.input.dt_frames) return std::nullopt;
  const FrameLayout layout = frame_layout(config, grids, basis, {});
  std::size_t record = sizeof(double);
  for (const DatasetSpec& d : frame_datasets(layout, config, angular, basis)) record += dataset_bytes(d);
  const int analysis_steps = config.input.n_steps / config.derived.steps_per_analysis;
  const int records = 1 + analysis_steps / config.derived.analyses_per_frame;  // tStart and every dtFrames
  const double values = double(layout.nz) * layout.nr * layout.nxy * layout.n_rot;
  const double steps_between = double(config.derived.analyses_per_frame) * config.derived.steps_per_analysis;
  const double extras = extras_seconds(FrameExtraSwitches::of(config.input), values, layout.nr);
  const double cpu = frame_seconds_per_value_cpu * values + extras;
  const double gpu = frame_seconds_per_value_gpu * values + extras;
  const double cpu_share = cpu / (steps_between * step_seconds_per_value_cpu * values);
  const double gpu_share = gpu / (steps_between * step_seconds_per_value_gpu * values);
  std::string cost = "about " + fixed(cpu, 2) + " s (CPU) or " + fixed(gpu, 2) + " s (GPU) per frame, " +
                     percent(cpu_share) + " or " + percent(gpu_share) + " of the propagation between frames";
  if (std::max(cpu_share, gpu_share) > cost_warning_fraction) cost += "; above 2 %: a larger dtFrames is cheaper";
  if (extras > 0.0) cost += "; FrameExtras " + fixed(extras, 2) + " s of it (on the host)";
  return PlotDataSize{"frames/", "dtFrames", static_cast<std::uintmax_t>(records) * record, cost};
}

void print_frame_extras_cost(std::ostream& out, int frames_written, const std::vector<double>& seconds) {
  if (frames_written == 0 || seconds.empty()) return;
  const std::vector<std::string>& names = frame_extra_names();  // the order of FrameExtras::seconds
  out << "  Frame extras (FrameExtras), seconds per frame:";
  for (std::size_t k = 0; k < names.size(); ++k)
    if (seconds[k] > 0.0) out << " " << names[k] << " " << fixed(seconds[k] / frames_written, 3);
  out << '\n';
}

void print_frame_cost(std::ostream& out, int frames_written, double seconds, double propagation_seconds) {
  if (frames_written == 0) return;
  out << "  Frames (dtFrames): " << frames_written << " written, " << fixed(seconds, 2) << " s in all, "
      << fixed(seconds / frames_written, 3) << " s per frame (" << percent(seconds / propagation_seconds)
      << " of the propagation time)\n";
}

}  // namespace spodvr

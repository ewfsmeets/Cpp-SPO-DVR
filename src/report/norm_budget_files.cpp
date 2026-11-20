// norm_budget/ (NormBudget): the losses of every absorber per analysis window, a series directory, and its projected
// size and cost. C++ only.
#include "report/norm_budget_files.hpp"

#include <algorithm>
#include <array>
#include <span>
#include <string>

#include "report/format.hpp"

namespace spodvr {

namespace {

// Measured on the thesis DZ grid with 8 threads (the plot plan's measurements): seconds per complex value of psi
// (nZ nR nX nY nrot) for one block-density reduction.
constexpr double reduction_seconds_per_value_cpu = 1.25e-9;
constexpr double reduction_seconds_per_value_gpu = 0.28e-9;

// The absorbers of the budget as the configuration has them (NormBudget::absorbers lists them in this number).
int absorber_count(const InputParameters& p) {
  return 4 + (p.z_opt_low.enabled() ? 1 : 0) + (p.zsp_opt_low.enabled() ? 1 : 0);
}

// The datasets of norm_budget/, in the order NormBudgetFiles::append writes them.
std::vector<DatasetSpec> budget_datasets(int n_absorbers) {
  return {
      {"loss.bin",
       "<f8",
       {n_absorbers},
       {"absorber"},
       "probability",
       "what each absorber removes in the window, summed over its cells, each exact for a constant inflow: "
       "rho_start (1 - F) + (rho_end - F rho_start) (gamma / (1 - F) - 1), F = exp(-gamma) the density's factor"},
      {"loss_error.bin",
       "<f8",
       {n_absorbers},
       {"absorber"},
       "probability",
       "summed over the cells, gamma^2 / 6 of the inflow (rho_end - rho_start + loss), the leading distance of the "
       "loss from the mean of rho_start (1 - F) and rho_end (1/F - 1): falls with dtAna^2"},
      {"norm_start.bin", "<f8", {}, {}, "probability", "the probability on both grids at the window's start"},
      {"norm_end.bin", "<f8", {}, {}, "probability", "the probability on both grids at the window's end"},
  };
}

nlohmann::ordered_json budget_manifest_fields(const Config& config, const std::vector<BudgetAbsorber>& absorbers) {
  nlohmann::ordered_json list = nlohmann::ordered_json::array();
  for (const BudgetAbsorber& a : absorbers) list.push_back({{"name", a.name}, {"grid", a.grid}});
  return {{"record_interval_au", config.input.dt_analysis},
          {"window", {{"steps", config.derived.steps_per_analysis}, {"dt_step", config.input.time_step}}},
          {"absorbers", list}};
}

}  // namespace

NormBudgetFiles::NormBudgetFiles(const std::filesystem::path& output_dir, const Config& config,
                                 const std::vector<BudgetAbsorber>& absorbers)
    : series_(output_dir / "norm_budget", "norm_budget", budget_datasets(static_cast<int>(absorbers.size())),
              budget_manifest_fields(config, absorbers)) {}

void NormBudgetFiles::append(double time, const NormBudgetWindow& window) {
  const std::array<double, 1> norm_start{window.norm_start};
  const std::array<double, 1> norm_end{window.norm_end};
  series_.append(time, {std::as_bytes(std::span(window.loss)), std::as_bytes(std::span(window.loss_error)),
                        std::as_bytes(std::span(norm_start)), std::as_bytes(std::span(norm_end))});
}

void print_norm_budget_cost(std::ostream& out, int windows_written, double seconds, double propagation_seconds) {
  if (windows_written == 0) return;
  out << "  Norm budget (NormBudget): " << windows_written << " windows written, " << fixed(seconds, 2)
      << " s in all for the densities of every analysis step (" << percent(seconds / propagation_seconds)
      << " of the propagation time, the averages' own included)\n";
}

std::optional<PlotDataSize> norm_budget_data_size(const Config& config) {
  if (!config.input.norm_budget) return std::nullopt;
  std::size_t record = sizeof(double);
  for (const DatasetSpec& d : budget_datasets(absorber_count(config.input))) record += dataset_bytes(d);
  const int records = config.input.n_steps / config.derived.steps_per_analysis;  // one per analysis window
  // A reduction per analysis step against the steps' propagation (both scale with psi's size); the averages rows
  // (every analyses_per_average-th step) reduce anyway.
  const double extra = 1.0 - 1.0 / config.derived.analyses_per_average;
  const double steps = config.derived.steps_per_analysis;
  const double cpu_share = extra * reduction_seconds_per_value_cpu / (steps * step_seconds_per_value_cpu);
  const double gpu_share = extra * reduction_seconds_per_value_gpu / (steps * step_seconds_per_value_gpu);
  std::string cost = "a block-density reduction every analysis step, about " + percent(cpu_share) + " (CPU) or " +
                     percent(gpu_share) + " (GPU) of the propagation";
  if (std::max(cpu_share, gpu_share) > cost_warning_fraction) cost += "; above 2 %: a larger dtAna is cheaper";
  return PlotDataSize{"norm_budget/", "NormBudget", static_cast<std::uintmax_t>(records) * record, cost};
}

}  // namespace spodvr

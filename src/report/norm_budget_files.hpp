#pragma once
// norm_budget/ (NormBudget, C++ only): what every absorber removes per analysis window, a series directory
// (io/series_writer.hpp), and its projected size and cost (docs/output_files.md, "norm_budget/").

#include <filesystem>
#include <optional>
#include <ostream>
#include <vector>

#include "analysis/norm_budget.hpp"
#include "input/config.hpp"
#include "io/series_writer.hpp"
#include "report/plot_data.hpp"

namespace spodvr {

class NormBudgetFiles {
 public:
  NormBudgetFiles(const std::filesystem::path& output_dir, const Config& config,
                  const std::vector<BudgetAbsorber>& absorbers);
  // The window that ends at `time`.
  void append(double time, const NormBudgetWindow& window);

 private:
  SeriesWriter series_;
};

// The measured cost of the norm budget at the end of the run: the windows written and the wall-clock time of the
// block densities of every analysis step and the budget's work on them (the averages rows' own reductions included).
void print_norm_budget_cost(std::ostream& out, int windows_written, double seconds, double propagation_seconds);

// The size of norm_budget/ (a record per analysis step) and its projected cost, or nothing without NormBudget. The
// cost is the block-density reduction of every analysis step that has no averages row (measured on the thesis DZ
// grid: 1.46 ns per complex value of psi on the CPU, 0.17 ns on the GPU).
[[nodiscard]] std::optional<PlotDataSize> norm_budget_data_size(const Config& config);

}  // namespace spodvr

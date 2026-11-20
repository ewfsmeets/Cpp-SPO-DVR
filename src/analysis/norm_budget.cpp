// The norm budget (NormBudget): per analysis window, the probability each absorber removes, from the cell
// probabilities at the window's start and end and the absorbers' factors per time step. C++ only.
#include "analysis/norm_budget.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace spodvr {

namespace {

// Appends an absorber to the list and returns its index.
int add_absorber(std::vector<BudgetAbsorber>& absorbers, const std::string& name, const std::string& grid) {
  absorbers.push_back({name, grid});
  return static_cast<int>(absorbers.size()) - 1;
}

bool inside(const std::optional<IndexRange>& range, int i) { return range && i >= range->first && i <= range->last; }

// A cell's loss over a window (gamma = 2n V_opt dt, its density's exponent) for a constant inflow during the window,
// and its error indicator, gamma^2 / 6 of the inflow (NormBudgetWindow).
struct CellLoss {
  double loss;
  double error;
};
CellLoss cell_loss(double gamma, double rho_start, double rho_end) {
  const double one_minus_f = -std::expm1(-gamma);  // 1 - F, accurate also where F is close to 1
  const double f = std::exp(-gamma);
  const double inflow_weight = gamma / one_minus_f - 1.0;  // gamma / 2 for a weak absorber
  const double loss = rho_start * one_minus_f + (rho_end - f * rho_start) * inflow_weight;
  const double inflow = rho_end - rho_start + loss;  // what flowed into the cell during the window
  return {loss, gamma * gamma / 6.0 * std::abs(inflow)};
}

}  // namespace

NormBudget::NormBudget(const Config& config, const SpatialGrids& grids, const PropagationTables& tables)
    : nr_(config.input.nr), steps_per_window_(config.derived.steps_per_analysis) {
  const int z_opt = add_absorber(absorbers_, "ZOpt", "scattering");
  const int z_opt_low = grids.z_optical_low ? add_absorber(absorbers_, "ZOptLow", "scattering") : -1;
  const int r_opt = add_absorber(absorbers_, "ROpt", "scattering");
  const int zsp_opt = add_absorber(absorbers_, "ZSpOpt", "specular");
  const int zsp_opt_low = grids.z_specular_optical_low ? add_absorber(absorbers_, "ZSpOptLow", "specular") : -1;
  const int r_opt_specular = add_absorber(absorbers_, "ROpt", "specular");
  // The Z tables hold the high- and the low-end absorber of their grid; the r table acts on both grids.
  scattering_ = absorption(tables.z_optical, tables.r_optical, z_opt, z_opt_low, grids.z_optical_low, r_opt);
  specular_ = absorption(tables.z_specular_optical, tables.r_optical, zsp_opt, zsp_opt_low,
                         grids.z_specular_optical_low, r_opt_specular);
}

NormBudget::GridAbsorption NormBudget::absorption(const std::vector<double>& z_factors,
                                                  const std::vector<double>& r_factors, int z_high, int z_low,
                                                  const std::optional<IndexRange>& low_range, int r_index) const {
  const std::size_t cells = z_factors.size() * nr_;
  GridAbsorption g{std::vector<double>(cells, 0.0), std::vector<double>(cells, 0.0), std::vector<int>(cells, -1),
                   std::vector<int>(cells, -1)};
  for (int iz = 0; iz < static_cast<int>(z_factors.size()); ++iz) {
    for (int ir = 0; ir < nr_; ++ir) {
      const double v_z = -std::log(z_factors[iz]);  // V_opt dt of the Z absorber (0 outside it)
      const double v_r = -std::log(r_factors[ir]);
      if (v_z == 0.0 && v_r == 0.0) continue;
      const std::size_t k = ir + std::size_t(nr_) * iz;
      g.window_exponent[k] = 2.0 * steps_per_window_ * (v_z + v_r);  // the density shrinks by exp(-it)
      if (v_z != 0.0) {
        g.z_absorber[k] = inside(low_range, iz) ? z_low : z_high;
        g.z_share[k] = v_z / (v_z + v_r);
      }
      if (v_r != 0.0) g.r_absorber[k] = r_index;
    }
  }
  return g;
}

void NormBudget::add_losses(const GridAbsorption& g, const std::vector<double>& start, const std::vector<double>& end,
                            NormBudgetWindow& window) const {
  if (start.size() != g.window_exponent.size() || end.size() != g.window_exponent.size())
    throw std::logic_error("NormBudget: cell probabilities of another grid size");
  for (std::size_t k = 0; k < g.window_exponent.size(); ++k) {
    if (g.z_absorber[k] < 0 && g.r_absorber[k] < 0) continue;
    const CellLoss cell = cell_loss(g.window_exponent[k], start[k], end[k]);
    if (g.z_absorber[k] >= 0) {
      window.loss[g.z_absorber[k]] += g.z_share[k] * cell.loss;
      window.loss_error[g.z_absorber[k]] += g.z_share[k] * cell.error;
    }
    if (g.r_absorber[k] >= 0) {
      const double r_share = 1.0 - g.z_share[k];
      window.loss[g.r_absorber[k]] += r_share * cell.loss;
      window.loss_error[g.r_absorber[k]] += r_share * cell.error;
    }
  }
}

std::optional<NormBudgetWindow> NormBudget::add(const CellProbabilities& cells) {
  std::optional<NormBudgetWindow> window;
  if (previous_) {
    NormBudgetWindow w{std::vector<double>(absorbers_.size(), 0.0), std::vector<double>(absorbers_.size(), 0.0), 0.0,
                       0.0};
    add_losses(scattering_, previous_->scattering, cells.scattering, w);
    add_losses(specular_, previous_->specular, cells.specular, w);
    for (const double p : previous_->scattering) w.norm_start += p;
    for (const double p : previous_->specular) w.norm_start += p;
    for (const double p : cells.scattering) w.norm_end += p;
    for (const double p : cells.specular) w.norm_end += p;
    window = std::move(w);
  }
  previous_ = cells;
  return window;
}

std::size_t NormBudget::memory_bytes(const SpatialGrids& grids, int nr) {
  const std::size_t cells = std::size_t(grids.z_scattering.size() + grids.z_specular.size()) * nr;
  const std::size_t scattering_cells = std::size_t(grids.z_scattering.size()) * nr;
  return cells * (4 * sizeof(double) + 2 * sizeof(int)) + 2 * scattering_cells * sizeof(double);
}

}  // namespace spodvr

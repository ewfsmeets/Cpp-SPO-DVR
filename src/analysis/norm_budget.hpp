#pragma once
// The norm budget (NormBudget): per analysis window, the probability each absorber removes, estimated from the cell
// probabilities at the window's start and end. C++ only (the factors are those of Absorption.f90 EXPOpt*).

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/averages.hpp"
#include "propagation/propagation_tables.hpp"

namespace spodvr {

// One absorber of the budget: its input name (ZOpt, ZOptLow, ROpt, ZSpOpt, ZSpOptLow) and the grid it acts on
// ("scattering" or "specular"; ROpt acts on both and is one entry per grid).
struct BudgetAbsorber {
  std::string name;
  std::string grid;
};

// The losses of one analysis window (t - dtAna, t], one per absorber in the order of NormBudget::absorbers. In a
// cell whose absorber factor per time step is f (exp(-V_opt dt)), the density decays at the rate kappa = -2 ln f / dt
// while probability flows in and out; over the window, gamma = kappa dtAna and F = exp(-gamma). With the inflow J
// constant during the window (d rho/dt = -kappa rho + J), the cell loses exactly
//   loss = rho_start (1 - F) + (rho_end - F rho_start) (gamma / (1 - F) - 1),
// which is exact without motion too (rho_end = F rho_start). loss_error: gamma^2 / 6 of the inflow
// (rho_end - rho_start + loss) summed over the cells, the leading distance between the loss and the mean of the two
// one-sided estimates rho_start (1 - F) and rho_end (1/F - 1): it falls with dtAna squared and indicates how well
// the window resolves the motion (it stays finite for long windows, where 1/F does not). norm_start, norm_end: the
// probability on both grids at the window's start and end (their difference is what the absorbers took, so it
// checks the losses).
struct NormBudgetWindow {
  std::vector<double> loss;
  std::vector<double> loss_error;
  double norm_start = 0.0;
  double norm_end = 0.0;
};

// The norm budget of a propagation: fed the cell probabilities after every analysis step (and once at tStart), it
// returns the losses of each window. Where a Z absorber and the r absorber act on the same cell, the cell's loss is
// shared in the ratio of their V_opt (V_opt dt = -ln f).
class NormBudget {
 public:
  NormBudget(const Config& config, const SpatialGrids& grids, const PropagationTables& tables);

  [[nodiscard]] const std::vector<BudgetAbsorber>& absorbers() const { return absorbers_; }

  // The window that ends at `cells` and starts at the cells of the previous call; nothing at the first call.
  [[nodiscard]] std::optional<NormBudgetWindow> add(const CellProbabilities& cells);

  // The bytes the budget needs during the propagation: the absorption of every cell of both grids and the cells of
  // the previous and the current call (2 + 2 doubles and 2 ints per cell), and the time loop's block densities of
  // tStart and of the current analysis step (2 doubles per scattering cell).
  [[nodiscard]] static std::size_t memory_bytes(const SpatialGrids& grids, int nr);

 private:
  // How one grid's cells absorb, per cell (izr = ir + nr iz): the density's exponent over the window, 2n V_opt dt
  // (F = exp(-exponent) = f^2n), the absorber the Z part goes to (-1: none) and its share, and the absorber of the r
  // part (-1: none; it takes the rest).
  struct GridAbsorption {
    std::vector<double> window_exponent;
    std::vector<double> z_share;
    std::vector<int> z_absorber;
    std::vector<int> r_absorber;
  };

  [[nodiscard]] GridAbsorption absorption(const std::vector<double>& z_factors, const std::vector<double>& r_factors,
                                          int z_high, int z_low, const std::optional<IndexRange>& low_range,
                                          int r_index) const;
  void add_losses(const GridAbsorption& grid, const std::vector<double>& start, const std::vector<double>& end,
                  NormBudgetWindow& window) const;

  int nr_;
  int steps_per_window_;
  std::vector<BudgetAbsorber> absorbers_;
  GridAbsorption scattering_;
  GridAbsorption specular_;
  std::optional<CellProbabilities> previous_;
};

}  // namespace spodvr

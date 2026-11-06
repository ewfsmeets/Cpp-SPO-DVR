#pragma once

#include <chrono>
#include <filesystem>
#include <ostream>
#include <string>
#include <vector>

#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/averages.hpp"
#include "propagation/propagation.hpp"
#include "report/format.hpp"

namespace spodvr {

// Rotational basis, initial wave packet and analysis surface of a propagation run.
void print_propagation_setup(std::ostream& out, const Style& style, const SpatialGrids& grids,
                             const WavePacketPropagation& propagation);

// The averages table, printed row by row while the propagation runs (Fortran: averages.out), and
// with colour (a terminal) a live progress line below it, redrawn every analysis step.
class AveragesTable {
 public:
  // The run: n_steps time steps of time_step from t_start (a restart counts from its own start).
  AveragesTable(std::ostream& out, const Style& style, int n_steps, double time_step, double t_start);
  void print_header() const;
  void print_row(const Averages& a);
  void print_note(const std::string& text, double time);
  // Progress bar, step, speed and estimated remaining time; only with colour (not into logs).
  void print_status(double time);
  // Removes the progress line and prints the propagation time and the time per step.
  void finish();

 private:
  void clear_status();
  [[nodiscard]] double elapsed_seconds() const;

  std::ostream& out_;
  const Style& style_;
  int n_steps_;
  double time_step_;
  double t_start_;
  std::chrono::steady_clock::time_point start_;
  bool status_shown_ = false;
};

// Full-precision averages data file (gate 3).
void write_averages_file(const std::filesystem::path& path, const std::vector<Averages>& averages);

}  // namespace spodvr

#pragma once

#include <filesystem>
#include <ostream>
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

// The averages table, printed row by row while the propagation runs (Fortran: averages.out).
class AveragesTable {
 public:
  AveragesTable(std::ostream& out, const Style& style) : out_(out), style_(style) {}
  void print_header() const;
  void print_row(const Averages& a) const;
  void print_note(const std::string& text) const;

 private:
  std::ostream& out_;
  const Style& style_;
};

// Full-precision averages data file (gate 3).
void write_averages_file(const std::filesystem::path& path, const std::vector<Averages>& averages);

}  // namespace spodvr

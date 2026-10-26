// Shared code of the propagation regression tests (gate 3): parsing the Fortran stdout of a
// full run (tools/fortran_reference/run_propagation_reference.sh) and running the C++
// propagation of the same case.
#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/propagation.hpp"

namespace spodvr::testing {

// The Fortran prints these values with 8 decimals (E15.8): compare to that precision.
constexpr double printed_tolerance = 1e-8;  // relative

struct FortranPropagationReference {
  double norm = 0.0;           // '* Norm initial wave packet (Z,r) ='
  double k_at_zero = 0.0;      // '- at k = 0'
  double k_at_edge = 0.0;      // '- at the end/beginning of the grid'
  double k_maximum = 0.0;      // '- maximum value ='
  int i_low = 0;               // 1-based
  int i_high = 0;              // 1-based
  double norm_inside = 0.0;    // '* Norm of the initial wave packet between iLow and iHigh'
  std::vector<Averages> averages;
};

inline double parse_fortran_real(std::string token) {
  for (char& c : token)
    if (c == 'D' || c == 'd') c = 'E';
  return std::stod(token);
}

// Last whitespace-separated token of a line.
inline std::string last_token(const std::string& line) {
  std::istringstream in(line);
  std::string token, last;
  while (in >> token) last = token;
  return last;
}

// Five numbers (time, norm, <Z>, <r>, rest norm), or nothing if the line is not an averages row.
inline std::optional<Averages> parse_averages_row(const std::string& line) {
  std::istringstream in(line);
  std::vector<std::string> tokens;
  std::string token;
  while (in >> token) tokens.push_back(token);
  if (tokens.size() != 5) return std::nullopt;
  try {
    return Averages{parse_fortran_real(tokens[0]), parse_fortran_real(tokens[1]), parse_fortran_real(tokens[2]),
                    parse_fortran_real(tokens[3]), parse_fortran_real(tokens[4])};
  } catch (const std::exception&) {
    return std::nullopt;
  }
}

inline std::optional<FortranPropagationReference> read_fortran_reference(const std::filesystem::path& stdout_file) {
  std::ifstream file(stdout_file);
  if (!file) return std::nullopt;
  FortranPropagationReference ref;
  bool in_averages = false;
  std::string line;
  while (std::getline(file, line)) {
    const auto has = [&](const char* text) { return line.find(text) != std::string::npos; };
    if (has("* Norm initial wave packet (Z,r) =")) ref.norm = parse_fortran_real(last_token(line));
    else if (has("- at k = 0")) ref.k_at_zero = parse_fortran_real(last_token(line));
    else if (has("- at the end/beginning of the grid")) ref.k_at_edge = parse_fortran_real(last_token(line));
    else if (has("- maximum value =")) ref.k_maximum = parse_fortran_real(last_token(line));
    else if (has("iLow =")) ref.i_low = std::stoi(last_token(line));
    else if (has("iHigh =")) ref.i_high = std::stoi(last_token(line));
    else if (has("* Norm of the initial wave packet between iLow and iHigh")) ref.norm_inside = parse_fortran_real(last_token(line));
    else if (has("time") && has("norm") && has("rest norm")) in_averages = true;
    else if (has("CPU TIMES")) in_averages = false;
    else if (in_averages)
      if (const auto row = parse_averages_row(line)) ref.averages.push_back(*row);
  }
  return ref;
}

inline double relative_difference(double a, double b) {
  if (a == b) return 0.0;
  return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
}

// The C++ run of a reference case, up to tEnd.
struct PropagationRun {
  InitialWavePacketInfo initial;
  std::vector<Averages> averages;
};

inline PropagationRun run_cpp_propagation(const std::filesystem::path& input_json) {
  std::ifstream file(input_json);
  const Config config = build_config(parse_input(nlohmann::json::parse(file)));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear,
                            config.input.j_in, config.input.mj_in);
  PotentialGrid potential = compute_potential_grid(*make_potential(config.input.system_id), config, grids, angular);
  const FftwFft fft;
  const GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver());
  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft);
  PropagationRun run;
  run.initial = propagation.initial_wave_packet();
  propagation.run([&](const Averages& a) { run.averages.push_back(a); }, [](double) {});
  return run;
}

// Gate 3: initial wave packet block and the averages table against the Fortran stdout.
inline void compare_with_fortran(const PropagationRun& cpp, const FortranPropagationReference& ref) {
  EXPECT_LE(relative_difference(cpp.initial.norm, ref.norm), printed_tolerance);
  EXPECT_LE(relative_difference(cpp.initial.k_space_at_zero, ref.k_at_zero), printed_tolerance);
  EXPECT_LE(relative_difference(cpp.initial.k_space_maximum, ref.k_maximum), printed_tolerance);
  EXPECT_NEAR(cpp.initial.k_space_at_edge, ref.k_at_edge, 1e-15);  // ~0: only absolute makes sense
  EXPECT_EQ(cpp.initial.i_low + 1, ref.i_low);
  EXPECT_EQ(cpp.initial.i_high + 1, ref.i_high);
  EXPECT_LE(relative_difference(cpp.initial.norm_inside, ref.norm_inside), printed_tolerance);

  ASSERT_EQ(cpp.averages.size(), ref.averages.size());
  for (std::size_t i = 0; i < ref.averages.size(); ++i) {
    const Averages& c = cpp.averages[i];
    const Averages& f = ref.averages[i];
    EXPECT_LE(relative_difference(c.time, f.time), printed_tolerance) << "row " << i;
    EXPECT_LE(relative_difference(c.norm, f.norm), printed_tolerance) << "t = " << f.time;
    EXPECT_LE(relative_difference(c.z_mean, f.z_mean), printed_tolerance) << "t = " << f.time;
    EXPECT_LE(relative_difference(c.r_mean, f.r_mean), printed_tolerance) << "t = " << f.time;
    EXPECT_LE(relative_difference(c.rest_norm, f.rest_norm), printed_tolerance) << "t = " << f.time;
  }
}

}  // namespace spodvr::testing

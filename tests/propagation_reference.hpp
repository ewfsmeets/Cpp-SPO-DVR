// Shared code of the propagation regression tests (gate 3): parsing the Fortran stdout of a
// full run (tools/fortran_reference/run_propagation_reference.sh) and running the C++
// propagation of the same case.
#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "comparison.hpp"
#include "core/unit_conversion.hpp"
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
#include "reference_case_input.hpp"

namespace spodvr::testing {

#define EXPECT_PRINTED_EQ(cpp_value, fortran_value) \
  EXPECT_LE(std::abs((cpp_value) - (fortran_value)), printed_half_unit(fortran_value))

// The momentum-space values of the initial wave packet (|psi(k)| at k = 0, at the grid edge, its
// maximum): |c - f| <= max(half a printed unit of f, 1e-14). At k = 0 or at the edge the value can be
// FFT rounding noise, depending on the energy range and dZ (scm_cu111_srp48_small: 2.7e-15 at k = 0),
// where the printed digits of the noise mean nothing: differences below 1e-14 absolute are accepted
// (docs/validation.md, gate 3). NaN on either side fails.
inline constexpr double momentum_noise_floor = 1e-14;
inline double momentum_tolerance(double fortran_value) {
  return nan_aware_max(printed_half_unit(fortran_value), momentum_noise_floor);
}
#define EXPECT_MOMENTUM_EQ(cpp_value, fortran_value) \
  EXPECT_LE(std::abs((cpp_value) - (fortran_value)), momentum_tolerance(fortran_value))

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
    if (c == 'D' || c == 'd' || c == 'e') c = 'E';
  // gfortran drops the letter of a three-digit exponent (0.100000-119): put it back, else std::stod stops there.
  if (const auto sign = token.find_first_of("+-", 1); token.find('E') == std::string::npos && sign != std::string::npos)
    token.insert(sign, "E");
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

// Averages rows of a Fortran averages.out (only the table, no other output).
inline std::vector<Averages> read_fortran_averages(const std::filesystem::path& averages_file) {
  std::ifstream file(averages_file);
  std::vector<Averages> rows;
  std::string line;
  while (std::getline(file, line))
    if (const auto row = parse_averages_row(line)) rows.push_back(*row);
  return rows;
}

// Initial wave-packet block from the Fortran stdout, the averages from averages.out
// (the stdout also has five-column probability tables after the propagation).
inline std::optional<FortranPropagationReference> read_fortran_reference(const std::filesystem::path& stdout_file,
                                                                         const std::filesystem::path& averages_file) {
  std::ifstream file(stdout_file);
  if (!file || !std::filesystem::exists(averages_file)) return std::nullopt;
  FortranPropagationReference ref;
  std::string line;
  while (std::getline(file, line)) {
    const auto has = [&](std::string_view text) { return line.find(text) != std::string::npos; };
    if (has("* Norm initial wave packet (Z,r) ="))
      ref.norm = parse_fortran_real(last_token(line));
    else if (has("- at k = 0"))
      ref.k_at_zero = parse_fortran_real(last_token(line));
    else if (has("- at the end/beginning of the grid"))
      ref.k_at_edge = parse_fortran_real(last_token(line));
    else if (has("- maximum value ="))
      ref.k_maximum = parse_fortran_real(last_token(line));
    else if (has("iLow ="))
      ref.i_low = std::stoi(last_token(line));
    else if (has("iHigh ="))
      ref.i_high = std::stoi(last_token(line));
    else if (has("* Norm of the initial wave packet between iLow and iHigh"))
      ref.norm_inside = parse_fortran_real(last_token(line));
  }
  ref.averages = read_fortran_averages(averages_file);
  return ref;
}

// The labels of one analysed rotation-diffraction state: (j, mJ) and the diffraction channel's order and (n, m).
struct StateLabel {
  int j = 0, mj = 0, order = 0, n = 0, m = 0;
};

// The C++ run of a reference case, up to tEnd.
struct PropagationRun {
  InitialWavePacketInfo initial;
  std::vector<Averages> averages;
  std::vector<std::vector<std::complex<double>>> bk;  // BK coefficients per analysis step
  ScatteringResults scattering;
  std::vector<std::pair<int, int>> rotation_labels;   // (J, mJ) of the analysed rotational states
  int initial_state = 0;                              // index is of the initial state in the analysis basis
  int n_states = 0;                                   // analysed states (the s_matrix layout is + n_states * v)
  std::optional<double> specular_line;                // Z of the fixed specular analysis line, if it became fixed
  std::vector<StateLabel> states;                     // per analysed state is
  int n_v = 0;                                        // analysed vibrational levels
  std::vector<std::vector<double>> level_energy_ev;   // E(v, j) of the gas phase (eV), v < n_v, j <= JMax
  ScatteringResults scattering_zana;                  // the same S-matrix with every phase referred to Zana
};

// The C++ run of a configuration with the PES pes in place of its own (a negative test's modified
// potential); backend: the propagation backend (CPU by default; the GPU tests pass the CUDA backend).
inline PropagationRun run_cpp_propagation(const Config& config, const Potential& pes,
                                          const BackendFactory& backend = make_cpu_backend) {
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  PotentialGrid potential = compute_potential_grid(pes, config, grids, angular);
  const FftwFft fft;
  const GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver());
  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft, backend);
  PropagationRun run;
  run.initial = propagation.initial_wave_packet();
  PropagationObserver observer;
  observer.on_averages = [&](const Averages& a) { run.averages.push_back(a); };
  observer.on_analysis_step = [&](double, const BkCoefficients& bk) { run.bk.push_back(bk.values); };
  propagation.run(observer);
  run.scattering = propagation.scattering_results();
  run.scattering_zana = propagation.scattering_results(std::nullopt);
  run.specular_line = propagation.analysis_surface().fixed_line();
  const AnalysisBasis& basis = propagation.analysis_basis();
  run.initial_state = basis.initial_state();
  run.n_states = basis.n_states();
  run.n_v = basis.n_v();
  for (int is = 0; is < basis.n_states(); ++is) {
    const DiffractionChannel& c = basis.channel(is);
    run.states.push_back({basis.j(is), basis.mj(is), c.order, c.n, c.m});
  }
  const double chaev = units::legacy.energy_factor();
  for (int v = 0; v < basis.n_v(); ++v) {
    run.level_energy_ev.emplace_back();
    for (int j = 0; j <= config.input.j_max; ++j) run.level_energy_ev.back().push_back(gas.states.energy(v, j) * chaev);
  }
  for (int irotana = 0; irotana < basis.n_rotation(); ++irotana)
    run.rotation_labels.emplace_back(basis.j(irotana * basis.n_diffraction()),
                                     basis.mj(irotana * basis.n_diffraction()));
  return run;
}

// The C++ run of the JSON input (a reference case's: case_config, reference_case_input.hpp), with
// its own PES.
inline PropagationRun run_cpp_propagation(const nlohmann::json& input,
                                          const BackendFactory& backend = make_cpu_backend) {
  const CaseSetup setup = case_setup(input);
  return run_cpp_propagation(setup.config, *setup.pes, backend);
}

// The C++ run of the JSON input file input_json.
inline PropagationRun run_cpp_propagation(const std::filesystem::path& input_json,
                                          const BackendFactory& backend = make_cpu_backend) {
  std::ifstream file(input_json);
  return run_cpp_propagation(nlohmann::json::parse(file), backend);
}

// Largest |C++ - Fortran| of the averages (norm, <Z>, <r>, rest norm) in units of the printed half
// unit, over the rows both runs have (the callers check the row counts). NaN if a value is NaN.
inline double averages_deviation(const PropagationRun& cpp, const FortranPropagationReference& ref) {
  double worst = 0.0;
  for (std::size_t i = 0; i < std::min(cpp.averages.size(), ref.averages.size()); ++i) {
    const Averages& c = cpp.averages[i];
    const Averages& f = ref.averages[i];
    for (const auto& [a, b] :
         {std::pair{c.norm, f.norm}, {c.z_mean, f.z_mean}, {c.r_mean, f.r_mean}, {c.rest_norm, f.rest_norm}})
      worst = nan_aware_max(worst, std::abs(a - b) / printed_half_unit(b));
  }
  return worst;
}

// Gate 3: initial wave packet block and the averages table against the Fortran output.
inline void compare_with_fortran(const PropagationRun& cpp, const FortranPropagationReference& ref) {
  EXPECT_PRINTED_EQ(cpp.initial.norm, ref.norm);
  EXPECT_MOMENTUM_EQ(cpp.initial.k_space_at_zero, ref.k_at_zero);
  EXPECT_MOMENTUM_EQ(cpp.initial.k_space_maximum, ref.k_maximum);
  EXPECT_MOMENTUM_EQ(cpp.initial.k_space_at_edge, ref.k_at_edge);
  EXPECT_EQ(cpp.initial.i_low + 1, ref.i_low);
  EXPECT_EQ(cpp.initial.i_high + 1, ref.i_high);
  EXPECT_PRINTED_EQ(cpp.initial.norm_inside, ref.norm_inside);

  ASSERT_EQ(cpp.averages.size(), ref.averages.size());
  for (std::size_t i = 0; i < ref.averages.size(); ++i) {
    const Averages& c = cpp.averages[i];
    const Averages& f = ref.averages[i];
    EXPECT_PRINTED_EQ(c.time, f.time) << "row " << i;
    EXPECT_PRINTED_EQ(c.norm, f.norm) << "t = " << f.time;
    EXPECT_PRINTED_EQ(c.z_mean, f.z_mean) << "t = " << f.time;
    EXPECT_PRINTED_EQ(c.r_mean, f.r_mean) << "t = " << f.time;
    EXPECT_PRINTED_EQ(c.rest_norm, f.rest_norm) << "t = " << f.time;
  }
  std::cout << "  largest deviation from the Fortran averages: " << averages_deviation(cpp, ref)
            << " x (half printed unit)\n";
}

}  // namespace spodvr::testing

#pragma once
// Gate 6 (docs/validation.md): the flux analysis of a reference case with doFluxAnalysis against the
// Fortran's flux files of the same case (fortran_flux_*, tools/fortran_reference/
// run_propagation_reference.sh), in three comparisons:
//   (a) phi and phi', recomputed here from the run's flux surfaces with the Fortran's eviroi, and the
//       run's own (C++ eviroi), element by element within tau_a B_n, where B_n = |c_n| sum_t max |psi_c(t)|
//       bounds every |phi| of energy n (B'_n with d psi / dr for phi');
//   (b) the run's own P(E) and P(x, y; E) against flux_reaction.dat and flux_density.dat, within
//       max(half a printed unit, tau_b prefac sum |phi| |phi'|);
//   (c) the Fortran's phi and phi' through the C++ flux_results: the Fortran's E15.8 text exactly (the
//       flux analysis is compiled with strict floating point in every build).
// Negative controls on the same data: realistic porting mistakes and errors just above a tolerance fail
// (a) to (c). Used by run_gate (regression_gate.hpp), so on the CPU and the GPU. Include it only in test
// programs linked with spodvr_strict_fp (as flux_propagation_run.hpp, whose sums it uses).

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "comparison.hpp"
#include "core/array_index.hpp"
#include "core/constants.hpp"
#include "flux_propagation_run.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/fortran_files.hpp"
#include "io/fortran_unformatted.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation_reference.hpp"

namespace spodvr::testing {

// (a): |phi_C++ - phi_Fortran| <= tau_a B_n. B_n bounds every |phi| of energy n over the whole run, so
// tau_a is an accuracy of the time integral, not of each (possibly strongly cancelling) element.
// Measured on cu111_dailight_short (2026-11-13), the larger of phi and phi', in units of B_n: with the
// Fortran's eviroi 2.9e-13 (CPU), 2.0e-13 (GPU), 6.0e-13 (build-fast and build-cuda-fast, whose psi
// carries fast-math and FMA rounding); with the C++ eviroi (1.1e-15 to 1.2e-15 hartree off) 1.8e-12
// to 1.9e-12 in all four builds. One tolerance for every build: margins 17 to 35 and 5.
inline constexpr double flux_coefficient_tolerance = 1e-11;
// (b): P may differ by this share of prefac sum |phi| |phi'|, the size of its terms before they cancel.
// On cu111_dailight_short the printed half unit is the larger term at 189 of 195 values (worst 0.977
// of the tolerance, CPU and GPU, 2026-11-13): this term matters only for values near zero.
inline constexpr double flux_probability_tolerance = 1e-10;

// The four Fortran flux files of a reference case (tools/fortran_reference/run_propagation_reference.sh).
inline const std::vector<std::string> fortran_flux_files = {"fortran_flux_coeffs.bin", "fortran_flux_coeffs_deriv.bin",
                                                            "fortran_flux_reaction.dat", "fortran_flux_density.dat"};

// flux_density.dat: 4 header lines (energies; skewing angle and point group; nx, dx; ny, dy), then one
// line per (E, x, y) with iy fastest. Fortran: FluxAnalyzor.f90 computeFluxProbability
inline constexpr std::size_t fortran_density_header_lines = 4;

// The line of flux_density.dat that holds P(x(ix), y(iy); E_n) (0-based).
inline std::size_t fortran_density_line(int n, int ix, int iy, const FluxGrid& grid, int n_energies) {
  const Index3D rows{std::size_t(grid.ny), std::size_t(grid.nx), std::size_t(n_energies)};
  return fortran_density_header_lines + rows(iy, ix, n);
}

// eviroi, the energy of the initial rovibrational state (hartree): record 2 of a Fortran BK file
// holds deltat, rkz, z0, alfwv, ZAnalysis, eviroi. Fortran: AsympInfo.f90 ReadOrWriteAsympInfo
inline double read_fortran_initial_state_energy(const std::filesystem::path& bk_file) {
  FortranRecordReader file(bk_file);
  (void)file.read("record 1 (sizes)");
  FortranRecordParser record = file.read("record 2 (deltat, rkz, z0, alfwv, ZAnalysis, eviroi)");
  for (int k = 0; k < 5; ++k) (void)record.read<double>();
  const double eviroi = record.read<double>();
  record.expect_end();
  return eviroi;
}

// The whitespace-separated tokens of each line of a text file.
inline std::vector<std::vector<std::string>> read_tokens(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::vector<std::vector<std::string>> lines;
  for (std::string line; std::getline(file, line);) {
    std::istringstream words(line);
    std::vector<std::string> tokens;
    for (std::string word; words >> word;) tokens.push_back(word);
    lines.push_back(tokens);
  }
  return lines;
}

// x as gfortran writes it with E15.8, without the leading blank of a positive number: "0.", the 8
// correctly rounded digits, 'E' and a signed 2-digit exponent (3 digits, without the 'E', beyond 99).
inline std::string fortran_e15_8(double x) {
  const std::string sign = std::signbit(x) ? "-" : "";
  if (x == 0.0) return sign + "0.00000000E+00";
  const std::string scientific = std::format("{:.7e}", std::abs(x));  // d.ddddddde±xx
  const std::string digits = scientific.substr(0, 1) + scientific.substr(2, 7);
  const int exponent = std::stoi(scientific.substr(scientific.find('e') + 1)) + 1;
  const char exponent_sign = exponent < 0 ? '-' : '+';
  if (std::abs(exponent) > 99) return std::format("{}0.{}{}{:03}", sign, digits, exponent_sign, std::abs(exponent));
  return std::format("{}0.{}E{}{:02}", sign, digits, exponent_sign, std::abs(exponent));
}

// The largest |a - b| / bound(n) of two coefficient arrays (layout (block, iz, n)), the number of
// elements where it exceeds `tolerance` (a NaN counts), and the number with |b| above it (those
// that a wrong value cannot pass by being close to zero), out of `compared`.
struct CoefficientDeviation {
  double worst = 0.0;
  int above = 0;
  int significant = 0;
  int compared = 0;
};

inline CoefficientDeviation coefficient_deviation(std::span<const std::complex<double>> a,
                                                  std::span<const std::complex<double>> b, const Index3D& index,
                                                  const std::vector<double>& bound, double tolerance) {
  CoefficientDeviation d;
  for (std::size_t n = 0; n < index.n2; ++n)
    for (std::size_t iz = 0; iz < index.n1; ++iz)
      for (std::size_t block = 0; block < index.n0; ++block) {
        const std::size_t k = index(block, iz, n);
        const double deviation = std::abs(a[k] - b[k]) / bound[n];
        d.worst = nan_aware_max(d.worst, deviation);
        if (!(deviation <= tolerance)) ++d.above;
        if (std::abs(b[k]) / bound[n] > tolerance) ++d.significant;
        ++d.compared;
      }
  return d;
}

// The largest modulus of a set of values (NaN if one is NaN).
inline double largest_modulus(std::span<const std::complex<double>> values) {
  double largest = 0.0;
  for (const std::complex<double>& v : values) largest = nan_aware_max(largest, std::abs(v));
  return largest;
}

// prefac_n sum_iz sum_irot |phi| |phi'| at every (ix, iy, n), ix fastest (FluxResults::density_index):
// the size of the terms of P(x, y; E_n) before they cancel. prefac_n = ((2 pi M) / (k_n mu)) dz.
inline std::vector<double> flux_term_magnitudes(const FluxCoefficients& c, const FluxEnergies& energies,
                                                const FluxGrid& grid) {
  const Index3D density_index{std::size_t(grid.nx), std::size_t(grid.ny), std::size_t(c.n_energies())};
  std::vector<double> sums(density_index.size(), 0.0);
  for (int n = 0; n < c.n_energies(); ++n) {
    const double prefac =
        ((Constants::two_pi * grid.total_mass) / (energies.momentum[n] * grid.reduced_mass)) * grid.dz;
    for (int iz = 0; iz < c.n_planes(); ++iz)
      for (int block = 0; block < c.n_blocks(); ++block) {
        const auto [ix, iy] = flux_block_xy(block, grid);
        const std::size_t k = c.index(block, iz, n);
        sums[density_index(ix, iy, n)] += prefac * std::abs(c.phi[k]) * std::abs(c.phi_derivative[k]);
      }
  }
  return sums;
}

// (b): the largest |C++ - Fortran| / tolerance of P(E) and P(x, y; E), with the tolerance
// max(half a printed unit, tau_b x term magnitude) per value; the number of values above 1, of values
// where the printed unit is the larger term, and whether it is at the worst value.
struct ProbabilityDeviation {
  double worst = 0.0;
  bool printed_unit_binds_at_worst = false;
  int above = 0;
  int compared = 0;
  int printed_unit_binds = 0;
};

// Adds one value of (b): the C++ value, the Fortran's E15.8 text and the size of its terms.
inline void add_probability(ProbabilityDeviation& d, double cpp, const std::string& fortran_text, double magnitude) {
  const double fortran = parse_fortran_real(fortran_text);
  const bool printed_unit_binds = printed_half_unit(fortran) >= flux_probability_tolerance * magnitude;
  const double tolerance = std::max(printed_half_unit(fortran), flux_probability_tolerance * magnitude);
  const double deviation = std::abs(cpp - fortran) / tolerance;
  if (!(deviation <= d.worst)) d.printed_unit_binds_at_worst = printed_unit_binds;  // a new worst (or NaN)
  d.worst = nan_aware_max(d.worst, deviation);
  if (!(deviation <= 1.0)) ++d.above;
  if (printed_unit_binds) ++d.printed_unit_binds;
  ++d.compared;
}

// The Fortran files of gate 6 of one case, read once.
struct FortranFluxFiles {
  double initial_state_energy = 0.0;  // eviroi (hartree)
  std::vector<std::vector<std::string>> reaction;  // flux_reaction.dat: E (eV), P per line
  std::vector<std::vector<std::string>> density;   // flux_density.dat: header, then E, x, y, P per line
};

// Gate 6 for one reference case: the case's set-up, the Fortran files and the sums this test makes
// from the run's flux surfaces (with the Fortran's eviroi, so that (a) compares the time integrals of
// the two wave functions, not two rounded eviroi). Built by for_case only for a case with
// doFluxAnalysis and the Fortran flux files.
class FluxGate {
 public:
  // The gate of the case in `dir`, or nothing for a case without doFluxAnalysis. A case with it but
  // without its four Fortran flux files fails (a missing reference must not look green).
  static std::optional<FluxGate> for_case(const std::filesystem::path& dir) {
    std::ifstream file(dir / "input.json");
    const Config config = build_config(parse_input(nlohmann::json::parse(file)));
    std::vector<std::string> missing;
    for (const std::string& name : fortran_flux_files)
      if (!std::filesystem::exists(dir / name)) missing.push_back(name);
    if (!config.input.do_flux_analysis) {
      EXPECT_EQ(missing.size(), fortran_flux_files.size())
          << "Fortran flux files in " << dir << ", but the input has no doFluxAnalysis";
      return std::nullopt;
    }
    if (!missing.empty()) {
      std::string names;
      for (const std::string& name : missing) names += " " + name;
      ADD_FAILURE() << "gate 6: doFluxAnalysis in " << dir << ", but no" << names
                    << "; run tools/fortran_reference/run_propagation_reference.sh " << dir.filename().string();
      return std::nullopt;
    }
    return FluxGate(dir, config);
  }

  // Adds each analysis step's flux surface to this test's sums: pass to run_cpp_propagation. It refers
  // to this gate, which must stay in place (not be moved) until the run has finished.
  [[nodiscard]] FluxSurfaceObserver observer() {
    return [this](double time, const FluxSurfaceWave& surface) { add_step(time, surface); };
  }

  // Gate 6 (a), (b), (c) and its negative controls, after the run.
  void compare(const PropagationRun& run) const {
    ASSERT_TRUE(run.flux && run.flux_coefficients) << "the C++ run has no flux analysis";
    ASSERT_EQ(std::size_t(steps_), run.bk.size()) << "one flux surface per analysis step";
    compare_coefficients(run);
    compare_probabilities(run);
    compare_fortran_coefficients_text();
    expect_mistakes_fail(run);
  }

 private:
  // Reads the case's set-up and the Fortran files; the sums start at zero.
  FluxGate(const std::filesystem::path& dir, const Config& config)
      : config_(config),
        grids_(make_spatial_grids(config_)),
        n_rot_(RotationalBasis(make_angular_grid(config_), config_.input.j_in, config_.input.mj_in).size()),
        grid_(make_flux_grid(config_, grids_, n_rot_)),
        files_{read_fortran_initial_state_energy(dir / "fortran_fort21.bin"),
               read_tokens(dir / "fortran_flux_reaction.dat"), read_tokens(dir / "fortran_flux_density.dat")},
        energies_(make_flux_energies(config_, parallel_kinetic_energies(config_)[0], files_.initial_state_energy)),
        phase_slip_energies_(with_phase_slip(energies_)),
        fortran_(n_rot_ * grid_.nx * grid_.ny, grids_.nz_potential, energies_.size()),
        sums_(fortran_.n_blocks(), fortran_.n_planes(), fortran_.n_energies()),
        phase_slip_sums_(fortran_.n_blocks(), fortran_.n_planes(), fortran_.n_energies()) {
    read_fortran_flux_coefficients(dir / "fortran_flux_coeffs.bin", dir / "fortran_flux_coeffs_deriv.bin", fortran_);
  }

  // The energies with the sign of E_tot,n flipped: c_n with exp(-i E t), a sign slip in the phase.
  static FluxEnergies with_phase_slip(FluxEnergies energies) {
    for (double& e : energies.total) e = -e;
    return energies;
  }

  // S and S' (the bounds' sums over time), and phi, phi' with c_n from the Fortran's eviroi and with the
  // phase slip, in the Fortran's order (add_flux_surface_by_hand).
  void add_step(double time, const FluxSurfaceWave& surface) {
    psi_sum_ += largest_modulus(surface.psi);
    derivative_sum_ += largest_modulus(surface.derivative);
    add_flux_surface_by_hand(sums_, energies_, time, config_.derived.analysis_time_step, surface);
    add_flux_surface_by_hand(phase_slip_sums_, phase_slip_energies_, time, config_.derived.analysis_time_step, surface);
    ++steps_;
  }

  // B_n = |c_n| S (phi) or |c_n| S' (phi'); |c_n| = |A_n| dt at every step.
  [[nodiscard]] std::vector<double> bounds(double time_sum) const {
    std::vector<double> b;
    for (int n = 0; n < energies_.size(); ++n)
      b.push_back(std::abs(energies_.transform_factor[n]) * config_.derived.analysis_time_step * time_sum);
    return b;
  }

  // (a), for this test's sums (the Fortran's eviroi) and the run's own coefficients (the C++ eviroi).
  void compare_coefficients(const PropagationRun& run) const {
    const double tau = flux_coefficient_tolerance;
    const std::vector<double> b = bounds(psi_sum_);
    const std::vector<double> b_derivative = bounds(derivative_sum_);
    const CoefficientDeviation phi = coefficient_deviation(sums_.phi, fortran_.phi, sums_.index, b, tau);
    const CoefficientDeviation derivative =
        coefficient_deviation(sums_.phi_derivative, fortran_.phi_derivative, sums_.index, b_derivative, tau);
    const CoefficientDeviation own =
        coefficient_deviation(run.flux_coefficients->phi, fortran_.phi, sums_.index, b, tau);
    const CoefficientDeviation own_derivative = coefficient_deviation(
        run.flux_coefficients->phi_derivative, fortran_.phi_derivative, sums_.index, b_derivative, tau);
    EXPECT_LE(phi.worst, tau);
    EXPECT_LE(derivative.worst, tau);
    EXPECT_LE(own.worst, tau);
    EXPECT_LE(own_derivative.worst, tau);
    std::cout << "  gate 6 (a): " << phi.compared << " values each of phi and phi', worst |C++ - Fortran| " << phi.worst
              << " and " << derivative.worst << " of B_n (tolerance " << tau
              << "); the run's own, with the C++ eviroi (" << run.initial_state_energy - files_.initial_state_energy
              << " hartree from the Fortran's): " << own.worst << " and " << own_derivative.worst << '\n';
  }

  // (b) for the results `results` and the coefficients `c` they come from.
  [[nodiscard]] ProbabilityDeviation probability_deviation(const FluxResults& results,
                                                           const FluxCoefficients& c) const {
    const std::vector<double> magnitude = flux_term_magnitudes(c, energies_, grid_);
    const double area = grid_.dx * grid_.dy * std::sin(grid_.skewing_angle);
    ProbabilityDeviation d;
    for (int n = 0; n < c.n_energies(); ++n) {
      double total = 0.0;
      for (int ix = 0; ix < grid_.nx; ++ix)
        for (int iy = 0; iy < grid_.ny; ++iy) {
          const std::size_t k = results.density_index(ix, iy, n);
          total += magnitude[k];
          // at(): a short file throws (a test failure) instead of reading past its end.
          const std::vector<std::string>& line =
              files_.density.at(fortran_density_line(n, ix, iy, grid_, c.n_energies()));
          add_probability(d, results.density[k], line.at(3), magnitude[k]);
        }
      add_probability(d, results.probability[n], files_.reaction.at(n).at(1), total * area);
    }
    return d;
  }

  // (b) for the run: its P(E) and P(x, y; E) against the Fortran's files.
  void compare_probabilities(const PropagationRun& run) const {
    ASSERT_EQ(files_.reaction.size(), std::size_t(energies_.size()));
    // One past the last line: the line of the energy after the last one.
    ASSERT_EQ(files_.density.size(), fortran_density_line(energies_.size(), 0, 0, grid_, energies_.size()));
    const ProbabilityDeviation d = probability_deviation(*run.flux, *run.flux_coefficients);
    EXPECT_LE(d.worst, 1.0);
    std::cout << "  gate 6 (b): " << d.compared << " values of P(E) and P(x, y; E), worst " << d.worst
              << " x tolerance max(half a printed unit, " << flux_probability_tolerance
              << " prefac sum |phi| |phi'|), where the "
              << (d.printed_unit_binds_at_worst ? "printed unit" : "second term")
              << " is the larger; the printed unit is the larger term at " << d.printed_unit_binds << '\n';
  }

  // The tokens of flux_reaction.dat that `results` would give.
  [[nodiscard]] std::vector<std::vector<std::string>> reaction_text(const FluxResults& results) const {
    std::vector<std::vector<std::string>> lines;
    for (int n = 0; n < energies_.size(); ++n)
      lines.push_back({fortran_e15_8(results.energy[n]), fortran_e15_8(results.probability[n])});
    return lines;
  }

  // The tokens of flux_density.dat that `results` would give, its header included.
  [[nodiscard]] std::vector<std::vector<std::string>> density_text(const FluxResults& results) const {
    std::vector<std::vector<std::string>> lines = {
        {std::to_string(energies_.size())},
        {fortran_e15_8(grid_.skewing_angle), config_.derived.system.point_group},
        {std::to_string(grid_.nx), fortran_e15_8(grid_.dx)},
        {std::to_string(grid_.ny), fortran_e15_8(grid_.dy)}};
    for (int n = 0; n < energies_.size(); ++n)
      for (int ix = 0; ix < grid_.nx; ++ix)
        for (int iy = 0; iy < grid_.ny; ++iy)
          lines.push_back({fortran_e15_8(results.energy[n]), fortran_e15_8(grids_.x[ix]), fortran_e15_8(grids_.y[iy]),
                           fortran_e15_8(results.density[results.density_index(ix, iy, n)])});
    return lines;
  }

  // Tokens of the Fortran files that the C++ text of `results` does not reproduce (a missing or
  // extra line counts all its tokens).
  [[nodiscard]] int text_differences(const FluxResults& results) const {
    int differences = 0;
    auto count = [&](const std::vector<std::vector<std::string>>& cpp, const std::vector<std::vector<std::string>>& f) {
      for (std::size_t line = 0; line < std::max(cpp.size(), f.size()); ++line) {
        const std::vector<std::string> none;
        const std::vector<std::string>& a = line < cpp.size() ? cpp[line] : none;
        const std::vector<std::string>& b = line < f.size() ? f[line] : none;
        for (std::size_t k = 0; k < std::max(a.size(), b.size()); ++k)
          if (k >= a.size() || k >= b.size() || a[k] != b[k]) ++differences;
      }
    };
    count(reaction_text(results), files_.reaction);
    count(density_text(results), files_.density);
    return differences;
  }

  // (c): the Fortran's phi and phi' through the C++ flux_results give the Fortran's text.
  void compare_fortran_coefficients_text() const {
    const int differences = text_differences(flux_results(fortran_, energies_, grid_));
    EXPECT_EQ(differences, 0);
    std::size_t tokens = 0;
    for (const auto& line : files_.reaction) tokens += line.size();
    for (const auto& line : files_.density) tokens += line.size();
    std::cout << "  gate 6 (c): the Fortran's phi and phi' through flux_results: " << differences << " of " << tokens
              << " values of the two files differ in their E15.8 text\n";
  }

  // The element (block, iz, n) of the Fortran's phi with the largest modulus, and its energy n.
  struct LargestElement {
    std::size_t k = 0;
    int n = 0;
  };
  [[nodiscard]] LargestElement largest_fortran_phi() const {
    LargestElement largest;
    for (int n = 0; n < fortran_.n_energies(); ++n)
      for (int iz = 0; iz < fortran_.n_planes(); ++iz)
        for (int block = 0; block < fortran_.n_blocks(); ++block) {
          const std::size_t k = fortran_.index(block, iz, n);
          if (std::abs(fortran_.phi[k]) > std::abs(fortran_.phi[largest.k])) largest = {k, n};
        }
    return largest;
  }

  // Negative controls: on the data of this run, realistic porting mistakes and errors just above a
  // tolerance fail (a), (b) or (c), and by how much. Not here: A_n as k conj(d) / |d|^2 instead of the
  // complex division moves c_n by about an ulp, far below (a); FluxTransformFactor.IsSmithsComplexDivision
  // catches it.
  void expect_mistakes_fail(const PropagationRun& run) const {
    const double tau = flux_coefficient_tolerance;
    const std::vector<double> b = bounds(psi_sum_);
    const std::vector<double> b_derivative = bounds(derivative_sum_);

    // c_n with exp(-i E t) (a sign slip in the phase): (a) fails at nearly every value that matters.
    const CoefficientDeviation slip = coefficient_deviation(phase_slip_sums_.phi, fortran_.phi, sums_.index, b, tau);
    EXPECT_GT(slip.significant, 0);
    EXPECT_GE(slip.above, slip.significant * 9 / 10);
    EXPECT_GT(slip.worst, 1e6 * tau);

    // D(ir) with the distance ir - ir_c (the Fortran's jr - ir reversed): D and so phi' change sign.
    FluxCoefficients reversed = *run.flux_coefficients;
    for (std::complex<double>& v : reversed.phi_derivative) v = -v;
    const CoefficientDeviation derivative =
        coefficient_deviation(reversed.phi_derivative, fortran_.phi_derivative, sums_.index, b_derivative, tau);
    EXPECT_GT(derivative.significant, 0);
    EXPECT_GE(derivative.above, derivative.significant);
    const ProbabilityDeviation sign = probability_deviation(flux_results(reversed, energies_, grid_), reversed);
    EXPECT_EQ(sign.above, sign.compared);

    // One value of phi off by 2 tau B_n fails (a) there only; by tau B_n / 2 it passes.
    const LargestElement largest = largest_fortran_phi();
    auto shifted = [&](double share) {
      FluxCoefficients changed = sums_;
      changed.phi[largest.k] += share * tau * b[largest.n];
      return coefficient_deviation(changed.phi, fortran_.phi, sums_.index, b, tau).above;
    };
    const int shifted_above = shifted(2.0);
    EXPECT_EQ(shifted_above, 1);
    EXPECT_EQ(shifted(0.5), 0);

    // (c) with the masses swapped (M for mu in prefac): every probability's text differs.
    FluxGrid swapped = grid_;
    std::swap(swapped.total_mass, swapped.reduced_mass);
    const int masses = text_differences(flux_results(fortran_, energies_, swapped));
    const int probabilities = energies_.size() * (1 + grid_.nx * grid_.ny);
    EXPECT_EQ(masses, probabilities);

    std::cout << "  gate 6 negative controls: c_n with exp(-i E t) fails (a) at " << slip.above << " of "
              << slip.compared << " values of phi (" << slip.significant << " of them larger than tau B_n; worst "
              << slip.worst << " of B_n); D reversed fails (a) at " << derivative.above << " of " << derivative.compared
              << " values of phi' and (b) at " << sign.above << " of " << sign.compared
              << "; one phi off by 2 tau B_n: " << shifted_above << " failure; M and mu swapped: " << masses << " of "
              << probabilities << " probabilities' text differs\n";
  }

  Config config_;                    // the case's input
  SpatialGrids grids_;               // its grids
  int n_rot_;                        // rotational functions per (x, y) point
  FluxGrid grid_;                    // the grid data of flux_results
  FortranFluxFiles files_;           // eviroi and the two text files
  FluxEnergies energies_;            // with the Fortran's eviroi
  FluxEnergies phase_slip_energies_; // the same with -E_tot,n: a negative control
  FluxCoefficients fortran_;         // flux_coeffs.dat, flux_coeffs_deriv.dat
  FluxCoefficients sums_;            // (a): this test's phi and phi' from the run's flux surfaces
  FluxCoefficients phase_slip_sums_; // the same with the phase slip
  double psi_sum_ = 0.0;             // S = sum_t max |psi_c(t)|
  double derivative_sum_ = 0.0;      // S' = sum_t max |d psi_c(t) / dr|
  int steps_ = 0;                    // analysis steps added
};

}  // namespace spodvr::testing

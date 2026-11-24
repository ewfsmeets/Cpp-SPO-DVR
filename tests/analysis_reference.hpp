// Shared code of the analysis regression tests (gates 4 and 5): the analysis limits and the per-energy blocks of the
// Fortran stdout (RotDifBasis.f90, Analyzor.f90 SetupAnalysis and smatrix_combine) and the Fortran BK coefficient
// file fort.21.
#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "comparison.hpp"
#include "io/fortran_unformatted.hpp"
#include "propagation_reference.hpp"

namespace spodvr::testing {

// Half a unit in the last printed digit of a Fortran number as printed (E14.7, D14.7, F10.7, ...),
// with a 1% margin for the decimal conversion. Values below `floor` are numerical noise
// (|S|^2 of channels that are not populated: ~1e-26 in both codes) and compared absolutely.
inline double printed_text_tolerance(const std::string& text, double floor = 1e-20) {
  std::string t = text;
  for (char& c : t)
    if (c == 'D' || c == 'd') c = 'E';
  const auto e = t.find('E');
  const std::string mantissa = t.substr(0, e);
  const auto point = mantissa.find('.');
  const int digits = static_cast<int>(mantissa.size() - point - 1);
  const int exponent = e == std::string::npos ? 0 : std::stoi(t.substr(e + 1));
  return std::max(1.01 * 0.5 * std::pow(10.0, exponent - digits), floor);
}

// Half a unit in the last printed digit of a Fortran real (printed_text_tolerance, at least floor), but floor for an
// E-format zero, which is exact: a nonzero value always prints a nonzero mantissa there. An F-format zero is not
// exact (f10.6 prints 3e-7 as 0.000000) and keeps its half unit.
inline double printed_half_unit(const std::string& text, double floor) {
  const std::string t = fortran_text_with_e(text);
  const auto e = t.find('E');
  if (e != std::string::npos && t.substr(0, e).find_first_of("123456789") == std::string::npos) return floor;
  return printed_text_tolerance(t, floor);
}

inline std::vector<std::string> split(const std::string& line) {
  std::istringstream in(line);
  std::vector<std::string> tokens;
  std::string token;
  while (in >> token) tokens.push_back(token);
  return tokens;
}

inline bool all_digits(const std::string& s) {
  return !s.empty() && s.find_first_not_of("0123456789") == std::string::npos;
}

// The value of a printed integer (an optional minus sign), or nothing.
inline std::optional<int> printed_integer(const std::string& text) {
  if (!all_digits(text.starts_with('-') ? text.substr(1) : text)) return std::nullopt;
  return std::stoi(text);
}

// --- the analysis limits ----------------------------------------------------------------------------------------

// The analysis limits a Fortran stdout prints, every line found: name -> the printed values, in order. The lines:
//   * number of rotational states analysed, nana =   15      (RotDifBasis.f90)
//   * Maximum j in analysis:                 jana =    4      (Analyzor.f90 SetupAnalysis)
using FortranAnalysisLimits = std::map<std::string, std::vector<std::string>>;

// The names a stdout prints: a run that analyses (propagation, restart, fine analysis) all five; the set-up of a
// potential-only run (fortran_setup_stdout.txt) stops before SetupAnalysis, which prints jana and nviba.
inline const std::vector<std::string> limits_of_a_run = {"jana", "nviba", "nana", "nstata", "inchna"};
inline const std::vector<std::string> limits_of_a_setup = {"nana", "nstata", "inchna"};

inline FortranAnalysisLimits read_fortran_analysis_limits(const std::filesystem::path& stdout_file) {
  std::ifstream file(stdout_file);
  if (!file) throw std::runtime_error("cannot open " + stdout_file.string());
  FortranAnalysisLimits printed;
  std::string line;
  while (std::getline(file, line))
    for (const std::string& name : limits_of_a_run) {
      const auto at = line.find(name + " =");
      if (at != std::string::npos && at > 0 && (line[at - 1] == ' ' || line[at - 1] == ','))
        printed[name].push_back(last_token(line));
    }
  return printed;
}

// Lines of text, each on a line of its own (a failure message's list).
inline std::string listed(const std::vector<std::string>& lines) {
  std::string text;
  for (const std::string& line : lines) text = text + "\n  " + line;
  return text;
}

// The outcome of a comparison of the analysis limits; failures in words.
struct LimitsComparison {
  int compared = 0;
  int failed = 0;
  std::vector<std::string> failures;
};

// Every limit the stdout must print (expected_names) printed once and equal to the C++ value; any other limit not
// printed.
inline LimitsComparison compare_analysis_limits(const AnalysisLimits& cpp, const FortranAnalysisLimits& printed,
                                                const std::vector<std::string>& expected_names) {
  const std::map<std::string, int> cpp_values = {
      {"jana", cpp.jana}, {"nviba", cpp.nviba}, {"nana", cpp.nana}, {"nstata", cpp.nstata}, {"inchna", cpp.inchna}};
  LimitsComparison c;
  auto fail = [&](const std::string& what) {
    c.failed = c.failed + 1;
    c.failures.push_back(what);
  };
  for (const auto& [name, cpp_value] : cpp_values) {
    const auto found = printed.find(name);
    const std::vector<std::string> values = found == printed.end() ? std::vector<std::string>{} : found->second;
    const bool expected = std::ranges::find(expected_names, name) != expected_names.end();
    if (static_cast<int>(values.size()) != (expected ? 1 : 0))
      fail(name + " printed " + std::to_string(values.size()) + " times, expected " + (expected ? "once" : "never"));
    for (const std::string& text : values) {
      c.compared = c.compared + 1;
      if (printed_integer(text) != cpp_value) fail(name + ": C++ " + std::to_string(cpp_value) + ", Fortran " + text);
    }
  }
  return c;
}

// The gate's assertions on the analysis limits: none failed, every expected one compared.
inline void expect_analysis_limits(const LimitsComparison& c, const std::vector<std::string>& expected_names) {
  EXPECT_EQ(c.failed, 0) << "failed analysis limits:" << listed(c.failures);
  EXPECT_EQ(c.compared, static_cast<int>(expected_names.size())) << "analysis limits compared";
}

// --- the per-energy blocks (Analyzor.f90 smatrix_combine) ------------------------------------------------------

// The rows of the blocks, as printed: the 'rotation-diffraction summed probabilities' (v, P(v), av edif, av erot;
// eV), the 'diffraction summed probabilities' (v, irot (1-based), j, mJ, P(v, j, mJ)) and the 'probabilities as
// function of diffraction order' (the v, j, mJ of the header above, the order, P).
struct FortranVibrationRow {
  int v = 0;
  std::string p, e_diffraction, e_rotation;
};
struct FortranRotationRow {
  int v = 0, irot = 0, j = 0, mj = 0;
  std::string p;
};
struct FortranOrderRow {
  int v = 0, j = 0, mj = 0, order = 0;
  std::string p;
};

// One per-energy block of the Fortran stdout; problems: lines the reader could not use (an unreadable label, an
// order row before its header).
struct FortranEnergyBlock {
  std::string total_energy;      // ' total energy =  0.45837825D+00'
  std::string collision_energy;  // ' trans energy ='
  std::string z_energy;          // ' trans energy in z ='
  std::string total_momentum;    // ' total momentum =' (a.u.)
  std::vector<FortranVibrationRow> vibration;
  std::vector<FortranRotationRow> rotation;
  int order_headers = 0;  // ' v, j, mj =   0   0   0'
  std::vector<FortranOrderRow> order;
  std::string p_total;        // ' sum of the probabilities ='
  std::string e_vibration;    // ' average vibrational energy (ev) =' (f12.6)
  std::string e_rotation;     // ' average rotational  energy (ev) ='
  std::string e_diffraction;  // ' average diffraction energy (ev) ='
  std::string p_reaction;     // ' the reaction probability ='
  std::vector<std::string> problems;
};

// The per-energy blocks of a Fortran stdout: a block starts at its ' total energy =' line, a section at its title;
// a row is a line of the section's number of tokens that starts with an integer, until 'sum of the probabilities'.
inline std::vector<FortranEnergyBlock> read_fortran_energy_blocks(const std::filesystem::path& stdout_file) {
  std::ifstream file(stdout_file);
  if (!file) throw std::runtime_error("cannot open " + stdout_file.string());
  std::vector<FortranEnergyBlock> blocks;
  enum class Section { none, vibration, rotation, order };
  Section section = Section::none;
  std::optional<std::tuple<int, int, int>> order_state;  // the last order header's (v, j, mJ)
  std::string line;
  int number = 0;
  while (std::getline(file, line)) {
    number = number + 1;
    const auto has = [&](std::string_view text) { return line.find(text) != std::string::npos; };
    if (has(" total energy =")) {
      blocks.emplace_back();
      blocks.back().total_energy = last_token(line);
      section = Section::none;
      continue;
    }
    if (blocks.empty()) continue;
    FortranEnergyBlock& b = blocks.back();
    const std::string at = "line " + std::to_string(number);
    if (has(" trans energy =")) {
      b.collision_energy = last_token(line);
    } else if (has(" trans energy in z =")) {
      b.z_energy = last_token(line);
    } else if (has(" total momentum =")) {
      b.total_momentum = last_token(line);
    } else if (has("rotation-diffraction summed probabilities")) {
      section = Section::vibration;
    } else if (has("diffraction summed probabilities")) {
      section = Section::rotation;
    } else if (has("probabilities as function of diffraction order")) {
      section = Section::order;
      order_state.reset();
    } else if (has("sum of the probabilities =")) {  // the end of the tables
      b.p_total = last_token(line);
      section = Section::none;
    } else if (has("average vibrational energy (ev) =")) {
      b.e_vibration = last_token(line);
    } else if (has("average rotational  energy (ev) =")) {
      b.e_rotation = last_token(line);
    } else if (has("average diffraction energy (ev) =")) {
      b.e_diffraction = last_token(line);
    } else if (has("the reaction probability =")) {
      b.p_reaction = last_token(line);
    } else if (section == Section::order && has(" v, j, mj =")) {
      const std::vector<std::string> t = split(line);
      const auto v = printed_integer(t[t.size() - 3]), j = printed_integer(t[t.size() - 2]),
                 mj = printed_integer(t.back());
      b.order_headers = b.order_headers + 1;
      if (v && j && mj) {
        order_state = {*v, *j, *mj};
      } else {
        order_state.reset();  // its rows are problems too
        b.problems.push_back(at + ": an unreadable order header");
      }
    } else {
      const std::vector<std::string> t = split(line);
      if (t.empty() || !all_digits(t[0])) continue;
      const int first = std::stoi(t[0]);
      if (section == Section::vibration && t.size() == 4) {
        b.vibration.push_back({first, t[1], t[2], t[3]});
      } else if (section == Section::rotation && t.size() == 5) {
        const auto irot = printed_integer(t[1]), j = printed_integer(t[2]), mj = printed_integer(t[3]);
        if (irot && j && mj)
          b.rotation.push_back({first, *irot, *j, *mj, t[4]});
        else
          b.problems.push_back(at + ": an unreadable rotation row");
      } else if (section == Section::order && t.size() == 2) {
        if (order_state)
          b.order.push_back(
              {std::get<0>(*order_state), std::get<1>(*order_state), std::get<2>(*order_state), first, t[1]});
        else
          b.problems.push_back(at + ": an order row before its header");
      }
    }
  }
  return blocks;
}

// The noise threshold of the average energies: an average is a ratio of two sums of |S|^2 (Analyzor.f90
// smatrix_combine: av edif and av erot of P(v), the three averages of the sum of the probabilities), compared only
// where the printed probability it divides by exceeds this; below, both codes print ratios of rounding noise.
// Measured on every gate-4 case (both Long ones and the restart run included) and gate-5 case (2026-11-24): the
// printed P(v) are either rounding noise, at most 9.4e-30 (hardwall_short, whose |S| are all <= 2.4e-15; its
// averages differ by up to 7.2e5 half printed units, av erot by 4.4 %), or at least 2.8e-14 (cu111_dailight_oblique,
// gate 5), where every average passes within its printed rounding (worst 0.98 of half a unit below P(v) = 1e-11).
// The threshold lies in that gap. (The largest |P_C++(v) - P_F(v)| where P(v) < 1e-10, about 5e-18, is half a
// printed unit of an E14.7 value just below 1e-10: the Fortran's print resolution, not the codes' noise.)
inline constexpr double average_energy_noise_threshold = 5e-15;

// The counts of one section: rows (values) read and those the Fortran's print rule gives for the C++ run (summed
// over the blocks; the comparison also checks each block), rows compared, and rows skipped (the averages: the
// probability below the noise threshold).
struct SectionCounts {
  int read = 0;
  int expected = 0;
  int compared = 0;
  int skipped = 0;
};

// The outcome of a comparison of the per-energy blocks with the C++ results; failures holds the first few.
struct EnergyBlockComparison {
  int blocks = 0;
  SectionCounts scalars;  // per block: total energy, trans energy (and in z), total momentum, sum, reaction probability
  SectionCounts averages;  // per block: the average vibrational, rotational, diffraction energy
  SectionCounts vibration;  // P(v)
  SectionCounts vibration_means;  // av edif and av erot of a P(v) row
  SectionCounts rotation;  // P(v, j, mJ)
  SectionCounts order_headers;  // the (v, j, mJ) headers of the order block
  SectionCounts order;  // P(v, j, mJ, order)
  int failed = 0;  // values outside their tolerance, rows without a C++ state, counts that differ
  double worst = 0.0;  // the largest deviation in units of its tolerance
  std::vector<std::string> failures;
};

// The rows the Fortran's print rules give for one block of the C++ run at total_energy, added to c's expected
// counts: 6 scalars and 3 averages; the v block nviba rows; the rotation block nviba * nana; with writeRDProbs the
// order block one header per (v, j, mJ) with E >= E(v, j) (smatrix_combine: `cycle` if E < E(v, j)) and one row per
// diffraction order under it.
inline void add_expected_rows(const PropagationRun& cpp, double total_energy, EnergyBlockComparison& c) {
  c.scalars.expected = c.scalars.expected + 6;
  c.averages.expected = c.averages.expected + 3;
  c.vibration.expected = c.vibration.expected + cpp.n_v;
  c.vibration_means.expected = c.vibration_means.expected + cpp.n_v;
  c.rotation.expected = c.rotation.expected + cpp.n_v * static_cast<int>(cpp.rotation_labels.size());
  if (!cpp.order_block_printed) return;
  for (int v = 0; v < cpp.n_v; ++v)
    for (const auto& [j, mj] : cpp.rotation_labels)
      if (!(total_energy < cpp.level_energy_ev[v][j])) {
        c.order_headers.expected = c.order_headers.expected + 1;
        c.order.expected = c.order.expected + cpp.n_orders;
      }
}

// Compares the per-energy blocks of a Fortran stdout with the C++ run, every printed value to half a unit of its
// last printed digit (printed_text_tolerance; the reaction probability at least 1e-11 absolute, see below), the
// averages where their probability exceeds noise_threshold. Every row must match a C++ state; per block the rows
// must number what the print rule gives. A NaN on either side fails.
inline EnergyBlockComparison compare_energy_blocks(const PropagationRun& cpp,
                                                   const std::vector<FortranEnergyBlock>& ref,
                                                   double noise_threshold = average_energy_noise_threshold) {
  EnergyBlockComparison c;
  c.blocks = static_cast<int>(ref.size());
  auto fail = [&](const std::string& what) {
    c.failed = c.failed + 1;
    if (c.failures.size() < 20) c.failures.push_back(what);
  };
  // absolute_floor: the smallest tolerance, for values obtained by cancellation.
  auto check = [&](double value, const std::string& text, const std::string& what, double absolute_floor = 0.0) {
    const std::optional<double> printed = fortran_number(text);
    if (!printed) {
      c.worst = std::numeric_limits<double>::quiet_NaN();
      fail(what + (text.empty() ? ": not printed" : ": Fortran printed '" + text + "'"));
      return;
    }
    const double tolerance = std::max(printed_half_unit(text, 1e-20), absolute_floor);
    const double deviation = std::abs(value - *printed) / tolerance;
    c.worst = nan_aware_max(c.worst, deviation);
    if (!(deviation <= 1.0)) fail(what + ": C++ " + std::format("{:.16e}", value) + ", Fortran " + text);
  };
  // Whether the probability an average divides by is known to be noise (a NaN or an unreadable one is not).
  auto noise = [&](const std::string& probability) {
    const std::optional<double> p = fortran_number(probability);
    return p && *p <= noise_threshold;
  };
  if (ref.size() != cpp.scattering.energies.size())
    fail("the stdout has " + std::to_string(ref.size()) + " energy blocks, the run " +
         std::to_string(cpp.scattering.energies.size()));
  const int n_rot = static_cast<int>(cpp.rotation_labels.size());
  auto rotation_index = [&](int j, int mj) {
    const auto found = std::ranges::find(cpp.rotation_labels, std::pair{j, mj});
    return found == cpp.rotation_labels.end() ? -1 : static_cast<int>(found - cpp.rotation_labels.begin());
  };
  for (std::size_t ie = 0; ie < std::min(ref.size(), cpp.scattering.energies.size()); ++ie) {
    const EnergyResult& e = cpp.scattering.energies[ie];
    const FortranEnergyBlock& f = ref[ie];
    const std::string at = "energy " + std::to_string(ie);
    for (const std::string& problem : f.problems) fail(at + ", " + problem);

    // The rows of each section against the print rule.
    EnergyBlockComparison block;
    add_expected_rows(cpp, e.total_energy, block);
    const auto count = [&](const std::string& section, int read, int expected, SectionCounts& total) {
      total.read = total.read + read;
      total.expected = total.expected + expected;
      if (read != expected)
        fail(at + ": " + std::to_string(read) + " " + section + " rows read, the print rule gives " +
             std::to_string(expected));
    };
    const auto printed = [](const std::vector<std::string>& texts) {
      return static_cast<int>(std::ranges::count_if(texts, [](const std::string& t) { return !t.empty(); }));
    };
    count("scalar",
          printed({f.total_energy, f.collision_energy, f.z_energy, f.total_momentum, f.p_total, f.p_reaction}),
          block.scalars.expected, c.scalars);
    count("average energy", printed({f.e_vibration, f.e_rotation, f.e_diffraction}), block.averages.expected,
          c.averages);
    count("P(v)", static_cast<int>(f.vibration.size()), block.vibration.expected, c.vibration);
    c.vibration_means.read = c.vibration_means.read + static_cast<int>(f.vibration.size());
    c.vibration_means.expected = c.vibration_means.expected + block.vibration_means.expected;
    count("P(v, j, mJ)", static_cast<int>(f.rotation.size()), block.rotation.expected, c.rotation);
    count("order header", f.order_headers, block.order_headers.expected, c.order_headers);
    count("P(v, j, mJ, order)", static_cast<int>(f.order.size()), block.order.expected, c.order);
    // Per header: each (v, j, mJ) printed once, with one row per diffraction order (0 .. mdifa).
    std::map<std::tuple<int, int, int>, int> rows_per_header;
    for (const FortranOrderRow& r : f.order) rows_per_header[{r.v, r.j, r.mj}] += 1;
    if (static_cast<int>(rows_per_header.size()) != f.order_headers)
      fail(at + ": " + std::to_string(f.order_headers) + " order headers for " +
           std::to_string(rows_per_header.size()) + " (v, j, mJ) with rows");
    for (const auto& [key, rows] : rows_per_header)
      if (rows != cpp.n_orders)
        fail(at + ": " + std::to_string(rows) + " order rows under (v, j, mJ) = (" + std::to_string(std::get<0>(key)) +
             ", " + std::to_string(std::get<1>(key)) + ", " + std::to_string(std::get<2>(key)) + "), " +
             std::to_string(cpp.n_orders) + " orders");

    // The scalars. The reaction probability is 1 - sum (Fortran: 1.0 - sump). When almost everything
    // scatters it is the small difference of two numbers near 1, so it is known only to the *absolute*
    // accuracy of the sum, not to the 7 significant digits the Fortran prints. Example
    // (cu111_dailight_small, E = 0.42 eV): Fortran prints 0.2943414D-05, i.e. a half-unit
    // tolerance of 5e-14, while the sum itself (~0.999997) is only reproducible to ~1e-13..1e-12
    // (different FFT library and summation order than the Fortran, rounding accumulated over
    // 7000 steps). Hence an absolute floor instead of the printed-digit tolerance.
    //
    // The floor is 1e-11 (it was 1e-12). Deviations from the Fortran (cu111_dailight_small):
    // - default build (no fast-math): 3.6e-13 at E = 0.42 eV (2.9434149e-6, printed half unit
    //   5e-14) and 4.1e-12 at E = 0.44 eV (-3.6025106e-5, printed half unit 5e-12);
    // - SPODVR_FAST_MATH=ON (-ffast-math: FMA contraction, reassociation, vectorised math,
    //   measured 2026-11-05 with -march=native and PGO): the norm moves by ~2e-12 relative and
    //   the scattered sum by ~1.5e-12 absolute against the default build, giving 2.0e-12 at
    //   E = 0.42 eV and 5.3e-12 at E = 0.44 eV (just over the printed half unit 5e-12: the
    //   Fortran value itself lies near a rounding boundary of its 7th digit).
    // These are rounding-level differences of the sum, far below any physical significance
    // (the probabilities themselves are 3e-6 and 4e-5), and the sum, P(v) and P(v, j, mJ) are
    // still checked to half a printed unit. 1e-11 gives a margin of ~2 over the fast-math
    // build; a deviation beyond it points to a real change, not to rounding.
    check(e.total_energy, f.total_energy, at + " total energy");
    check(e.collision_energy, f.collision_energy, at + " trans energy");
    check(e.z_energy, f.z_energy, at + " trans energy in z");
    check(e.total_momentum, f.total_momentum, at + " total momentum");
    check(e.p_total, f.p_total, at + " sum of the probabilities");
    check(e.p_reaction, f.p_reaction, at + " reaction probability", 1e-11);
    c.scalars.compared = c.scalars.compared + 6;
    // The averages over every scattered molecule (f12.6, eV), ratios to the sum of the probabilities: noise if the
    // printed P(v) (E14.7; the printed sum has 7 decimals only) add up to at most the threshold (a NaN is not).
    double printed_sum = 0.0;
    for (const FortranVibrationRow& r : f.vibration)
      printed_sum = printed_sum + fortran_number(r.p).value_or(std::numeric_limits<double>::quiet_NaN());
    if (printed_sum <= noise_threshold) {
      c.averages.skipped = c.averages.skipped + 3;
      // Not compared, but a C++ average must still be a number.
      if (!std::isfinite(e.e_vibration) || !std::isfinite(e.e_rotation) || !std::isfinite(e.e_diffraction))
        fail(at + ": a non-finite C++ average energy");
    } else {
      check(e.e_vibration, f.e_vibration, at + " average vibrational energy");
      check(e.e_rotation, f.e_rotation, at + " average rotational energy");
      check(e.e_diffraction, f.e_diffraction, at + " average diffraction energy");
      c.averages.compared = c.averages.compared + 3;
    }

    // P(v) and the average energies of the molecules scattered into v (ratios to P(v)).
    for (const FortranVibrationRow& r : f.vibration) {
      const std::string row = at + " v = " + std::to_string(r.v);
      if (r.v < 0 || r.v >= cpp.n_v) {
        fail(row + ": no C++ vibrational level");
        continue;
      }
      check(e.p_vibration[r.v], r.p, row + " P(v)");
      c.vibration.compared = c.vibration.compared + 1;
      if (noise(r.p)) {
        c.vibration_means.skipped = c.vibration_means.skipped + 1;
        if (!std::isfinite(e.e_diffraction_v[r.v]) || !std::isfinite(e.e_rotation_v[r.v]))
          fail(row + ": a non-finite C++ av edif or av erot");
        continue;
      }
      check(e.e_diffraction_v[r.v], r.e_diffraction, row + " av edif");
      check(e.e_rotation_v[r.v], r.e_rotation, row + " av erot");
      c.vibration_means.compared = c.vibration_means.compared + 1;
    }
    // P(v, j, mJ): the Fortran's irot is the C++ rotational index + 1.
    for (const FortranRotationRow& r : f.rotation) {
      const std::string row = at + " (v, irot, j, mJ) = (" + std::to_string(r.v) + ", " + std::to_string(r.irot) +
                              ", " + std::to_string(r.j) + ", " + std::to_string(r.mj) + ")";
      const int irotana = r.irot - 1;
      if (r.v < 0 || r.v >= cpp.n_v || irotana < 0 || irotana >= n_rot ||
          cpp.rotation_labels[irotana] != std::pair{r.j, r.mj}) {
        fail(row + ": no C++ rotational state");
        continue;
      }
      check(e.p_rotation[irotana + n_rot * r.v], r.p, row + " P(v, j, mJ)");
      c.rotation.compared = c.rotation.compared + 1;
    }
    // P(v, j, mJ, order), printed for the open (v, j) only.
    for (const FortranOrderRow& r : f.order) {
      const std::string row = at + " (v, j, mJ, order) = (" + std::to_string(r.v) + ", " + std::to_string(r.j) + ", " +
                              std::to_string(r.mj) + ", " + std::to_string(r.order) + ")";
      const int irotana = rotation_index(r.j, r.mj);
      if (r.v < 0 || r.v >= cpp.n_v || irotana < 0 || r.order < 0 || r.order >= cpp.n_orders) {
        fail(row + ": no C++ state");
        continue;
      }
      if (e.total_energy < cpp.level_energy_ev[r.v][r.j]) {
        fail(row + ": a closed channel (total energy below E(v, j))");
        continue;
      }
      check(e.p_order[r.order + cpp.n_orders * (irotana + n_rot * r.v)], r.p, row + " P(v, j, mJ, order)");
      c.order.compared = c.order.compared + 1;
    }
  }
  return c;
}

// The gate's assertions on a comparison of the energy blocks: no failure, every section's rows as the print rule
// gives, each compared (the averages: compared or skipped as noise).
inline void expect_energy_blocks_gate(const EnergyBlockComparison& c) {
  EXPECT_EQ(c.failed, 0) << "failed values, rows or counts in the energy blocks (the first 20):" << listed(c.failures);
  EXPECT_GT(c.blocks, 0) << "no energy block read";
  const auto expect_section = [](const std::string& name, const SectionCounts& s) {
    EXPECT_EQ(s.read, s.expected) << name << " rows read vs the print rule";
    EXPECT_EQ(s.compared + s.skipped, s.read) << name << " rows compared or skipped vs rows read";
  };
  expect_section("scalar", c.scalars);
  expect_section("average energy", c.averages);
  expect_section("P(v)", c.vibration);
  expect_section("av edif / av erot", c.vibration_means);
  expect_section("P(v, j, mJ)", c.rotation);
  expect_section("P(v, j, mJ, order)", c.order);
  EXPECT_EQ(c.order_headers.read, c.order_headers.expected) << "order headers read vs the print rule";
}

// Gate 4 (and 5) on the per-energy blocks of a Fortran stdout: compare_energy_blocks, its counts printed and asserted.
inline void compare_energy_blocks_with_fortran(const PropagationRun& cpp, const std::vector<FortranEnergyBlock>& ref) {
  const EnergyBlockComparison c = compare_energy_blocks(cpp, ref);
  std::cout << "  " << c.blocks << " energy blocks: " << c.scalars.compared << " scalars, " << c.averages.compared
            << " average energies (" << c.averages.skipped << " skipped), " << c.vibration.compared << " P(v), "
            << c.vibration_means.compared << " av edif/erot rows (" << c.vibration_means.skipped << " skipped), "
            << c.rotation.compared << " P(v, j, mJ), " << c.order.compared << " P per order under "
            << c.order_headers.read << " headers (skipped: probability <= " << average_energy_noise_threshold
            << "); largest deviation " << c.worst << " x (half printed unit)\n";
  expect_energy_blocks_gate(c);
}

// The BK coefficient records of a Fortran fort.21 (after its ten header records).
inline std::vector<std::vector<std::complex<double>>> read_fortran_bk_records(const std::filesystem::path& path) {
  FortranRecordReader file(path);
  for (int k = 1; k <= 10; ++k) (void)file.read("header record " + std::to_string(k));
  std::vector<std::vector<std::complex<double>>> records;
  while (!file.at_end()) {
    FortranRecordParser record = file.read("BK record");
    std::vector<std::complex<double>> values(record.size() / sizeof(std::complex<double>));
    record.read_array(std::span<std::complex<double>>(values));
    record.expect_end();
    records.push_back(std::move(values));
  }
  return records;
}

// The BK deviation of a C++ run from the Fortran fort.21: the largest reference coefficient, the
// gate-4 scale (that coefficient, at least 1e-4) and the largest |C++ - Fortran|, over the analysis
// steps and coefficients both runs have (the callers check the counts).
// A NaN coefficient on either side makes `worst` (in the reference also the scale) NaN.
struct BkDeviation {
  double largest = 0.0;
  double scale = 0.0;
  double worst = 0.0;
};

inline BkDeviation bk_deviation(const PropagationRun& cpp, const std::vector<std::vector<std::complex<double>>>& ref) {
  BkDeviation d;
  for (const auto& record : ref)
    for (const auto& c : record) d.largest = nan_aware_max(d.largest, std::abs(c));
  d.scale = nan_aware_max(d.largest, 1e-4);
  for (std::size_t step = 0; step < std::min(cpp.bk.size(), ref.size()); ++step)
    for (std::size_t k = 0; k < std::min(cpp.bk[step].size(), ref[step].size()); ++k)
      d.worst = nan_aware_max(d.worst, std::abs(cpp.bk[step][k] - ref[step][k]));
  return d;
}

// Gate 4 (BK coefficients): every coefficient of every analysis step, relative to the largest
// coefficient of the run. The coefficients are amplitudes of at most ~1; before the wave packet
// reaches the analysis line they are rounding noise (hardwall_short: all below 2e-15), so the
// reference scale is at least 1e-4, i.e. noise is compared absolutely at tolerance * 1e-4.
// A NaN coefficient on either side makes the deviation NaN: it fails.
inline void compare_bk_with_fortran(const PropagationRun& cpp,
                                    const std::vector<std::vector<std::complex<double>>>& ref, double tolerance) {
  ASSERT_EQ(cpp.bk.size(), ref.size());
  for (std::size_t step = 0; step < ref.size(); ++step) ASSERT_EQ(cpp.bk[step].size(), ref[step].size());
  const BkDeviation d = bk_deviation(cpp, ref);
  EXPECT_LE(d.worst, tolerance * d.scale);
  std::cout << "  " << ref.size() << " BK records compared; max|c| = " << d.largest
            << ", largest |C++ - Fortran| = " << d.worst << " (" << d.worst / d.scale << " of the reference scale)\n";
}

// The specular analysis line as the Fortran stdout prints it when the line becomes fixed (Analyzor.f90
// FindNewAnalysisPoint, 5.19 :254 and :263: '==> ZAnalysis on specular channel fixed at Z = ', e15.8): the printed
// value of every such line.
inline std::vector<std::string> read_fortran_fixed_line(const std::filesystem::path& stdout_file) {
  std::ifstream file(stdout_file);
  if (!file) throw std::runtime_error("cannot open " + stdout_file.string());
  std::vector<std::string> printed;
  std::string line;
  while (std::getline(file, line))
    if (line.find("ZAnalysis on specular channel fixed at Z =") != std::string::npos)
      printed.push_back(last_token(line));
  return printed;
}

// Gate 3: the C++ specular analysis line becomes fixed if and only if the Fortran's does (printed at most once; the
// Fortran prints no time), at the same Z to half a printed unit (a line that never fixes, as in hardwall_short, is
// absent in both).
inline void expect_fixed_line(const std::optional<double>& cpp_line, const std::vector<std::string>& printed) {
  EXPECT_LE(printed.size(), 1u) << "the Fortran printed its fixed specular line " << printed.size() << " times";
  EXPECT_EQ(cpp_line.has_value(), !printed.empty())
      << "specular line fixed: C++ " << (cpp_line ? "yes" : "no") << ", Fortran " << (printed.empty() ? "no" : "yes");
  if (printed.size() != 1 || !cpp_line) return;
  const std::optional<double> fortran = fortran_number(printed[0]);
  EXPECT_TRUE(fortran.has_value()) << "Fortran fixed line '" << printed[0] << "' is not a finite number";
  if (!fortran) return;
  EXPECT_LE(std::abs(*cpp_line - *fortran), printed_text_tolerance(printed[0]))
      << "fixed specular line: C++ Z = " << *cpp_line << ", Fortran " << printed[0];
}

}  // namespace spodvr::testing

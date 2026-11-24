#pragma once

// Gate 4, the S-matrix elements: the Fortran's on-the-fly fort.25 (Analyzor.f90 smatrix_combine, unit 25: Re S,
// Im S, P and the phase of every printed state) against the C++ ScatteringResults, element by element, with the
// counts the comparison must reach (tests/regression_gate.hpp asserts them).

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <optional>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include "analysis_reference.hpp"
#include "comparison.hpp"
#include "propagation_reference.hpp"

namespace spodvr::testing {

// One printed row of an on-the-fly fort.25: its energy block (0-based, in the order of the analysis energies),
// the block header's (v, j, mJ, diffraction order), the row's diffraction indices (n, m) and Re S, Im S, P and
// the phase as printed, and its line (1-based).
struct Fort25Row {
  int block = 0;
  int v = 0, j = 0, mj = 0, order = 0, n = 0, m = 0;
  std::string re, im, p, phase;
  int line = 0;
};

// An on-the-fly fort.25: the rows, the total energy of each block as printed (eV), an independent count of
// the lines that look like a row (a pattern of their own, not the reader's), and the lines the reader could
// not use (a number that is not finite or not readable, a row outside a block, a repeated state).
struct Fort25File {
  std::vector<Fort25Row> rows;
  std::vector<std::string> total_energy;
  int row_like_lines = 0;
  std::vector<std::string> problems;
};

// A printed Fortran real with an E exponent: D becomes E, and a three-digit exponent without its letter (gfortran
// writes 0.100000-119) gets it back, so that printed_text_tolerance reads its exponent.
inline std::string fortran_text_with_e(const std::string& text) {
  static const std::regex no_letter(R"(([-+]?\d*\.\d*)([-+]\d{3}))");
  std::string t = text;
  for (char& c : t)
    if (c == 'D' || c == 'd') c = 'E';
  std::smatch parts;
  if (std::regex_match(t, parts, no_letter)) t = parts[1].str() + "E" + parts[2].str();
  return t;
}

// The value of a printed Fortran real, or nothing if the text is not a finite number (NaN, Infinity, ****).
inline std::optional<double> fortran_number(const std::string& text) {
  static const std::regex number(R"([-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?)");
  const std::string t = fortran_text_with_e(text);
  if (!std::regex_match(t, number)) return std::nullopt;
  try {
    const double x = std::stod(t);
    return std::isfinite(x) ? std::optional<double>(x) : std::nullopt;
  } catch (const std::out_of_range&) {
    return std::nullopt;
  }
}

inline Fort25File read_fortran_fort25_rows(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path.string());
  static const std::regex energy(R"( total energy = +(\S+) *)");
  static const std::regex header(R"( v j mj dif\.order are: +(-?\d+) +(-?\d+) +(-?\d+) +(-?\d+) *)");
  static const std::regex row(R"( +(-?\d+) +(-?\d+) +(\S+) +(\S+) +(\S+) +(\S+) *)");
  static const std::regex row_like(R"( +-?[0-9]+ +-?[0-9]+ +[-0-9.]+(D[-+][0-9]{2}|[-+][0-9]{3}) .*)");
  Fort25File file;
  std::optional<std::tuple<int, int, int, int>> state;  // the last block header's (v, j, mJ, order)
  std::set<std::tuple<int, int, int, int, int, int>> seen;
  std::string line;
  int number = 0;
  std::smatch match;
  while (std::getline(in, line)) {
    number = number + 1;
    if (std::regex_match(line, row_like)) file.row_like_lines = file.row_like_lines + 1;
    if (std::regex_match(line, match, energy)) {
      file.total_energy.push_back(match[1].str());
      state.reset();
    } else if (std::regex_match(line, match, header)) {
      state = {std::stoi(match[1].str()), std::stoi(match[2].str()), std::stoi(match[3].str()),
               std::stoi(match[4].str())};
    } else if (std::regex_match(line, match, row)) {
      const std::string at = "line " + std::to_string(number);
      if (file.total_energy.empty() || !state) {
        file.problems.push_back(at + ": a row outside a block");
        continue;
      }
      const auto [v, j, mj, order] = *state;
      Fort25Row r{static_cast<int>(file.total_energy.size()) - 1,
                  v,
                  j,
                  mj,
                  order,
                  std::stoi(match[1].str()),
                  std::stoi(match[2].str()),
                  match[3].str(),
                  match[4].str(),
                  match[5].str(),
                  match[6].str(),
                  number};
      for (const std::string& token : {r.re, r.im, r.p, r.phase})
        if (!fortran_number(token)) file.problems.push_back(at + ": '" + token + "' is not a finite number");
      if (!seen.insert({r.block, v, j, mj, r.n, r.m}).second) file.problems.push_back(at + ": the state repeats");
      file.rows.push_back(r);
    }
  }
  return file;
}

// Half a unit in the last printed digit of a Fortran real (printed_text_tolerance, at least floor), but floor
// for a printed zero: an E-format zero is exact (a nonzero value always prints a nonzero mantissa).
inline double printed_half_unit(const std::string& text, double floor) {
  const std::string t = fortran_text_with_e(text);
  const std::string mantissa = t.substr(0, t.find('E'));
  if (mantissa.find_first_of("123456789") == std::string::npos) return floor;
  return printed_text_tolerance(t, floor);
}

// The tolerances of the S-matrix elements (measured, tests/regression_gate.hpp). Two errors add: the printed
// digits of the Fortran (half a unit of Re and Im, 0 for a printed zero) and the difference of the two codes, at
// most code = s_scale_fraction * max(the case's largest |S_F|, s_floor), an absolute error. |S_C - S_F| within
// their sum; P within half its printed unit (noise floor 1e-20) plus the code difference's effect on |S|^2; the
// phase within half its printed unit plus code / |S| (the Fortran prints the phase of its unrounded S, so the
// digits of Re and Im do not enter), compared only where that bound is at most phase_cap (|S| >~ 1e3 code).
struct Fort25Tolerance {
  double s_scale_fraction = 0.0;
  // The floor of the scale: in a case whose elements are all rounding noise (hardwall_short, |S| <= 2.4e-15) the
  // code term is s_scale_fraction * s_floor; above it, S elements are amplitudes of a unit-norm packet.
  double s_floor = 1e-3;
  double phase_cap = 1e-3;
};

// The counts gate 4 must reach in a case (measured, tests/regression_gate.hpp): the rows whose phase is
// compared, and the rows that fail when the C++ refers the initial state's phase to Zana for every state.
struct Fort25Pins {
  int phase_compared;
  int zana_failures;
};

// The outcome of one comparison; failures holds the first few, with their lines.
struct Fort25Comparison {
  int rows = 0;               // rows read
  int expected_rows = 0;      // the Fortran's print rule applied to the C++ run (independent of the file)
  int row_like_lines = 0;     // the independent count of the file
  int compared = 0;           // rows matched to a C++ state and compared
  int phase_compared = 0;     // of those, rows whose phase was compared
  int failed = 0;             // rows outside a tolerance, unmatched, or with a problem
  double s_scale = 0.0;       // the largest |S| of the file
  double worst_ratio = 0.0;   // max (|S_C - S_F| - printed digits) / max(s_scale, s_floor): the code difference
  std::vector<std::string> failures;
};

// Compares the C++ run with a read fort.25. A NaN on either side fails its row (the comparisons are written so).
inline Fort25Comparison compare_fort25(const PropagationRun& cpp, const Fort25File& file,
                                       const Fort25Tolerance& tolerance) {
  Fort25Comparison c;
  c.rows = static_cast<int>(file.rows.size());
  c.row_like_lines = file.row_like_lines;
  auto fail = [&](const std::string& what) {
    c.failed = c.failed + 1;
    if (c.failures.size() < 20) c.failures.push_back(what);
  };
  for (const std::string& problem : file.problems) fail(problem);
  if (file.total_energy.size() != cpp.scattering.energies.size())
    fail("the file has " + std::to_string(file.total_energy.size()) + " energy blocks, the run " +
         std::to_string(cpp.scattering.energies.size()));
  for (const EnergyResult& e : cpp.scattering.energies)
    for (int v = 0; v < cpp.n_v; ++v)
      for (int is = 0; is < cpp.n_states; ++is)
        if (e.total_energy > cpp.level_energy_ev[v][cpp.states[is].j]) c.expected_rows = c.expected_rows + 1;
  for (const Fort25Row& r : file.rows)
    if (const auto re = fortran_number(r.re), im = fortran_number(r.im); re && im)
      c.s_scale = nan_aware_max(c.s_scale, std::abs(std::complex<double>(*re, *im)));
  const double scale = nan_aware_max(c.s_scale, tolerance.s_floor);
  // Each block's total energy, as printed (E15.8: half a unit of the 8th digit).
  for (std::size_t b = 0; b < file.total_energy.size() && b < cpp.scattering.energies.size(); ++b) {
    const std::optional<double> printed = fortran_number(file.total_energy[b]);
    const double cpp_energy = cpp.scattering.energies[b].total_energy;
    if (!printed ||
        !(std::abs(cpp_energy - *printed) <= printed_text_tolerance(fortran_text_with_e(file.total_energy[b]), 0.0)))
      fail("block " + std::to_string(b + 1) + ": total energy C++ " + std::to_string(cpp_energy) + ", Fortran " +
           file.total_energy[b]);
  }
  // With a single diffraction order the Fortran prints one block header per v for all rotational states
  // (Analyzor.f90 r1727 lines 888-896: the header follows a change of the order), so rows carry no (j, mJ).
  bool several_orders = false;
  for (const StateLabel& s : cpp.states) several_orders = several_orders || s.order != cpp.states.front().order;
  if (!several_orders && !file.rows.empty())
    fail("one diffraction order: the Fortran's fort.25 does not label the rotational states (not comparable)");
  for (const Fort25Row& r : file.rows) {
    const std::string at = "line " + std::to_string(r.line);
    if (r.block >= static_cast<int>(cpp.scattering.energies.size()) || r.v < 0 || r.v >= cpp.n_v) {
      fail(at + ": no C++ energy or level for it");
      continue;
    }
    int is = -1;
    for (int k = 0; k < cpp.n_states; ++k) {
      const StateLabel& s = cpp.states[k];
      if (s.j == r.j && s.mj == r.mj && s.n == r.n && s.m == r.m) is = k;
    }
    if (is < 0 || cpp.states[is].order != r.order) {
      fail(at + ": no C++ state (v, j, mJ, order, n, m) = (" + std::to_string(r.v) + ", " + std::to_string(r.j) + ", " +
           std::to_string(r.mj) + ", " + std::to_string(r.order) + ", " + std::to_string(r.n) + ", " +
           std::to_string(r.m) + ")");
      continue;
    }
    const EnergyResult& e = cpp.scattering.energies[r.block];
    if (!(e.total_energy > cpp.level_energy_ev[r.v][r.j])) {  // the Fortran prints open channels only
      fail(at + ": a closed channel (total energy not above E(v, j))");
      continue;
    }
    const int k = is + cpp.n_states * r.v;
    const std::optional<double> re = fortran_number(r.re), im = fortran_number(r.im), p = fortran_number(r.p),
                                phase = fortran_number(r.phase);
    if (!re || !im || !p || !phase) continue;  // counted as a problem of the file
    c.compared = c.compared + 1;
    const std::complex<double> s_f(*re, *im);
    const std::complex<double> s_c = e.s_matrix[k];
    const double digits = std::hypot(printed_half_unit(r.re, 0.0), printed_half_unit(r.im, 0.0));
    const double code = tolerance.s_scale_fraction * scale;  // the allowed difference of the two codes
    const double s_error = std::abs(s_c - s_f);
    if (s_error > digits) c.worst_ratio = nan_aware_max(c.worst_ratio, (s_error - digits) / scale);
    const double s_abs = std::abs(s_f);
    const double p_tolerance = printed_half_unit(r.p, 1e-20) + 2.0 * (s_abs + digits) * code + code * code;
    // The C++ phase must be finite on every row, also where it is not compared.
    bool ok = s_error <= digits + code && std::abs(e.probability[k] - *p) <= p_tolerance && std::isfinite(e.phase[k]);
    if (*re != 0.0 && s_abs > 0.0) {
      const double phase_tolerance = printed_half_unit(r.phase, 0.0) + code / s_abs;
      if (phase_tolerance <= tolerance.phase_cap) {
        c.phase_compared = c.phase_compared + 1;
        const double difference = std::remainder(e.phase[k] - *phase, 2.0 * std::numbers::pi);
        ok = ok && std::abs(difference) <= phase_tolerance;
      }
    }
    if (!ok)
      fail(at + ": S C++ " + std::to_string(s_c.real()) + " " + std::to_string(s_c.imag()) + ", Fortran " + r.re + " " +
           r.im + "; P C++ " + std::to_string(e.probability[k]) + ", Fortran " + r.p + "; phase C++ " +
           std::to_string(e.phase[k]) + ", Fortran " + r.phase);
  }
  return c;
}

// Gate 4's assertions on one comparison: no failed row, the rows read equal to the Fortran's print rule and to the
// independent count of the file, every row compared, and the pinned counts of the case (zana_failures: the failed
// rows of the same comparison with the initial state's phase referred to Zana).
inline void expect_fort25_gate(const Fort25Comparison& c, int zana_failures, const Fort25Pins& pins) {
  for (const std::string& failure : c.failures) ADD_FAILURE() << "fort.25 " << failure;
  EXPECT_EQ(c.failed, 0) << "failed rows";
  EXPECT_GT(c.rows, 0) << "no rows read";
  EXPECT_EQ(c.rows, c.expected_rows) << "rows read vs the Fortran's print rule";
  EXPECT_EQ(c.rows, c.row_like_lines) << "rows read vs an independent count of the file";
  EXPECT_EQ(c.compared, c.rows) << "rows compared vs rows read";
  EXPECT_EQ(c.phase_compared, pins.phase_compared) << "phases compared vs the case's pin";
  EXPECT_EQ(zana_failures, pins.zana_failures)
      << "failures with the initial state's phase referred to Zana (the Fortran up to r1719) vs the case's pin";
}

}  // namespace spodvr::testing

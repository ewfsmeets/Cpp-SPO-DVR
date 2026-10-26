// Shared code of the analysis regression test (gate 4): the probability blocks of the Fortran
// stdout (Analyzor.f90 smatrix_combine) and the Fortran BK coefficient file fort.21.
#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

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

// One per-energy block of the Fortran stdout.
struct FortranEnergyBlock {
  std::string total_energy;                                 // ' total energy =  0.45837825D+00'
  std::map<int, std::string> p_vibration;                   // v -> probability
  std::map<std::tuple<int, int, int>, std::string> p_rotation;  // (v, j, mJ) -> probability
  std::string p_total;                                      // ' sum of the probabilities ='
  std::string p_reaction;                                   // ' the reaction probability ='
};

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

inline std::vector<FortranEnergyBlock> read_fortran_probabilities(const std::filesystem::path& stdout_file) {
  std::ifstream file(stdout_file);
  std::vector<FortranEnergyBlock> blocks;
  std::string section;
  std::string line;
  while (std::getline(file, line)) {
    const auto has = [&](const char* text) { return line.find(text) != std::string::npos; };
    if (has(" total energy =")) {
      blocks.emplace_back();
      blocks.back().total_energy = last_token(line);
      section.clear();
      continue;
    }
    if (blocks.empty()) continue;
    FortranEnergyBlock& b = blocks.back();
    if (has("rotation-diffraction summed probabilities")) section = "v";
    else if (has("diffraction summed probabilities")) section = "rotation";
    else if (has("probabilities as function of diffraction order")) section.clear();
    else if (has("sum of the probabilities =")) b.p_total = last_token(line);
    else if (has("the reaction probability =")) b.p_reaction = last_token(line);
    const std::vector<std::string> t = split(line);
    if (section == "v" && t.size() == 4 && all_digits(t[0])) b.p_vibration[std::stoi(t[0])] = t[1];
    if (section == "rotation" && t.size() == 5 && all_digits(t[0]))
      b.p_rotation[{std::stoi(t[0]), std::stoi(t[2]), std::stoi(t[3])}] = t[4];
  }
  return blocks;
}

// Gate 4 (probabilities): every printed probability to all printed digits.
inline void compare_probabilities_with_fortran(const PropagationRun& cpp, const std::vector<FortranEnergyBlock>& ref) {
  ASSERT_EQ(cpp.scattering.energies.size(), ref.size());
  double worst = 0.0;
  int n_compared = 0;
  auto check = [&](double value, const std::string& text, const std::string& what) {
    const double tolerance = printed_text_tolerance(text);
    const double deviation = std::abs(value - parse_fortran_real(text)) / tolerance;
    EXPECT_LE(deviation, 1.0) << what << ": C++ " << value << ", Fortran " << text;
    worst = std::max(worst, deviation);
    n_compared = n_compared + 1;
  };
  const int n_rot = static_cast<int>(cpp.rotation_labels.size());
  for (std::size_t ie = 0; ie < ref.size(); ++ie) {
    const EnergyResult& e = cpp.scattering.energies[ie];
    const FortranEnergyBlock& f = ref[ie];
    const std::string at = "energy " + std::to_string(ie);
    check(e.total_energy, f.total_energy, at + " total energy");
    check(e.p_total, f.p_total, at + " sum of the probabilities");
    check(e.p_reaction, f.p_reaction, at + " reaction probability");
    for (const auto& [v, text] : f.p_vibration) check(e.p_vibration[v], text, at + " P(v)");
    for (int irotana = 0; irotana < n_rot; ++irotana)
      for (const auto& [key, text] : f.p_rotation) {
        const auto& [v, j, mj] = key;
        if (cpp.rotation_labels[irotana] == std::pair{j, mj})
          check(e.p_rotation[irotana + n_rot * v], text, at + " P(v, j, mJ)");
      }
  }
  std::cout << "  " << n_compared << " printed probabilities compared; largest deviation " << worst
            << " x (half printed unit)\n";
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

// Gate 4 (BK coefficients): every coefficient of every analysis step, relative to the largest
// coefficient of the run.
inline void compare_bk_with_fortran(const PropagationRun& cpp, const std::vector<std::vector<std::complex<double>>>& ref,
                                    double tolerance) {
  ASSERT_EQ(cpp.bk.size(), ref.size());
  double largest = 0.0;
  for (const auto& record : ref)
    for (const auto& c : record) largest = std::max(largest, std::abs(c));
  ASSERT_GT(largest, 0.0);
  double worst = 0.0;
  for (std::size_t step = 0; step < ref.size(); ++step) {
    ASSERT_EQ(cpp.bk[step].size(), ref[step].size());
    for (std::size_t k = 0; k < ref[step].size(); ++k)
      worst = std::max(worst, std::abs(cpp.bk[step][k] - ref[step][k]) / largest);
  }
  EXPECT_LE(worst, tolerance);
  std::cout << "  " << ref.size() << " BK records compared; largest |C++ - Fortran| / max|c| = " << worst << '\n';
}

}  // namespace spodvr::testing

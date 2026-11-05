// Readers for the Fortran fine-analysis output (gate 5): the stdout of ReadOrWriteAsympInfo and
// probaf, fil10 (fort.10), unit 3 (fort.3) and fort.25. Values are kept as printed text, so they
// can be compared to half a unit in the last printed digit (printed_text_tolerance).
#pragma once

#include <gtest/gtest.h>

#include <cctype>
#include <complex>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "analysis_reference.hpp"

namespace spodvr::testing {

inline bool is_number(const std::string& s) {
  if (s.empty()) return false;
  for (const char c : s)
    if (!(std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '+' || c == 'E' || c == 'D'))
      return false;
  return std::isdigit(static_cast<unsigned char>(s.back())) != 0;
}

// The lines of a formatted file: numeric lines (all tokens numbers) and text lines. Numbers are
// read across line boundaries; integer label lines can be recognised with peek_labels().
class NumberStream {
 public:
  explicit NumberStream(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path.string());
    std::string line;
    while (std::getline(in, line)) {
      std::vector<std::string> tokens = split(line);
      if (tokens.empty()) continue;
      bool numeric = true;
      for (const auto& t : tokens) numeric = numeric && is_number(t);
      if (numeric)
        lines_.push_back({std::move(tokens), ""});
      else
        lines_.push_back({{}, line});
    }
  }
  [[nodiscard]] bool at_end() const { return line_ >= lines_.size(); }
  std::string number() {
    if (at_end() || lines_[line_].tokens.empty())
      throw std::runtime_error("expected a number at line " + std::to_string(line_));
    std::string v = lines_[line_].tokens[token_];
    token_ = token_ + 1;
    if (token_ == lines_[line_].tokens.size()) next_line();
    return v;
  }
  std::vector<std::string> numbers(int n) {
    std::vector<std::string> v;
    for (int i = 0; i < n; ++i) v.push_back(number());
    return v;
  }
  int integer() { return std::stoi(number()); }
  void text(const std::string& contains) {
    if (at_end() || !lines_[line_].tokens.empty() || lines_[line_].text.find(contains) == std::string::npos)
      throw std::runtime_error("expected the text '" + contains + "' at line " + std::to_string(line_));
    next_line();
  }
  // Number of tokens of the current line if it is a line of integers (labels), else 0.
  [[nodiscard]] std::size_t peek_labels() const {
    if (at_end() || token_ != 0) return 0;
    for (const auto& t : lines_[line_].tokens)
      if (!all_digits(t[0] == '-' ? t.substr(1) : t)) return 0;
    return lines_[line_].tokens.size();
  }

 private:
  struct Line {
    std::vector<std::string> tokens;  // empty for a text line
    std::string text;
  };
  void next_line() {
    line_ = line_ + 1;
    token_ = 0;
  }
  std::vector<Line> lines_;
  std::size_t line_ = 0;
  std::size_t token_ = 0;
};

// ReadOrWriteAsympInfo (nread = 2) and the probaf table of the Fortran stdout.
struct FortranAsympInfo {
  int n_records = 0;           // '# timesteps for which bk coefs available is:'
  int n_used = 0;              // '# data used, ndused ='
  std::string spres_mev;       // 'spectral resolution (mev) =' (E14.7)
  std::string emin, emax;      // 'meaningful results in kin energy range (ev)' (F10.6)
  int nen95 = -1;
  // probaf table: idat -> (ekin, preact, phase)
  std::map<int, std::tuple<std::string, std::string, std::string>> table;
  int n_skip_messages = 0;     // 'etrans outside 95% conf int' + 'etrans(icomin)='
};

inline FortranAsympInfo read_fortran_asymp_info(const std::filesystem::path& stdout_file) {
  FortranAsympInfo info;
  std::ifstream in(stdout_file);
  std::string line;
  bool in_table = false;
  while (std::getline(in, line)) {
    const std::vector<std::string> t = split(line);
    if (line.find("bk coefs available is:") != std::string::npos) info.n_records = std::stoi(t.back());
    if (line.find("# data used, ndused =") != std::string::npos) info.n_used = std::stoi(t.back());
    if (line.find("spectral resolution (mev) =") != std::string::npos) info.spres_mev = t.back();
    if (line.find("meaningful results in kin energy range") != std::string::npos) {
      info.emin = t[t.size() - 2];
      info.emax = t.back();
    }
    if (line.find("nen95 energies,nen95 =") != std::string::npos) info.nen95 = std::stoi(t.back());
    if (line.find("etrans outside 95% conf int") != std::string::npos ||
        line.find("etrans(icomin)=") != std::string::npos)
      info.n_skip_messages = info.n_skip_messages + 1;
    if (line.find("idat    ekin(ev)    preact     phase") != std::string::npos) {
      in_table = true;
      continue;
    }
    if (in_table) {
      if (t.empty()) {
        if (!info.table.empty()) in_table = false;
        continue;
      }
      if (t.size() != 4 || !all_digits(t[0])) {
        in_table = false;
        continue;
      }
      info.table[std::stoi(t[0])] = {t[1], t[2], t[3]};
    }
  }
  return info;
}

// fil10 of the fine analysis (Analyzor.f90 probaf). Per-energy arrays of nen95 printed values.
struct FortranFil10 {
  int nen95 = 0;
  int n_v = 0;
  std::vector<std::string> energies;     // fetrz
  std::vector<std::string> p_initial;    // prben of the initial channel
  std::vector<std::string> phase_initial;  // phaen of the initial channel
  std::vector<std::string> p_reaction;   // prtot
  std::vector<std::string> e_vib, e_rot, e_dif;  // evav, roav, edav
  std::map<int, std::vector<std::string>> p_v;   // pv(:, v)
  std::map<int, std::vector<std::string>> e_dif_v, e_rot_v;  // edavv, eravv
  std::map<std::tuple<int, int, int>, std::vector<std::string>> p_rovib;  // (v, j, mJ) -> pvr
  std::map<std::tuple<int, int, int, int>, std::vector<std::string>> p_order;  // (order, v, j, mJ) -> pvrd
};

inline FortranFil10 read_fortran_fil10(const std::filesystem::path& path, bool with_orders) {
  NumberStream s(path);
  FortranFil10 f;
  (void)s.numbers(3);  // vin, jin, mjin
  f.nen95 = s.integer();
  (void)s.numbers(3);  // dE, Emin, Emax
  const int n_v = s.integer() + 1;
  f.n_v = n_v;
  (void)s.numbers(5);  // mdifa, jmin, jana, dJ, mjmax
  const int n = f.nen95;
  f.energies = s.numbers(n);
  f.p_initial = s.numbers(n);
  f.phase_initial = s.numbers(n);
  f.p_reaction = s.numbers(n);
  s.text("Average vibrational energy");
  f.e_vib = s.numbers(n);
  s.text("Average rotational energy");
  f.e_rot = s.numbers(n);
  s.text("Average diffractional energy");
  f.e_dif = s.numbers(n);
  s.text("Summed vibrational probabilities");
  for (int v = 0; v < n_v; ++v) {
    const int iv = s.integer();
    f.p_v[iv] = s.numbers(n);
  }
  s.text("Average diffractional and rotational energies");
  for (int v = 0; v < n_v; ++v) {
    const int iv = s.integer();
    f.e_dif_v[iv] = s.numbers(n);
    f.e_rot_v[iv] = s.numbers(n);
  }
  s.text("Rovibrational excitation probabilities");
  while (s.peek_labels() == 3) {
    const int v = s.integer();
    const int j = s.integer();
    const int mj = s.integer();
    f.p_rovib[{v, j, mj}] = s.numbers(n);
  }
  while (with_orders && s.peek_labels() == 4) {
    const int order = s.integer();
    const int v = s.integer();
    const int j = s.integer();
    const int mj = s.integer();
    f.p_order[{order, v, j, mj}] = s.numbers(n);
  }
  if (!s.at_end()) throw std::runtime_error("unexpected content at the end of " + path.string());
  return f;
}

// Unit 3 (writeSMatrix): the probability of every state on the fine grid,
// (order, n, m, v, j, mJ) -> nen95 values.
inline std::map<std::tuple<int, int, int, int, int, int>, std::vector<std::string>> read_fortran_unit3(
    const std::filesystem::path& path) {
  NumberStream s(path);
  (void)s.numbers(3);
  const int n = s.integer();
  (void)s.numbers(3);
  (void)s.numbers(6);
  (void)s.numbers(n);  // energies
  std::map<std::tuple<int, int, int, int, int, int>, std::vector<std::string>> states;
  while (!s.at_end()) {
    const int order = s.integer();
    const int nq = s.integer();
    const int mq = s.integer();
    const int v = s.integer();
    const int j = s.integer();
    const int mj = s.integer();
    states[{order, nq, mq, v, j, mj}] = s.numbers(n);
  }
  return states;
}

// fort.25 (smatrix_combine): per energy block the S-matrix element, probability and phase of
// every printed state (open rovibrational channels), (v, j, mJ, n, m) -> (Re S, Im S, P, phase).
struct Fort25Entry {
  std::string re, im, p, phase;
};
struct Fort25Block {
  std::string total_energy;
  std::map<std::tuple<int, int, int, int, int>, Fort25Entry> states;
};

inline std::vector<Fort25Block> read_fortran_fort25(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path.string());
  std::vector<Fort25Block> blocks;
  std::string line;
  int v = 0, j = 0, mj = 0;
  while (std::getline(in, line)) {
    const std::vector<std::string> t = split(line);
    if (line.find("total energy =") != std::string::npos) {
      blocks.push_back({t.back(), {}});
    } else if (line.find("v j mj dif.order are:") != std::string::npos) {
      v = std::stoi(t[t.size() - 4]);
      j = std::stoi(t[t.size() - 3]);
      mj = std::stoi(t[t.size() - 2]);
    } else if (t.size() == 6 && is_number(t[0]) && is_number(t[2]) && !blocks.empty()) {
      blocks.back().states[{v, j, mj, std::stoi(t[0]), std::stoi(t[1])}] = {t[2], t[3], t[4], t[5]};
    }
  }
  return blocks;
}

}  // namespace spodvr::testing

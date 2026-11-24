// Readers for the Fortran fine-analysis output (gate 5): the stdout of ReadOrWriteAsympInfo and
// probaf, fil10 (fort.10) and unit 3 (fort.3) (fort.25: fort25_reference.hpp), and the counts gate 5
// must reach. Values are kept as printed text, so they can be compared to half a unit in the last
// printed digit (printed_text_tolerance).
#pragma once

#include <gtest/gtest.h>

#include <cctype>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "analysis_reference.hpp"
#include "fort25_reference.hpp"

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

// The values gate 5 compares in each Fortran output: fil10 (fort.10), unit 3 (fort.3), the probabilities of fort.25
// and its Re S, Im S and phases.
struct FineCounts {
  int f10 = 0;
  int u3 = 0;
  int p25 = 0;
  int s25 = 0;
};

// The counts the Fortran's write rules give for the fine analysis of a run (Analyzor.f90 probaf, smatrix_combine),
// whose labels and interpolated results are in `run` (set_analysis_labels), on nen95 grid energies; s_skipped: the
// fort.25 rows whose S and phase are skipped (with a named reason, counted by the gate).
// fil10 per energy: the energy, P and phase of the initial state, P_reaction, the three average energies, P(v),
// av edif and av erot per v, P(v, j, mJ) and, with writeRDProbs, P(v, j, mJ, order). Unit 3: every state.
// fort.25: the print rule's rows, three values (Re S, Im S, phase) of those not skipped.
inline FineCounts expected_fine_counts(const PropagationRun& run, int nen95, int s_skipped) {
  const int n_rot = static_cast<int>(run.rotation_labels.size());
  const int orders = run.order_block_printed ? run.n_orders : 0;
  const int fort25_rows = fort25_print_rule_rows(run);
  return {nen95 * (7 + 3 * run.n_v + run.n_v * n_rot * (1 + orders)), nen95 * run.n_states * run.n_v, fort25_rows,
          3 * (fort25_rows - s_skipped)};
}

// Gate 5's assertion on its counts: each output compared in full.
inline void expect_fine_counts(const FineCounts& measured, const FineCounts& expected) {
  EXPECT_EQ(measured.f10, expected.f10) << "fil10 values compared vs the Fortran's write rule";
  EXPECT_EQ(measured.u3, expected.u3) << "unit 3 values compared vs the Fortran's write rule";
  EXPECT_EQ(measured.p25, expected.p25) << "fort.25 probabilities compared vs the Fortran's print rule";
  EXPECT_EQ(measured.s25, expected.s25) << "fort.25 S and phases compared vs the rows not skipped";
}

}  // namespace spodvr::testing

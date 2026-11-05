#pragma once
// Readers for the Fortran setup reference files written by
// tools/fortran_reference/run_setup_reference.sh (gates 1 and 2):
//   fortran_setup_stdout.txt  the vasy table ("The gas phase vibrational potential")
//   pes_reference.txt         the potential extracted from the .pes file (extract_pes.py)
//   debug_asyp.dat            rovibrational energies and eigenvectors (DoDump)

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace spodvr::test {

inline double relative_difference(double a, double b) {
  if (a == b) return 0.0;
  return std::abs(a - b) / std::max(std::abs(a), std::abs(b));
}

// Fortran list-directed output may use 'D' exponents.
inline double parse_fortran_real(std::string token) {
  for (char& c : token)
    if (c == 'D' || c == 'd') c = 'E';
  return std::stod(token);
}

// The Fortran "INFO FROM SETTING UP THE ROVIBRATIONAL ASYMPTOTIC POTENTIAL" table: (r, V in eV).
// Empty if the table is not found.
inline std::vector<std::pair<double, double>> read_vasy_table(const std::filesystem::path& stdout_file) {
  std::ifstream file(stdout_file);
  std::string line;
  while (std::getline(file, line) && line.find("The gas phase vibrational potential") == std::string::npos) {
  }
  std::vector<std::pair<double, double>> table;
  if (!file) return table;
  while (std::getline(file, line)) {
    std::istringstream in(line);
    std::string a, b;
    if (!(in >> a >> b)) {
      if (!table.empty()) break;
      continue;
    }
    if (a.find_first_of("0123456789") == std::string::npos) continue;  // header / rule
    table.emplace_back(parse_fortran_real(a), parse_fortran_real(b));
  }
  return table;
}

// pes_reference.txt in the 6D format of extract_pes.py: dimensions, angular grid and a subsample
// of V (hartree) at 0-based grid indices.
struct PesReference {
  int nzpot = 0, nr = 0, nx = 0, ny = 0, ntheta = 0, nphi = 0;
  bool use_symmetry = false;
  int n_xy_computed = 0;
  std::vector<double> thetas;
  std::vector<double> phis;
  struct Point {
    int itheta, iphi, ix, iy, iz, ir;
    double v;
  };
  std::vector<Point> points;
};

inline PesReference read_pes_reference_6d(const std::filesystem::path& path) {
  std::ifstream file(path);
  std::vector<std::string> lines;
  for (std::string line; std::getline(file, line);)
    if (!line.empty() && line[0] != '#') lines.push_back(line);
  PesReference ref;
  std::istringstream header(lines.at(0));
  header >> ref.nzpot >> ref.nr >> ref.nx >> ref.ny >> ref.ntheta >> ref.nphi;
  std::istringstream symmetry(lines.at(1));
  int use_symmetry = 0;
  symmetry >> use_symmetry >> ref.n_xy_computed;
  ref.use_symmetry = use_symmetry != 0;
  std::size_t pos = 2;
  for (int t = 0; t < ref.ntheta; ++t) ref.thetas.push_back(parse_fortran_real(lines.at(pos++)));
  for (int p = 0; p < ref.nphi; ++p) ref.phis.push_back(parse_fortran_real(lines.at(pos++)));
  for (; pos < lines.size(); ++pos) {
    std::istringstream in(lines[pos]);
    PesReference::Point p{};
    std::string v_text;
    in >> p.itheta >> p.iphi >> p.ix >> p.iy >> p.iz >> p.ir >> v_text;
    p.itheta -= 1;
    p.iphi -= 1;
    p.ix -= 1;
    p.iy -= 1;
    p.iz -= 1;
    p.ir -= 1;
    p.v = parse_fortran_real(v_text);
    ref.points.push_back(p);
  }
  return ref;
}

// debug_asyp.dat: records "v : iv j: j E: energy" followed by "Eigenstate: psi(1) ... psi(nr)".
struct RovibrationalReference {
  double energy = 0.0;
  std::vector<double> psi;
};

inline std::map<std::pair<int, int>, RovibrationalReference> read_rovibrational_states(
    const std::filesystem::path& path) {
  std::ifstream file(path);
  std::stringstream contents;
  contents << file.rdbuf();
  std::map<std::pair<int, int>, RovibrationalReference> reference;
  std::string token;
  std::istringstream in(contents.str());
  std::pair<int, int> current{-1, -1};
  bool reading_psi = false;
  while (in >> token) {
    if (token == "v") {
      reading_psi = false;
      std::string colon, j_label, e_label, e_text;
      int v = 0, j = 0;
      in >> colon >> v >> j_label >> j >> e_label >> e_text;
      current = {v, j};
      reference[current].energy = parse_fortran_real(e_text);
    } else if (token == "Eigenstate:") {
      reading_psi = true;
    } else if (token.rfind("===", 0) == 0) {
      reading_psi = false;
    } else if (reading_psi && current.first >= 0) {
      reference[current].psi.push_back(parse_fortran_real(token));
    }
  }
  return reference;
}

}  // namespace spodvr::test

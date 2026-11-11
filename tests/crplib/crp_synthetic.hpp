#pragma once
// Synthetic CRP data sets for the tests (no real data set is needed): the cut geometry of each
// symmetry (Appendix A of docs/crp_library.md, written down independently of the library code),
// an exactly symmetric pairwise model PES, and "impulse" data sets (one cut 1, everything else 0).

#include <array>
#include <cmath>
#include <format>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "crplib/dataset.hpp"

namespace crplib_test {

using crplib::DataSet;
using crplib::Symmetry;

struct Vec3 {
  double x, y, z;
};

// One molecular cut: site in lattice units (skewed (u, v) for C3v/C6v, Cartesian for C4v) and the
// orientation in degrees, in the order the Fortran interpolation expects the cuts.
struct CutGeometry {
  double u, v, theta, phi;
};

inline std::vector<CutGeometry> molecule_cut_geometry(Symmetry symmetry) {
  constexpr double third = 1.0 / 3.0, sixth = 1.0 / 6.0;
  const std::vector<CutGeometry> top_bridge = {{0, 0, 0, 0},    {0, 0, 90, 0},    {0, 0, 90, 30},  {0.5, 0, 0, 0},
                                               {0.5, 0, 90, 0}, {0.5, 0, 90, 60}, {0.5, 0, 90, 90}};
  std::vector<CutGeometry> hollow = {{third, third, 0, 0},
                                     {third, third, 90, 0},
                                     {third, third, 90, 30},
                                     {third, third, 45, 30},
                                     {third, third, 45, 210}};
  std::vector<CutGeometry> cuts = top_bridge;
  if (symmetry == Symmetry::c4v) {
    return {{0, 0, 0, 0},         {0, 0, 90, 0},        {0, 0, 90, 45},        {0.5, 0, 0, 0},
            {0.5, 0, 90, 0},      {0.5, 0, 90, 45},     {0.5, 0, 90, 90},      {0.5, 0.5, 0, 0},
            {0.5, 0.5, 90, 0},    {0.5, 0.5, 90, 45},   {0.25, 0.25, 0, 0},    {0.25, 0.25, 45, 225},
            {0.25, 0.25, 45, 45}, {0.25, 0.25, 90, 45}, {0.25, 0.25, 90, 135}, {0.25, 0.25, 45, 135}};
  }
  cuts.insert(cuts.end(), hollow.begin(), hollow.end());  // C6v hollow / C3v hcp
  if (symmetry == Symmetry::c6v) {
    for (const auto& [theta, phi] : {std::pair{0, 0}, {90, 30}, {90, 120}, {45, 30}, {45, 210}, {45, 120}})
      cuts.push_back({sixth, sixth, double(theta), double(phi)});  // t2h (InterpolateOnT2F of C6v)
    return cuts;
  }
  for (const auto& [theta, phi] : {std::pair{0, 0}, {90, 30}, {90, 120}, {45, 30}, {45, 120}, {45, 210}})
    cuts.push_back({sixth, sixth, double(theta), double(phi)});  // C3v t2h
  for (const auto& [theta, phi] : {std::pair{0, 0}, {90, 240}, {90, 330}, {45, 150}, {45, 240}, {45, 330}})
    cuts.push_back({third, 5.0 * sixth, double(theta), double(phi)});  // C3v t2f
  for (CutGeometry h : hollow) {
    h.u = 2 * third;  // C3v fcc hollow (2/3, 2/3)
    h.v = 2 * third;
    cuts.push_back(h);
  }
  return cuts;
}

// Atom-surface sites, Cartesian, in lattice units (as CRPData.3D), in the Fortran order.
inline std::vector<std::pair<double, double>> atom_sites(Symmetry symmetry) {
  constexpr double third = 1.0 / 3.0, sixth = 1.0 / 6.0;
  const std::vector<std::pair<double, double>> skewed_c6v = {{0, 0},         {0, 0.5},   {third, third}, {sixth, sixth},
                                                             {sixth, third}, {0, sixth}, {0, 2 * third}};
  const std::vector<std::pair<double, double>> skewed_c3v = {
      {0, 0},     {0.5, 0.5}, {third, third},     {sixth, 2 * third}, {sixth, 0.5},
      {0, sixth}, {0, third}, {third, 5 * sixth}, {0.5, 2 * third},   {2 * third, 2 * third}};
  if (symmetry == Symmetry::c4v) return {{0, 0}, {0.5, 0.5}, {0.5, 0}, {0.25, 0.25}, {0.5, 0.25}, {0.25, 0}};
  std::vector<std::pair<double, double>> cartesian;
  for (const auto& [u, v] : symmetry == Symmetry::c6v ? skewed_c6v : skewed_c3v)
    cartesian.emplace_back(u + v / 2.0, v * std::sqrt(3.0) / 2.0);
  return cartesian;
}

// Cartesian position (bohr) of a molecular cut's site.
inline std::pair<double, double> cut_site_cartesian(Symmetry symmetry, const CutGeometry& c, double a) {
  if (symmetry == Symmetry::c4v) return {c.u * a, c.v * a};
  return {(c.u + c.v / 2.0) * a, c.v * std::sqrt(3.0) / 2.0 * a};
}

// The two atom positions of a molecule with centre (x, y, z), bond length r and orientation (theta,
// phi) in radians, as the Fortran builds them (rab).
inline std::pair<Vec3, Vec3> atom_pair(double x, double y, double z, double r, double theta, double phi) {
  const double cz = 0.5 * r * std::cos(theta), cx = 0.5 * r * std::sin(theta) * std::cos(phi),
               cy = 0.5 * r * std::sin(theta) * std::sin(phi);
  return {{x - cx, y - cy, z - cz}, {x + cx, y + cy, z + cz}};
}

// Smooth cut-off: 1 below `on`, 0 above `off`, C2 in between.
inline double switch_off(double rho, double on, double off) {
  if (rho <= on) return 1.0;
  if (rho >= off) return 0.0;
  const double t = (rho - on) / (off - on);
  return (1 - t) * (1 - t) * (1 - t) * (1 + 3 * t + 6 * t * t);
}
inline double switch_off_slope(double rho, double on, double off) {
  if (rho <= on || rho >= off) return 0.0;
  const double t = (rho - on) / (off - on);
  return -30.0 * t * t * (1 - t) * (1 - t) / (off - on);
}

// An exactly symmetric model PES: H-surface atom Morse pairs summed over a slab with the lattice of
// the symmetry (C6v: one triangular layer; C3v: four ABC-stacked layers; C4v: two square layers),
// an H-H Morse term and a centre-of-mass coupling. Hartree, bohr.
class PairwiseModel {
 public:
  PairwiseModel(Symmetry symmetry, double a, std::vector<double> spacings)
      : symmetry_(symmetry), a_(a), spacings_(std::move(spacings)) {}

  // Pair potential and its slope (Morse, cut off smoothly at 8 bohr).
  static double pair(double rho) {
    const double e = std::exp(-1.0 * (rho - 3.2));
    return 0.01 * (e * e - 2 * e) * switch_off(rho, 6.0, 8.0);
  }
  static double pair_slope(double rho) {
    const double e = std::exp(-1.0 * (rho - 3.2));
    return 0.01 * (-2 * e * e + 2 * e) * switch_off(rho, 6.0, 8.0) +
           0.01 * (e * e - 2 * e) * switch_off_slope(rho, 6.0, 8.0);
  }
  static double coupling(double rho) { return 0.004 * std::exp(-0.9 * rho) * switch_off(rho, 6.0, 8.0); }
  // H2 Morse potential (free molecule).
  static double h2(double r) {
    const double e = std::exp(-1.028 * (r - 1.40));
    return 0.1745 * (e * e - 2 * e);
  }
  // Gas-phase Z dependence (only enters through extrapol.dat).
  static double z_dependence(double z) { return -0.001 * std::exp(-(z - 4.0)); }
  // The CRP corrugation function: deliberately not the pair potential (so the interpolated
  // remainder is not trivial); zero with zero slope at 6 bohr.
  static double corrugation(double rho) { return 0.03 * std::exp(-1.2 * rho) * switch_off(rho, 4.0, 6.0); }
  static double corrugation_slope(double rho) {
    return 0.03 * (-1.2 * std::exp(-1.2 * rho) * switch_off(rho, 4.0, 6.0) +
                   std::exp(-1.2 * rho) * switch_off_slope(rho, 4.0, 6.0));
  }

  // Sum of f(|p - R|) and of its z-gradient over the surface atoms R within 8 bohr.
  template <typename F>
  double sum_over_atoms(const Vec3& p, F f) const {
    double sum = 0.0;
    for_each_atom(p, [&](double dx, double dy, double dz) { sum += f(std::sqrt(dx * dx + dy * dy + dz * dz)); });
    return sum;
  }
  double atom(const Vec3& p) const { return sum_over_atoms(p, pair); }
  double atom_dz(const Vec3& p) const {
    double sum = 0.0;
    for_each_atom(p, [&](double dx, double dy, double dz) {
      const double rho = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (rho > 0.0) sum += pair_slope(rho) * dz / rho;
    });
    return sum;
  }
  double molecule(const Vec3& a, const Vec3& b) const {
    const Vec3 centre{(a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2};
    const double r = std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z - a.z) * (b.z - a.z));
    return atom(a) + atom(b) + h2(r) + (1.0 + 0.5 * (r - 1.4)) * sum_over_atoms(centre, coupling);
  }

 private:
  // Calls f(p - R) for every atom R of the slab within reach (8 bohr) of p.
  template <typename F>
  void for_each_atom(const Vec3& p, F f) const {
    const bool square = symmetry_ == Symmetry::c4v;
    const int layers = symmetry_ == Symmetry::c6v ? 1 : static_cast<int>(spacings_.size()) + 1;
    const int reach = static_cast<int>(8.0 / (a_ * (square ? 1.0 : std::sqrt(3.0) / 2.0))) + 3;
    // crystal coordinates of p
    const double v0 = square ? p.y / a_ : p.y / (a_ * std::sqrt(3.0) / 2.0);
    const double u0 = square ? p.x / a_ : p.x / a_ - v0 / 2.0;
    double layer_z = 0.0;
    for (int k = 0; k < layers; ++k) {
      if (k > 0) layer_z -= spacings_[k - 1];
      if (layer_z < p.z - 8.0) break;
      // layer shift: ABC stacking (k mod 3)/3 for C3v, (1/2, 1/2) for the second C4v layer
      const double shift = square ? (k % 2) * 0.5 : (k % 3) / 3.0;
      for (int i = static_cast<int>(std::floor(u0)) - reach; i <= static_cast<int>(std::floor(u0)) + reach; ++i) {
        for (int j = static_cast<int>(std::floor(v0)) - reach; j <= static_cast<int>(std::floor(v0)) + reach; ++j) {
          const double u = i + shift, v = j + shift;
          const double x = square ? u * a_ : (u + v / 2.0) * a_;
          const double y = square ? v * a_ : v * std::sqrt(3.0) / 2.0 * a_;
          const double dx = p.x - x, dy = p.y - y, dz = p.z - layer_z;
          if (dx * dx + dy * dy + dz * dz < 64.0) f(dx, dy, dz);
        }
      }
    }
  }

  Symmetry symmetry_;
  double a_;
  std::vector<double> spacings_;
};

// Grids of the synthetic data sets (bohr). The atom z range covers every atom position of the
// molecular cut knots: [Z_min - r_max / 2, Z_max + r_max / 2].
struct SyntheticGrids {
  std::vector<double> r, z, atom_z, corrugation_rho, z_dependence_z, h2_r;
  double zso_min = 5.0, zso_max = 6.5;
};

inline std::vector<double> grid(double first, double last, int n, double wiggle = 0.0) {
  std::vector<double> g(n);
  for (int i = 0; i < n; ++i) {
    const double t = static_cast<double>(i) / (n - 1);
    g[i] = first + (last - first) * t + (i > 0 && i < n - 1 ? wiggle * (last - first) / n * std::sin(7.0 * i) : 0.0);
  }
  return g;
}

inline SyntheticGrids default_grids() {
  SyntheticGrids g;
  g.r = grid(0.6, 3.2, 14, 0.2);
  g.z = grid(0.5, 7.0, 16, 0.2);
  g.atom_z = grid(-1.6, 9.0, 34, 0.2);
  g.corrugation_rho = grid(0.0, 6.0, 25);
  g.z_dependence_z = grid(4.0, 10.0, 13);
  g.h2_r = grid(0.4, 8.0, 39);
  return g;
}

// Small grids for the impulse data sets (their response does not depend on the grid); the atom z
// range is wide so that the (zero) atom PES never stops.
inline SyntheticGrids impulse_grids() {
  SyntheticGrids g;
  g.r = grid(0.6, 3.2, 5);
  g.z = grid(0.5, 7.0, 5);
  g.atom_z = grid(-10.0, 20.0, 5);
  g.corrugation_rho = grid(0.0, 6.0, 5);
  g.z_dependence_z = grid(4.0, 10.0, 5);
  g.h2_r = grid(0.4, 8.0, 5);
  return g;
}

inline double lattice_constant(Symmetry symmetry) { return symmetry == Symmetry::c6v ? 4.8 : 4.82; }
inline std::vector<double> spacings(Symmetry symmetry) {
  if (symmetry == Symmetry::c3v) return {3.94, 3.94, 3.94};
  if (symmetry == Symmetry::c4v) return {3.41};
  return {3.94};
}

// The parts every synthetic data set shares: CRPData, cut lists (names, geometry), grids.
inline DataSet skeleton(Symmetry symmetry, const SyntheticGrids& g) {
  DataSet d;
  d.symmetry = symmetry;
  d.units = crplib::internal_units;
  d.lattice_constant = lattice_constant(symmetry);
  d.interlayer_spacings = spacings(symmetry);
  d.info = crplib::DataSetInfo{"synthetic test data set", std::string(crplib::symmetry_name(symmetry)), "none"};
  d.zso_min = g.zso_min;
  d.zso_max = g.zso_max;
  d.corrugation.x = g.corrugation_rho;
  d.z_dependence.x = g.z_dependence_z;
  d.h2_potential.x = g.h2_r;
  int k = 0;
  for (const auto& [x, y] : atom_sites(symmetry)) {
    crplib::AtomCut cut;
    cut.x = x;
    cut.y = y;
    cut.file_name = std::format("site{}.dat", ++k);
    cut.potential.x = g.atom_z;
    d.atom_cuts.push_back(cut);
  }
  k = 0;
  for (const CutGeometry& c : molecule_cut_geometry(symmetry)) {
    crplib::MoleculeCut cut;
    cut.x = c.u;
    cut.y = c.v;
    cut.theta = c.theta;
    cut.phi = c.phi;
    cut.file_name = std::format("cut{:02d}.dat", ++k);
    cut.r = g.r;
    cut.z = g.z;
    d.molecule_cuts.push_back(cut);
  }
  return d;
}

// Every potential 0, molecular cut `impulse_cut` (0-based) equal to 1, CutsInterpolation true:
// the PES is then the angular and Fourier response to that one cut (VH2 last value 0).
inline DataSet impulse_data_set(Symmetry symmetry, int impulse_cut, const SyntheticGrids& g = default_grids()) {
  DataSet d = skeleton(symmetry, g);
  d.cuts_interpolation = true;
  d.corrugation.f.assign(d.corrugation.x.size(), 0.0);
  d.z_dependence.f.assign(d.z_dependence.x.size(), 0.0);
  d.h2_potential.f.assign(d.h2_potential.x.size(), 0.0);
  for (crplib::AtomCut& cut : d.atom_cuts) cut.potential.f.assign(cut.potential.x.size(), 0.0);
  for (std::size_t k = 0; k < d.molecule_cuts.size(); ++k) {
    crplib::MoleculeCut& cut = d.molecule_cuts[k];
    cut.v.assign(cut.r.size() * cut.z.size(), static_cast<int>(k) == impulse_cut ? 1.0 : 0.0);
  }
  return d;
}

// The pairwise model tabulated as a CutsInterpolation = false data set (total energies), in atomic
// units.
inline DataSet pairwise_data_set(Symmetry symmetry, const SyntheticGrids& g = default_grids()) {
  DataSet d = skeleton(symmetry, g);
  const PairwiseModel model(symmetry, d.lattice_constant, d.interlayer_spacings);
  d.cuts_interpolation = false;
  for (double rho : d.corrugation.x) d.corrugation.f.push_back(PairwiseModel::corrugation(rho));
  d.corrugation_slope_first = PairwiseModel::corrugation_slope(d.corrugation.x.front());
  d.corrugation_slope_last = PairwiseModel::corrugation_slope(d.corrugation.x.back());
  for (double z : d.z_dependence.x) d.z_dependence.f.push_back(PairwiseModel::z_dependence(z));
  for (double r : d.h2_potential.x) d.h2_potential.f.push_back(PairwiseModel::h2(r));
  for (crplib::AtomCut& cut : d.atom_cuts) {
    const double x = cut.x * d.lattice_constant, y = cut.y * d.lattice_constant;
    for (double z : cut.potential.x) cut.potential.f.push_back(model.atom({x, y, z}));
    cut.slope_first = model.atom_dz({x, y, cut.potential.x.front()});
    cut.slope_last = model.atom_dz({x, y, cut.potential.x.back()});
  }
  const auto geometry = molecule_cut_geometry(symmetry);
  for (std::size_t k = 0; k < d.molecule_cuts.size(); ++k) {
    crplib::MoleculeCut& cut = d.molecule_cuts[k];
    const auto [x, y] = cut_site_cartesian(symmetry, geometry[k], d.lattice_constant);
    for (double z : cut.z) {
      for (double r : cut.r) {
        const auto [a, b] = atom_pair(x, y, z, r, geometry[k].theta * M_PI / 180.0, geometry[k].phi * M_PI / 180.0);
        cut.v.push_back(model.molecule(a, b));
      }
    }
  }
  return d;
}

// The same data set with every length in Angstrom and every energy in eV (the CRP library's
// factors), as a data set written in those units would be.
inline DataSet in_angstrom_and_ev(DataSet d) {
  constexpr double bohr = 0.529177249, ev = 27.2113961;
  auto lengths = [](std::vector<double>& v) {
    for (double& x : v) x *= bohr;
  };
  auto energies = [](std::vector<double>& v) {
    for (double& x : v) x *= ev;
  };
  d.units = crplib::Units{crplib::LengthUnit::angstrom, crplib::EnergyUnit::ev};
  d.lattice_constant *= bohr;
  lengths(d.interlayer_spacings);
  lengths(d.corrugation.x);
  energies(d.corrugation.f);
  d.corrugation_slope_first *= ev / bohr;
  d.corrugation_slope_last *= ev / bohr;
  for (crplib::AtomCut& cut : d.atom_cuts) {
    lengths(cut.potential.x);
    energies(cut.potential.f);
    cut.slope_first *= ev / bohr;
    cut.slope_last *= ev / bohr;
  }
  d.zso_min *= bohr;
  d.zso_max *= bohr;
  lengths(d.z_dependence.x);
  energies(d.z_dependence.f);
  lengths(d.h2_potential.x);
  energies(d.h2_potential.f);
  for (crplib::MoleculeCut& cut : d.molecule_cuts) {
    lengths(cut.r);
    lengths(cut.z);
    energies(cut.v);
  }
  return d;
}

// The symmetry of a reference case (tests/crplib/reference/<case>/): from its name prefix.
inline Symmetry case_symmetry(const std::string& name) {
  if (name.starts_with("c3v")) return Symmetry::c3v;
  if (name.starts_with("c4v")) return Symmetry::c4v;
  return Symmetry::c6v;
}

// The C3v pairwise data set with two extra cuts after the 29 used ones: a copy of cut 1 on a shorter
// r grid (it narrows Rmax) and one without a file (skipped, as in the Fortran).
inline DataSet with_missing_extra_cut() {
  DataSet d = pairwise_data_set(Symmetry::c3v);
  crplib::MoleculeCut narrow = d.molecule_cuts[0];
  narrow.file_name = "cut30.dat";
  const std::size_t nr = narrow.r.size(), kept = nr - 4;
  std::vector<double> v;
  for (std::size_t j = 0; j < narrow.z.size(); ++j)
    for (std::size_t i = 0; i < kept; ++i) v.push_back(narrow.v[i + nr * j]);
  narrow.r.resize(kept);
  narrow.v = v;
  crplib::MoleculeCut absent = d.molecule_cuts[1];
  absent.file_name = "cut31.dat";
  absent.present = false;
  absent.r.clear();
  absent.z.clear();
  absent.v.clear();
  d.molecule_cuts.push_back(narrow);
  d.molecule_cuts.push_back(absent);
  return d;
}

// The data set of a reference case: <sym>_impulse_<k> (cut k = 1), <sym>_pairwise,
// c6v_pairwise_angstrom, c4v_negative_spacing and c3v_missing_extra_cut. Throws
// std::invalid_argument for other names.
inline DataSet case_data_set(const std::string& name) {
  const Symmetry symmetry = case_symmetry(name);
  if (name.find("_impulse_") != std::string::npos)
    return impulse_data_set(symmetry, std::stoi(name.substr(name.rfind('_') + 1)) - 1, impulse_grids());
  if (name == "c6v_pairwise_angstrom") return in_angstrom_and_ev(pairwise_data_set(symmetry));
  if (name == "c4v_negative_spacing") {
    DataSet d = pairwise_data_set(symmetry);
    d.interlayer_spacings = {-0.98};  // relative: z of the second layer = 0.98 a / sqrt(2)
    return d;
  }
  if (name.ends_with("_pairwise")) return pairwise_data_set(symmetry);
  if (name == "c3v_missing_extra_cut") return with_missing_extra_cut();
  throw std::invalid_argument("unknown reference case " + name);
}

}  // namespace crplib_test

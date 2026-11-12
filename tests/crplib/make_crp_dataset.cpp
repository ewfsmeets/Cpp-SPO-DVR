// crplib_make_dataset <case> <case_dir>: writes the synthetic data set of a reference case
// (<case_dir>/dataset/) and the point lists the Fortran reference driver evaluates for it
// (<case_dir>/<mode>[_<length>_<energy>].in, hexadecimal doubles). Used by
// tools/fortran_reference/run_crp_reference.sh --generate; the tests read the same files.
//
// Cases: <symmetry>_pairwise (c6v, c3v, c4v), c6v_pairwise_angstrom (the data set in Angstrom/eV),
// c4v_negative_spacing, <symmetry>_impulse_<k> (cut k, 1-based, equal to 1, all else 0).
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>
#include <vector>

#include "crp_synthetic.hpp"
#include "crplib/dataset.hpp"
#include "crplib/molecule_surface.hpp"
#include "crplib/units.hpp"

namespace {

using crplib::Symmetry;
using Point = std::vector<double>;

std::string hex(double x) { return std::format("{:016x}", std::bit_cast<std::uint64_t>(x)); }

void write_points(const std::filesystem::path& file, const std::vector<Point>& points) {
  std::ofstream out(file);
  out << points.size() << '\n';
  for (const Point& p : points) {
    for (std::size_t i = 0; i < p.size(); ++i) out << (i ? " " : "") << hex(p[i]);
    out << '\n';
  }
}

// Cartesian (bohr) of the lattice coordinates (u, v): skewed for C3v/C6v, Cartesian for C4v.
std::pair<double, double> cartesian(Symmetry symmetry, double u, double v, double a) {
  if (symmetry == Symmetry::c4v) return {u * a, v * a};
  return {(u + v / 2.0) * a, v * std::sqrt(3.0) / 2.0 * a};
}

// A deterministic sparse subset of a product grid: keeps about one in `stride` combinations.
bool keep(int index, int stride) { return (index * 7919 + 13) % stride == 0; }

// Atom positions: lattice points (incl. negative ones, the ORIGEN quirk), points outside the cell,
// fold boundaries, and z at the ends and inside the data range (first/last interval: end slopes).
std::vector<Point> atom_points(Symmetry symmetry, double a, const std::vector<double>& z) {
  const double z_first = z.front(), z_last = z.back();
  const std::vector<double> heights = {z_first, z_first + 0.05, 0.37, 1.9, 4.4, z_last - 0.05, z_last};
  const std::vector<double> coordinates = {-2.0, -1.0, -0.37, 0.0, 0.21, 0.25, 0.5, 0.75, 0.9, 1.0, 1.73};
  std::vector<Point> points;
  for (double u : coordinates) {
    for (double v : coordinates) {
      for (double h : heights) {
        const int key = static_cast<int>(std::lround(100 * u + 37 * v + 11 * h));
        if ((key % 3 + 3) % 3 != 0 && !(u == v)) continue;
        const auto [x, y] = cartesian(symmetry, u, v, a);
        points.push_back({x, y, h});
      }
    }
  }
  if (symmetry == Symmetry::c4v) {  // REGION boundaries: x = y, x = a/2, y = a/2, x + y = a
    for (const auto& [x, y] : {std::pair{0.3 * a, 0.3 * a},
                               {0.5 * a, 0.2 * a},
                               {0.2 * a, 0.5 * a},
                               {0.7 * a, 0.3 * a},
                               {0.5 * a, 0.5 * a},
                               {0.8 * a, 0.8 * a},
                               {0.6 * a, 0.9 * a},
                               {0.9 * a, 0.6 * a}}) {
      for (double h : {0.37, 2.2}) points.push_back({x, y, h});
    }
  }
  return points;
}

// The two atoms (driver order Za Xa Ya Zb Xb Yb) of a molecule with centre (x, y, z), bond length r
// and orientation (theta, phi), as the Fortran wrapper builds them.
Point molecule(double x, double y, double z, double r, double theta, double phi) {
  const double cz = 0.5 * r * std::cos(theta), cx = 0.5 * r * std::sin(theta) * std::cos(phi),
               cy = 0.5 * r * std::sin(theta) * std::sin(phi);
  return {z - cz, x - cx, y - cy, z + cz, x + cx, y + cy};
}

// Molecules in every branch: below ZsoMin, in the switching region, at and above ZsoMax, beyond
// Zasy; r and Z at and outside the cut grids (the silent V = 0 branch); theta = 0 (rho = 0);
// centres inside, outside and on the edges of the cell. Every atom of an evaluated near-surface
// point stays inside the atom z range (no Fortran STOP).
std::vector<Point> molecule_points(Symmetry symmetry, double a) {
  const std::vector<std::pair<double, double>> centres = {{0, 0},       {0.13, 0.29}, {0.5, 0},   {0.41, 0.07},
                                                          {0.71, 0.66}, {-0.37, 1.2}, {1.0, 1.0}, {2.3, -0.6}};
  const std::vector<double> heights = {0.4, 0.5, 0.8, 1.7, 3.1, 4.9, 5.0, 5.6, 6.4, 6.5, 6.9, 8.5, 11.0};
  const std::vector<double> bonds = {0.3, 0.55, 0.6, 0.95, 1.4, 2.2, 3.1, 3.2, 3.4};
  const std::vector<double> thetas = {0.0, 0.31, 0.9, std::numbers::pi / 2, 2.3};
  const std::vector<double> phis = {0.0, 0.7, 2.9, -1.3};
  std::vector<Point> points;
  int index = 0;
  for (const auto& [u, v] : centres) {
    const auto [x, y] = cartesian(symmetry, u, v, a);
    for (double z : heights)
      for (double r : bonds)
        for (double theta : thetas)
          for (double phi : phis)
            if (keep(index++, 23)) points.push_back(molecule(x, y, z, r, theta, phi));
  }
  return points;
}

// The QD wrapper's arguments u, v (bohr along the lattice vectors), Z, r, theta, phi, including
// points that are clamped in r (both ends) and Z (below Zmin).
std::vector<Point> wrapper_points(double a) {
  // u = 0, v != 0: x = cos(pi/3) v alone, so the folded skew constant shows up in V (test of the
  // mutation 0.5 for 0x1.0000000000001p-1; with u != 0 its ulp is rounded away).
  const std::vector<std::pair<double, double>> centres = {{0, 0},      {0.13, 0.29}, {0.5, 0.5},  {0.83, 0.41},
                                                          {1.7, -0.4}, {0.0, 0.37},  {0.0, -0.62}};
  const std::vector<double> heights = {-2.0, 0.5, 2.0, 5.5, 6.7, 9.0, 12.0};
  const std::vector<double> bonds = {0.3, 0.8, 1.4, 2.6, 4.0, 9.0};
  const std::vector<double> thetas = {0.05, 0.8, std::numbers::pi / 2};
  const std::vector<double> phis = {0.0, 1.1, 3.0};
  std::vector<Point> points;
  int index = 0;
  for (const auto& [u, v] : centres)
    for (double z : heights)
      for (double r : bonds)
        for (double theta : thetas)
          for (double phi : phis)
            if (keep(index++, 7)) points.push_back({u * a, v * a, z, r, theta, phi});
  return points;
}

// (Zcc, r) of the gas-phase routines (Zcc >= ZsoMin: below it the Fortran stops).
std::vector<Point> gas_phase_points() {
  std::vector<Point> points;
  for (double z : {5.0, 5.3, 6.0, 6.5, 8.0, 10.0, 12.0})
    for (double r : {0.5, 1.4, 3.0, 7.0}) points.push_back({z, r});
  return points;
}

// Atom pairs for CAMBCORD: general, vertical (dx = dy = 0), antiparallel.
std::vector<Point> cambcord_points() {
  return {{1.0, 0.3, 0.2, 2.1, 1.1, -0.4}, {1.0, 0.5, 0.5, 2.4, 0.5, 0.5}, {3.0, 2.0, 1.0, 1.0, 2.0, 1.0},
          {0.7, -1.2, 3.3, 0.9, 0.4, 2.2}, {2.5, 0.0, 0.0, 2.5, 1.4, 0.0}, {5.1, 7.3, -2.2, 4.4, 6.9, -2.9}};
}

std::vector<Point> scaled(std::vector<Point> points, double factor) {
  for (Point& p : points)
    for (double& x : p) x *= factor;
  return points;
}

// Molecules for a real data set, placed by its own ranges (bohr): below Zmin, in the cut grids,
// at ZsoMin, in the switching region, at ZsoMax, beyond Zasy; r below Rmin, inside, beyond Rmax.
// An evaluated molecule has Zmin <= Zcc <= ZsoMax and Rmin <= r <= Rmax, so its atoms stay within
// the atom cuts when those cover [Zmin - Rmax/2, ZsoMax + Rmax/2] (no Fortran STOP).
std::vector<Point> real_molecule_points(Symmetry symmetry, double a, const crplib::MoleculeSurface& pes) {
  const double z_min = pes.z_min(), zso_min = pes.zso_min(), zso_max = pes.zso_max(), zasy = pes.z_asymptotic();
  const double r_min = pes.r_min(), r_max = pes.r_max();
  const std::vector<std::pair<double, double>> centres = {{0, 0},       {0.13, 0.29}, {0.5, 0},   {0.41, 0.07},
                                                          {0.71, 0.66}, {-0.37, 1.2}, {1.0, 1.0}, {2.3, -0.6}};
  const std::vector<double> heights = {z_min - 0.2,
                                       z_min + 0.1,
                                       z_min + 0.6,
                                       0.5 * (z_min + zso_min),
                                       zso_min - 0.1,
                                       zso_min,
                                       0.5 * (zso_min + zso_max),
                                       zso_max,
                                       zso_max + 0.5,
                                       zasy,
                                       zasy + 2.0};
  const std::vector<double> bonds = {r_min - 0.1, r_min + 0.05, 1.4, 0.5 * (r_min + r_max), r_max - 0.1, r_max + 0.2};
  const std::vector<double> thetas = {0.0, 0.31, 0.9, std::numbers::pi / 2, 2.3};
  const std::vector<double> phis = {0.0, 0.7, 2.9, -1.3};
  std::vector<Point> points;
  int index = 0;
  for (const auto& [u, v] : centres) {
    const auto [x, y] = cartesian(symmetry, u, v, a);
    for (double z : heights)
      for (double r : bonds)
        for (double theta : thetas)
          for (double phi : phis)
            if (keep(index++, 19)) points.push_back(molecule(x, y, z, r, theta, phi));
  }
  return points;
}

// (Zcc, r) of the gas-phase routines for a real data set: Zcc >= ZsoMin (below it the Fortran
// stops), r inside the free-molecule table (bohr).
std::vector<Point> real_gas_phase_points(const crplib::MoleculeSurface& pes, double r_first, double r_last) {
  std::vector<Point> points;
  const double zso_min = pes.zso_min(), zso_max = pes.zso_max(), zasy = pes.z_asymptotic();
  for (double z : {zso_min + 0.05, zso_min + 0.3, 0.5 * (zso_min + zso_max), zso_max, zasy - 0.5, zasy, zasy + 2.0})
    for (double r : {r_first + 0.1, 1.4, 3.0, r_last - 0.1}) points.push_back({z, r});
  return points;
}

// The point lists of a reference case on a real data set (which stays where it is): atoms over the
// atom z range, molecules and gas-phase points by the data set's ranges, the QD wrapper, CAMBCORD,
// and the atom and molecule points again in Angstrom (the data set's own length unit).
void write_real_case(const std::filesystem::path& dir, const std::filesystem::path& data_set_dir) {
  const crplib::DataSet d = crplib::read_data_set(data_set_dir);
  const crplib::MoleculeSurface pes(d);
  const crplib::LengthUnit unit = d.units.length;
  const double a = crplib::convert_length(d.lattice_constant, unit, crplib::LengthUnit::bohr);
  std::vector<double> atom_z;
  for (double z : d.atom_cuts.front().potential.x)
    atom_z.push_back(crplib::convert_length(z, unit, crplib::LengthUnit::bohr));
  const double r_first = crplib::convert_length(d.h2_potential.x.front(), unit, crplib::LengthUnit::bohr);
  const double r_last = crplib::convert_length(d.h2_potential.x.back(), unit, crplib::LengthUnit::bohr);
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  write_points(dir / "atom_0_2.in", atom_points(d.symmetry, a, atom_z));
  write_points(dir / "molecule_0_2.in", real_molecule_points(d.symmetry, a, pes));
  write_points(dir / "wrapper.in", wrapper_points(a));
  write_points(dir / "cambcord.in", cambcord_points());
  write_points(dir / "gasphase.in", real_gas_phase_points(pes, r_first, r_last));
  constexpr double bohr = 0.529177249;  // the same atom and molecule points, given in Angstrom
  write_points(dir / "molecule_1_3.in", scaled(real_molecule_points(d.symmetry, a, pes), bohr));
  write_points(dir / "atom_1_3.in", scaled(atom_points(d.symmetry, a, atom_z), bohr));
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3 && argc != 4) {
    std::cerr << "usage: crplib_make_dataset <case> <case_dir> [<real data set>]\n"
                 "  with a real data set: only the point lists (the data set stays where it is)\n";
    return 1;
  }
  try {
    using namespace crplib_test;
    const std::string name = argv[1];
    const std::filesystem::path dir = argv[2];
    if (argc == 4) {
      write_real_case(dir, argv[3]);
      std::cout << "wrote " << dir.string() << "\n";
      return 0;
    }
    const Symmetry symmetry = case_symmetry(name);
    const double a = lattice_constant(symmetry);  // bohr
    const crplib::DataSet d = case_data_set(name);
    const bool impulse = name.find("_impulse_") != std::string::npos;
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    crplib::write_data_set(d, dir / "dataset");

    if (impulse) {  // pure interpolation: molecules below ZsoMin, any orientation
      std::vector<Point> points;
      for (const auto& [u, v] : {std::pair{0.0, 0.0},
                                 {0.5, 0.0},
                                 {1.0 / 3, 1.0 / 3},
                                 {1.0 / 6, 1.0 / 6},
                                 {0.5, 0.5},
                                 {0.0, 0.5},
                                 {2.0 / 3, 1.0 / 6},
                                 {1.0 / 6, 2.0 / 3},
                                 {2.0 / 3, 2.0 / 3},
                                 {5.0 / 6, 5.0 / 6},
                                 {1.0 / 3, 5.0 / 6},
                                 {5.0 / 6, 1.0 / 3},
                                 {0.25, 0.25},
                                 {0.75, 0.25},
                                 {0.75, 0.75},
                                 {0.25, 0.75},
                                 {0.37, 0.61}}) {
        const auto [x, y] = cartesian(symmetry, u, v, a);
        for (double theta : {0.0, 0.4, 1.2, std::numbers::pi / 2})
          for (double phi : {0.0, 0.6, 2.2, -1.9}) points.push_back(molecule(x, y, 2.0, 1.4, theta, phi));
      }
      write_points(dir / "molecule_0_2.in", points);
      write_points(dir / "wrapper.in", wrapper_points(a));  // compared bitwise (CutsInterpolation = T)
    } else {
      // The atom z range as the library reads it (in the Angstrom case the unit round trip moves
      // the ends by an ulp; a point just outside would stop the Fortran).
      const crplib::DataSet read_back = crplib::read_data_set(dir / "dataset");
      std::vector<double> atom_z;
      for (double z : read_back.atom_cuts.front().potential.x)
        atom_z.push_back(crplib::convert_length(z, read_back.units.length, crplib::LengthUnit::bohr));
      write_points(dir / "atom_0_2.in", atom_points(symmetry, a, atom_z));
      write_points(dir / "molecule_0_2.in", molecule_points(symmetry, a));
      write_points(dir / "wrapper.in", wrapper_points(a));
      write_points(dir / "cambcord.in", cambcord_points());
      if (name == "c6v_pairwise_angstrom") {
        constexpr double bohr = 0.529177249;  // the same points, given in Angstrom
        write_points(dir / "molecule_1_3.in", scaled(molecule_points(symmetry, a), bohr));
        write_points(dir / "atom_1_3.in", scaled(atom_points(symmetry, a, atom_z), bohr));
      } else {
        write_points(dir / "gasphase.in",
                     gas_phase_points());  // Zcc = ZsoMin exactly would stop after unit round trips
      }
    }
    std::cout << "wrote " << dir.string() << "\n";
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

// crplib_make_dataset <case> <case_dir>: writes the synthetic data set of a reference case
// (<case_dir>/dataset/) and the point lists the Fortran reference driver evaluates for it
// (<case_dir>/<mode>[_<length>_<energy>].in, hexadecimal doubles). Used by
// tools/fortran_reference/run_crp_reference.sh; the tests read the same files.
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "crp_synthetic.hpp"
#include "crplib/dataset.hpp"

namespace {

using crplib::Symmetry;

std::string hex(double x) { return std::format("{:016x}", std::bit_cast<std::uint64_t>(x)); }

void write_points(const std::filesystem::path& file, const std::vector<std::vector<double>>& points) {
  std::ofstream out(file);
  out << points.size() << '\n';
  for (const auto& p : points) {
    for (std::size_t i = 0; i < p.size(); ++i) out << (i ? " " : "") << hex(p[i]);
    out << '\n';
  }
}

// Cartesian (bohr) of the lattice coordinates (u, v): skewed for C3v/C6v, Cartesian for C4v.
std::pair<double, double> cartesian(Symmetry symmetry, double u, double v, double a) {
  if (symmetry == Symmetry::c4v) return {u * a, v * a};
  return {(u + v / 2.0) * a, v * std::sqrt(3.0) / 2.0 * a};
}

// Atom positions: lattice points (incl. negative ones, the ORIGEN quirk), points outside the cell,
// fold boundaries, and z at the ends and inside the data range (first/last interval: end slopes).
std::vector<std::vector<double>> atom_points(Symmetry symmetry, const crplib::DataSet& d) {
  const std::vector<double>& z = d.atom_cuts.front().potential.x;
  const double z_first = z.front(), z_last = z.back();
  const std::vector<double> heights = {z_first, z_first + 0.05, 0.37, 1.9, 4.4, z_last - 0.05, z_last};
  std::vector<double> coordinates = {-2.0, -1.0, -0.37, 0.0, 0.21, 0.25, 0.5, 0.75, 0.9, 1.0, 1.73};
  std::vector<std::vector<double>> points;
  for (double u : coordinates) {
    for (double v : coordinates) {
      for (double h : heights) {
        // a sparse but deterministic subset of the product grid
        const int key = static_cast<int>(std::lround(100 * u + 37 * v + 11 * h));
        if ((key % 3 + 3) % 3 != 0 && !(u == v)) continue;
        const auto [x, y] = cartesian(symmetry, u, v, d.lattice_constant);
        points.push_back({x, y, h});
      }
    }
  }
  if (symmetry == Symmetry::c4v) {  // REGION boundaries: x = y, x = a/2, y = a/2, x + y = a
    const double a = d.lattice_constant;
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

crplib::DataSet case_data_set(const std::string& name) {
  using namespace crplib_test;
  if (name == "c6v_pairwise") return pairwise_data_set(Symmetry::c6v);
  if (name == "c3v_pairwise") return pairwise_data_set(Symmetry::c3v);
  if (name == "c4v_pairwise") return pairwise_data_set(Symmetry::c4v);
  if (name == "c6v_pairwise_angstrom") return in_angstrom_and_ev(pairwise_data_set(Symmetry::c6v));
  if (name == "c4v_negative_spacing") {
    crplib::DataSet d = pairwise_data_set(Symmetry::c4v);
    d.interlayer_spacings = {-0.98};  // relative: z of the second layer = 0.98 a / sqrt(2)
    return d;
  }
  throw std::runtime_error("unknown case " + name);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: crplib_make_dataset <case> <case_dir>\n";
    return 1;
  }
  try {
    const std::string name = argv[1];
    const std::filesystem::path dir = argv[2];
    const crplib::DataSet d = case_data_set(name);
    std::filesystem::remove_all(dir / "dataset");
    crplib::write_data_set(d, dir / "dataset");
    const crplib::DataSet internal = crplib::read_data_set(dir / "dataset");  // the values as read
    // Points in bohr: convert the data set's lengths for the Angstrom case.
    crplib::DataSet in_bohr = internal;
    if (internal.units.length == crplib::LengthUnit::angstrom) {
      in_bohr.lattice_constant =
          crplib::convert_length(internal.lattice_constant, crplib::LengthUnit::angstrom, crplib::LengthUnit::bohr);
      for (auto& cut : in_bohr.atom_cuts)
        for (double& z : cut.potential.x)
          z = crplib::convert_length(z, crplib::LengthUnit::angstrom, crplib::LengthUnit::bohr);
    }
    write_points(dir / "atom_0_2.in", atom_points(d.symmetry, in_bohr));
    std::cout << "wrote " << dir.string() << "\n";
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  return 0;
}

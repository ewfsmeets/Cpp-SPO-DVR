// Validation helper (item 4 of docs/validation_thesis_fig5_3_dz): the C++ potential of a run,
// computed with the library of the build that ran, on the full grid.
// Usage: potential_check <input.json> <output dir>
// Writes into the output directory:
//   asymptotic_potential.dat, gas_phase_states.dat  as the run writes them (self-check: bitwise
//                                                    equal to the run's own files)
//   grids.txt              Z, Z specular, r, x, y, theta (+ weights), phi, k grids (%.17e)
//   v_stored.bin           V (hartree) as the run stores it: (itp, irep, izr), float64
//   symmetry_maps.txt      per ixy: representative, operation; per operation: the phi map
//   v_full.bin             V (hartree) evaluated at every full-grid point, no symmetry:
//                          (itp, ixy, izr), itp = iphi + nphi*itheta, ixy = iy + ny*ix
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <span>
#include <string>

#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/binary_file.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/kinetic_operators.hpp"
#include "report/sections.hpp"

using namespace spodvr;
namespace fs = std::filesystem;

namespace {

void write_values(std::ofstream& out, const std::string& name, std::span<const double> values) {
  out << name << ' ' << values.size() << '\n';
  char buffer[64];
  for (const double v : values) {
    std::snprintf(buffer, sizeof buffer, "%.17e\n", v);
    out << buffer;
  }
}

void write_grids(const fs::path& path, const SpatialGrids& g, const AngularGrid& a) {
  std::ofstream out(path);
  write_values(out, "z_scattering", g.z_scattering.points());
  write_values(out, "z_specular", g.z_specular.points());
  write_values(out, "r", g.r.points());
  write_values(out, "x", g.x.points());
  write_values(out, "y", g.y.points());
  write_values(out, "theta", a.thetas());
  write_values(out, "theta_weights", a.theta_weights());
  write_values(out, "phi", a.phis());
  write_values(out, "kz_scattering", g.kz_scattering.k());
  write_values(out, "kz_specular", g.kz_specular.k());
  write_values(out, "kr", g.kr.k());
  const double spacings[] = {g.kz_scattering.spacing(), g.kz_specular.spacing(), g.kr.spacing()};
  write_values(out, "dk_z_scattering_z_specular_r", spacings);
  out << "nz_potential " << g.nz_potential << "\niz_analysis " << g.iz_analysis << '\n';
}

void write_symmetry(const fs::path& path, const GridSymmetry& s) {
  std::ofstream out(path);
  out << "n_xy " << s.n_xy() << " n_representatives " << s.n_representatives() << " nphi " << s.nphi() << '\n';
  for (int ixy = 0; ixy < s.n_xy(); ++ixy)
    out << ixy << ' ' << s.representative_of(ixy) << ' ' << s.operation_of(ixy) << '\n';
  for (int irep = 0; irep < s.n_representatives(); ++irep)
    out << "rep " << irep << ' ' << s.representative_xy(irep) << '\n';
  for (int ixy = 0; ixy < s.n_xy(); ++ixy) {
    out << "phimap " << s.operation_of(ixy);
    for (const int i : s.phi_index(s.operation_of(ixy))) out << ' ' << i;
    out << '\n';
  }
}

void write_binary(const fs::path& path, std::span<const double> values) {
  BinaryWriter file(path);
  file.write_array(values);
  file.close();
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: potential_check <input.json> <output dir>\n";
    return 2;
  }
  const fs::path out_dir = argv[2];
  fs::create_directories(out_dir);
  const Config config = build_config(read_input_file(argv[1]));
  const auto pes = make_potential(config.input);
  const AngularGrid angular = make_angular_grid(config);
  const SpatialGrids grids = make_spatial_grids(config);
  write_grids(out_dir / "grids.txt", grids, angular);

  // As in app/run.cpp: the run's symmetry and work-unit groups
  GridSymmetry symmetry = make_grid_symmetry(config, angular);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  const ScatteringShape shape{basis.size(), config.input.nx, config.input.ny, config.input.nr, config.input.nz};
  const PotentialGrid stored =
      compute_potential_grid(*pes, config, grids, angular, std::move(symmetry), xy_group_points(shape));
  write_asymptotic_potential_file(out_dir / "asymptotic_potential.dat", grids, stored.asymptotic);
  const GasPhase gas = compute_gas_phase(config, grids, angular, stored.asymptotic, FftwFft(), LapackEigensolver());
  write_gas_phase_file(out_dir / "gas_phase_states.dat", grids, gas);
  write_symmetry(out_dir / "symmetry_maps.txt", stored.symmetry);
  write_binary(out_dir / "v_stored.bin", std::span<const double>(stored.values.data(), stored.values.size()));
  std::cout << "stored potential: " << stored.values.size() << " values, " << stored.symmetry.n_representatives()
            << " representatives\n";

  // The same PES at every full-grid point (identity symmetry)
  const PotentialGrid full =
      compute_potential_grid(*pes, config, grids, angular,
                             GridSymmetry(config.input.nx, config.input.ny, angular.n_phi()), xy_group_points(shape));
  write_binary(out_dir / "v_full.bin", std::span<const double>(full.values.data(), full.values.size()));
  std::cout << "full potential: " << full.values.size() << " values\n";
  return 0;
}

// Validation helper (item 18b of docs/validation_thesis_fig5_3_dz): feeds the BK coefficients of a
// run (native bk_coefficients.bin or Fortran fort.21) through the C++ on-the-fly analysis
// (SMatrixAccumulator) and writes the result tables, as the run writes them at its end.
// Usage: replay <input.json> <BK file> <native|fortran> <output dir>
// Self-check: replaying a C++ run's own file must reproduce its s_matrix.tsv bitwise.
// Precedent: the restart path (app/run.cpp restore_from_checkpoint, add_earlier_analysis_step).
#include <filesystem>
#include <iostream>
#include <string>

#include "analysis/s_matrix.hpp"
#include "app/fine_analysis_run.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/checkpoint.hpp"
#include "io/fortran_files.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/kinetic_operators.hpp"
#include "report/analysis_files.hpp"

using namespace spodvr;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
  if (argc != 5) {
    std::cerr << "usage: replay <input.json> <BK file> <native|fortran> <output dir>\n";
    return 2;
  }
  const fs::path bk_path = argv[2];
  const std::string format = argv[3];
  const fs::path out_dir = argv[4];
  fs::create_directories(out_dir);

  const Config config = build_config(read_input_file(argv[1]));
  const auto pes = make_potential(config.input);
  const AngularGrid angular = make_angular_grid(config);
  const SpatialGrids grids = make_spatial_grids(config);
  // V_asy as a propagation run obtains it (app/run.cpp: from the potential grid), not via
  // compute_asymptotic_potential (the onlyAnalysis path): the gas-phase energies, and through
  // them the phase factors exp(iEt) of the analysis, must be bitwise those of the run.
  const RotationalBasis rotational_basis(angular, config.input.j_in, config.input.mj_in);
  const ScatteringShape shape{rotational_basis.size(), config.input.nx, config.input.ny, config.input.nr,
                              config.input.nz};
  const std::vector<double> asymptotic =
      compute_potential_grid(*pes, config, grids, angular, make_grid_symmetry(config, angular), xy_group_points(shape))
          .asymptotic;
  const GasPhase gas = compute_gas_phase(config, grids, angular, asymptotic, FftwFft(), LapackEigensolver());
  const AnalysisSetup setup(config, angular, gas);

  SMatrixAccumulator accumulator(config, setup.basis, gas);
  int n_records = 0;
  double last_time = 0.0;
  auto add = [&](double time, const BkCoefficients& bk) {
    accumulator.add(bk, time);
    n_records = n_records + 1;
    last_time = time;
  };
  if (format == "native")
    (void)read_bk_file(bk_path, config, setup.basis, add);
  else if (format == "fortran")
    (void)read_fortran_bk_file(bk_path, config, setup.rotations, setup.basis, add);
  else {
    std::cerr << "unknown format " << format << '\n';
    return 2;
  }
  const ScatteringResults results = accumulator.results();
  const auto files = write_scattering_files(out_dir, config, setup.basis, results);
  std::cout << "replayed " << n_records << " records (last t = " << last_time << ") from " << bk_path.string()
            << "; wrote " << files.size() << " files to " << out_dir.string() << '\n';
  return 0;
}

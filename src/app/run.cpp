// Program driver: input -> grids -> potential -> gas-phase states -> propagation -> output.
// Fortran: Main.f90
#include "app/run.hpp"

#include <omp.h>

#include <chrono>
#include <filesystem>
#include <vector>

#include "core/error.hpp"
#include "core/numerical_libraries.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/propagation.hpp"
#include "report/propagation_report.hpp"
#include "report/sections.hpp"

namespace spodvr {

namespace {

// A new propagation (nread = 0) from tStart to tEnd; writes averages.dat.
// Fortran: Main.f90 from SetupRDBasis to the end of the propagation loop
void run_propagation(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                     const AngularGrid& angular, PotentialGrid potential, const GasPhase& gas, const Fft& fft,
                     const std::filesystem::path& output_dir) {
  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft);
  print_propagation_setup(out, style, grids, propagation);

  const AveragesTable table(out, style);
  table.print_header();
  std::vector<Averages> averages;
  const auto start = std::chrono::steady_clock::now();
  propagation.run(
      [&](const Averages& a) {
        table.print_row(a);
        averages.push_back(a);
      },
      [&](double z) { table.print_note("analysis surface on the specular grid fixed at Z = " + fixed(z, 6) + " bohr"); });
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;

  write_averages_file(output_dir / "averages.dat", averages);
  out << "\n  Propagation time: " << fixed(elapsed.count(), 1) << " s\n";
}

}  // namespace

// Fortran: Main.f90 (new runs and potential-only runs).
int run_program(const CommandLine& command_line, std::ostream& out, std::ostream& err, bool color) {
  const Style style{color};
  try {
    if (command_line.threads) omp_set_num_threads(*command_line.threads);

    out << style.emphasis("SPO-DVR (C++) version " SPODVR_VERSION) << "  -  " << omp_get_max_threads()
        << " OpenMP thread(s)\n";
    out << "  " << numerical_libraries_summary() << '\n';

    const Config config = build_config(read_input_file(command_line.input_file));
    print_input_summary(out, style, config);

    const auto potential = make_potential(config.input.system_id);
    const AngularGrid angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear,
                              config.input.j_in, config.input.mj_in);
    const SpatialGrids grids = make_spatial_grids(config);
    print_grids(out, style, config, grids, angular);

    if (config.derived.run_mode == 2)
      throw NotSupportedError("analysis-only runs (onlyAnalysis) are not implemented yet");

    PotentialGrid potential_grid = compute_potential_grid(*potential, config, grids, angular);
    print_asymptotic_potential(out, style, config, grids, potential_grid);

    const FftwFft fft;
    const LapackEigensolver eigensolver;
    const GasPhase gas = compute_gas_phase(config, grids, angular, potential_grid.asymptotic, fft, eigensolver);
    print_gas_phase(out, style, config, gas);

    std::filesystem::create_directories(command_line.output_dir);
    write_asymptotic_potential_file(command_line.output_dir / "asymptotic_potential.dat", grids, potential_grid);
    write_gas_phase_file(command_line.output_dir / "gas_phase_states.dat", grids, gas);
    out << "\n  Data files written to " << command_line.output_dir.string() << "\n";

    if (config.derived.run_mode == 6) {
      out << '\n' << style.emphasis("Potential-only run finished.") << '\n';
      return 0;
    }
    if (config.derived.run_mode == 1)
      throw NotSupportedError("restart runs (NewCalc = false) are not implemented yet");

    run_propagation(out, style, config, grids, angular, std::move(potential_grid), gas, fft, command_line.output_dir);
    out << '\n' << style.emphasis("Propagation finished (the asymptotic analysis is not implemented yet).") << '\n';
    return 0;
  } catch (const InputError& e) {
    err << style.warning("Input error") << ": " << e.what() << '\n';
    return 2;
  } catch (const NotSupportedError& e) {
    err << style.warning("Not supported") << ": " << e.what() << '\n';
    return 3;
  } catch (const std::exception& e) {
    err << style.warning("Error") << ": " << e.what() << '\n';
    return 1;
  }
}

}  // namespace spodvr

// Program driver: input -> grids -> potential -> gas-phase states -> output. Fortran: Main.f90
#include "app/run.hpp"

#include <omp.h>

#include <filesystem>

#include "core/error.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "report/sections.hpp"

namespace spodvr {

// Fortran: Main.f90, up to and including SetupAsympWF (the complete nread = 6 run).
int run_program(const CommandLine& command_line, std::ostream& out, std::ostream& err, bool color) {
  const Style style{color};
  try {
    if (command_line.threads) omp_set_num_threads(*command_line.threads);

    out << style.emphasis("SPO-DVR (C++) version " SPODVR_VERSION) << "  -  " << omp_get_max_threads()
        << " OpenMP thread(s)\n";

    const Config config = build_config(read_input_file(command_line.input_file));
    print_input_summary(out, style, config);

    const auto potential = make_potential(config.input.system_id);
    const AngularGrid angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear,
                              config.input.j_in, config.input.mj_in);
    const SpatialGrids grids = make_spatial_grids(config);
    print_grids(out, style, config, grids, angular);

    if (config.derived.run_mode == 2)
      throw NotSupportedError("analysis-only runs (onlyAnalysis) are not implemented yet");

    const PotentialGrid potential_grid = compute_potential_grid(*potential, config, grids, angular);
    print_asymptotic_potential(out, style, config, grids, potential_grid);

    const FftwFft fft;
    const LapackEigensolver eigensolver;
    const GasPhase gas = compute_gas_phase(config, grids, angular, potential_grid.asymptotic, fft, eigensolver);
    print_gas_phase(out, style, config, gas);

    std::filesystem::create_directories(command_line.output_dir);
    write_asymptotic_potential_file(command_line.output_dir / "asymptotic_potential.dat", grids, potential_grid);
    write_gas_phase_file(command_line.output_dir / "gas_phase_states.dat", grids, gas);
    out << "\n  Data files written to " << command_line.output_dir.string() << "\n";

    if (config.derived.run_mode != 6)
      throw NotSupportedError("propagation is not implemented yet; only potential-only runs (\"Potential\": true) "
                              "are supported");
    out << '\n' << style.emphasis("Potential-only run finished.") << '\n';
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

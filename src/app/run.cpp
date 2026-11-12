// Program driver: input -> grids -> potential -> gas-phase states -> propagation -> output.
// Fortran: Main.f90
#include "app/run.hpp"

#include <omp.h>

#include <chrono>
#include <filesystem>
#include <optional>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include "app/backend_factory.hpp"
#include "app/banner.hpp"
#include "app/fine_analysis_run.hpp"
#include "core/error.hpp"
#include "core/loop_schedule.hpp"
#include "core/numerical_libraries.hpp"
#include "core/real_compare.hpp"
#include "core/system_memory.hpp"
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
#include "potential/crp_potential.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"
#include "report/analysis_files.hpp"
#include "report/analysis_report.hpp"
#include "report/crp_report.hpp"
#include "report/fine_analysis_report.hpp"
#include "report/propagation_report.hpp"
#include "report/sections.hpp"

namespace spodvr {

namespace {

// Checkpoint file names in an output / restart directory: native, and as written by the
// Fortran code (fort.9, fort.21; a Fortran restart reads them as fort.8, fort.22).
constexpr std::string_view checkpoint_file = "checkpoint.bin";
constexpr std::string_view bk_file = "bk_coefficients.bin";
constexpr std::string_view fortran_wave_function_file = "fort.9";
constexpr std::string_view fortran_bk_file = "fort.21";

// BK coefficients of the earlier runs, in time order.
using BkRecords = std::vector<std::pair<double, BkCoefficients>>;

// The CRP data set and the grid points outside its range: the whole r grid and the Z points of the
// 6D potential (none in analysis-only runs).
void report_crp_data_set(std::ostream& out, const Style& style, const CrpPotential& crp, const SpatialGrids& grids,
                         bool analysis_only) {
  const std::span<const double> z_6d = grids.z_scattering.points().first(analysis_only ? 0 : grids.nz_potential);
  print_crp_data_set(out, style, crp.summary(), crp.range_report(z_6d, grids.r.points()));
}

// The first of the files that exists in dir (the Fortran input names fort.8 / fort.22 or the
// output names fort.9 / fort.21).
std::filesystem::path first_existing(const std::filesystem::path& dir, std::string_view a, std::string_view b) {
  if (std::filesystem::exists(dir / a)) return dir / a;
  if (std::filesystem::exists(dir / b)) return dir / b;
  throw CheckpointError("neither " + std::string(a) + " nor " + std::string(b) + " found in '" + dir.string() + "'");
}

// Restart (nread = 1): reads the wave packet and analysis surface of the earlier run and its BK
// coefficients, which are added to the S-matrix and returned (to be copied into the new BK files).
// Fortran: WaveFunction.f90 ReadWaveFunc, AsympInfo.f90 ReadOrWriteAsympInfo (nread = 1),
// Analyzor.f90 GetSMatrixElemsOfOldRun
BkRecords restore_from_checkpoint(std::ostream& out, const Config& config, WavePacketPropagation& propagation) {
  const InputParameters& p = config.input;
  const std::filesystem::path dir = p.restart_dir;
  BkRecords records;
  auto add_record = [&](double time, const BkCoefficients& bk) {
    propagation.add_earlier_analysis_step(time, bk);
    records.emplace_back(time, bk);
  };

  double checkpoint_time = 0.0;
  if (p.restart_format == "native") {
    CheckpointInfo info{};
    propagation.read_wave_packet([&](WavePacket& psi) { info = read_checkpoint(dir / checkpoint_file, config, psi); });
    propagation.restore_analysis_surface(info.surface);
    checkpoint_time = info.time;
    const auto last = read_bk_file(dir / bk_file, config, propagation.analysis_basis(), add_record);
    if (last && fortran_compare::rne(*last, checkpoint_time))
      throw CheckpointError("the BK coefficient file in '" + dir.string() + "' ends at t = " + full_precision(*last) +
                            ", the checkpoint at t = " + full_precision(checkpoint_time));
  } else {
    // The Fortran wave function file has no time: the BK records give it (n deltat).
    const AnalysisSurfaceState fresh = propagation.analysis_surface().state();
    AnalysisSurfaceState state{};
    propagation.read_wave_packet([&](WavePacket& psi) {
      state = read_fortran_wave_function(first_existing(dir, "fort.8", "fort.9"), psi, fresh);
    });
    propagation.restore_analysis_surface(state);
    read_fortran_bk_file(first_existing(dir, "fort.22", "fort.21"), config, propagation.rotational_basis(),
                         propagation.analysis_basis(), add_record);
    checkpoint_time = records.empty() ? 0.0 : records.back().first;
  }
  if (fortran_compare::rne(checkpoint_time, p.t_start))
    throw CheckpointError("tStart = " + full_precision(p.t_start) + " but the checkpoint in '" + dir.string() +
                          "' is at t = " + full_precision(checkpoint_time));
  out << "\n  Restarted from " << dir.string() << " (" << p.restart_format
      << " files) at t = " << fixed(checkpoint_time, 2) << " a.u.\n";
  return records;
}

// A restart must not overwrite the files it continues from.
void check_restart_directory(const Config& config, const std::filesystem::path& output_dir) {
  if (config.derived.run_mode != 1) return;
  const std::filesystem::path dir = config.input.restart_dir;
  if (std::filesystem::exists(dir) && std::filesystem::exists(output_dir) &&
      std::filesystem::equivalent(dir, output_dir))
    throw InputError("RestartDir",
                     "the output directory must differ from RestartDir (the checkpoint would be overwritten)");
}

// The BK coefficient files of this run: native always, Fortran fort.21 on request.
struct BkOutput {
  BkFileWriter native;
  std::optional<FortranBkFileWriter> fortran;

  void append(double time, const BkCoefficients& bk) {
    native.append(time, bk);
    if (fortran) fortran->append(bk);
  }
};

// A propagation from tStart to tEnd (new run or restart) followed by the asymptotic analysis;
// flux_slot_count: MemoryEstimate::flux_slot_count.
// Fortran: Main.f90 from SetupRDBasis to StateToStateProbs
void run_propagation(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                     const AngularGrid& angular, PotentialGrid potential, const GasPhase& gas, const Fft& fft,
                     const std::filesystem::path& output_dir, Device device, int flux_slot_count) {
  check_restart_directory(config, output_dir);
  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft, backend_factory(device),
                                    flux_slot_count);
  print_propagation_setup(out, style, grids, propagation);
  print_analysis_setup(out, style, config, propagation.analysis_basis(), propagation.analysis_energies());

  const BkRecords earlier =
      config.derived.run_mode == 1 ? restore_from_checkpoint(out, config, propagation) : BkRecords{};
  const bool write_native = config.input.checkpoint_format != "fortran";
  const bool write_fortran = config.input.checkpoint_format != "native";
  BkOutput bk_out{BkFileWriter(output_dir / bk_file, config, propagation.analysis_basis()), std::nullopt};
  if (write_fortran) {
    const int n_records = static_cast<int>(earlier.size()) + config.input.n_steps / config.derived.steps_per_analysis;
    bk_out.fortran.emplace(output_dir / fortran_bk_file, config, propagation.rotational_basis(),
                           propagation.analysis_basis(), gas, n_records);
  }
  for (const auto& [time, bk] : earlier) bk_out.append(time, bk);

  const InputParameters& p = config.input;
  AveragesTable table(out, style, p.n_steps, p.time_step, p.t_end - p.n_steps * p.time_step);
  table.print_header();
  std::vector<Averages> averages;
  PropagationObserver observer;
  observer.on_averages = [&](const Averages& a) {
    table.print_row(a);
    averages.push_back(a);
  };
  observer.on_surface_fixed = [&](double time, double z) {
    table.print_note("analysis surface on the specular grid fixed at Z = " + fixed(z, 6) + " bohr", time);
  };
  observer.on_analysis_step = [&](double time, const BkCoefficients& bk) {
    bk_out.append(time, bk);
    table.print_status(time);
  };
  propagation.run(observer);
  table.finish();

  // Checkpoint for a restart. Fortran: SaveWaveFunc
  const AnalysisSurfaceState surface = propagation.analysis_surface().state();
  if (write_native)
    write_checkpoint(output_dir / checkpoint_file, config, config.input.t_end, propagation.wave_packet_on_host(),
                     surface);
  if (write_fortran)
    write_fortran_wave_function(output_dir / fortran_wave_function_file, propagation.wave_packet_on_host(), surface);

  // Asymptotic analysis: S-matrix and probabilities at the analysis energies.
  const ScatteringResults results = propagation.scattering_results();
  print_scattering_results(out, style, config, propagation.analysis_basis(), results);

  write_averages_file(output_dir / "averages.dat", averages);
  const auto files = write_scattering_files(output_dir, config, propagation.analysis_basis(), results);
  out << "\n  Data files written to " << output_dir.string() << ": averages.dat";
  for (const auto& f : files) out << ", " << f.filename().string();
  out << ", " << bk_file;
  if (write_native) out << ", " << checkpoint_file;
  if (write_fortran) out << ", " << fortran_bk_file << ", " << fortran_wave_function_file;
  out << '\n';
}

// The fine analysis of a BK file (run mode 2, or alsoAnalysis after the propagation): grid,
// memory estimate, spectrum, interpolated results and the files in <output>/fine_analysis/.
// Fortran: Main.f90 (nread = 2: gtaife, probaf, StateToStateProbs)
void run_fine(std::ostream& out, const Style& style, const Config& config, const AngularGrid& angular,
              const GasPhase& gas, const Fft& fft, const std::filesystem::path& output_dir, const BkInput& input) {
  const AnalysisSetup setup(config, angular, gas);
  if (config.derived.run_mode == 2)
    print_analysis_setup(out, style, config, setup.basis, make_analysis_energies(config, gas));
  out << "\n  Fine analysis of " << input.path.string() << (input.fortran_format ? " (Fortran file)" : "") << '\n';
  const auto start = std::chrono::steady_clock::now();
  const FineAnalysisResults results = run_fine_analysis(config, setup, gas, input, fft, [&](const FineAnalysisGrid& g) {
    print_fine_analysis_grid(out, style, config, g);
    print_memory_estimate(out, style, estimate_fine_analysis_memory(g, setup.basis), physical_memory_bytes(),
                          "MEMORY ESTIMATE (FINE ANALYSIS)");
  });
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
  print_fine_analysis_results(out, style, config, setup.basis, results);
  out << "\n  Fine analysis time: " << fixed(elapsed.count(), 1) << " s\n";

  const auto files = write_fine_analysis_files(output_dir / "fine_analysis", config, setup.basis, results);
  out << "\n  Fine-analysis files written to " << (output_dir / "fine_analysis").string() << ":";
  for (const auto& f : files) out << ' ' << std::filesystem::relative(f, output_dir / "fine_analysis").string();
  out << '\n';
}

// The BK file of the run to analyse (RestartDir, RestartFormat). Fortran: fort.22
BkInput analysis_input(const Config& config) {
  const std::filesystem::path dir = config.input.restart_dir;
  if (config.input.restart_format == "native") return {dir / bk_file, false};
  return {first_existing(dir, "fort.22", "fort.21"), true};
}

}  // namespace

// Fortran: Main.f90 (new runs, restarts, analysis-only and potential-only runs).
int run_program(const CommandLine& command_line, std::ostream& out, std::ostream& err, bool color) {
  const Style style{color};
  try {
    if (command_line.threads) omp_set_num_threads(*command_line.threads);
    // The schedule of the (Z, r) group loops, set before any large array is first written.
    const int numa_nodes = numa_node_count();
    const LoopSchedule schedule =
        command_line.schedule ? parse_loop_schedule(*command_line.schedule) : default_loop_schedule(numa_nodes);
    apply_loop_schedule(schedule);

    out << banner() << '\n';
    out << style.emphasis("SPO-DVR (C++) version " SPODVR_VERSION) << "  -  " << omp_get_max_threads()
        << " OpenMP thread(s), loop schedule " << to_string(schedule) << " (" << numa_nodes << " NUMA node"
        << (numa_nodes > 1 ? "s" : "") << ")\n";
    out << "  " << numerical_libraries_summary() << '\n';

    const Config config = build_config(read_input_file(command_line.input_file));
    // Checked before the (possibly long) set-up: the GPU must be available in this build.
    const Device device = parse_device(command_line.device);
    const bool propagates = config.derived.run_mode != 2 && config.derived.run_mode != 6;
    if (device == Device::Cuda && !propagates)
      out << "  Note: --device cuda has no effect in this run mode (no propagation).\n";
    if (device == Device::Cuda && propagates) out << "  " << prepare_device(device) << '\n';
    print_input_summary(out, style, config);

    const auto potential = make_potential(config.input);
    const AngularGrid angular = make_angular_grid(config);
    const SpatialGrids grids = make_spatial_grids(config);
    // The flux analysis is not ported: a propagating run refuses it once its input checks have
    // passed. The developer tools and the tests' run_cpp_propagation bypass this on purpose.
    if (config.input.do_flux_analysis && propagates)
      throw NotSupportedError("doFluxAnalysis = true: the flux analysis is not implemented yet");
    print_grids(out, style, config, grids, angular);

    // Analysis-only runs (onlyAnalysis) need no 6D potential: only V(r) for the gas phase.
    const bool analysis_only = config.derived.run_mode == 2;
    if (fortran_compare::seq(potential->name(), "CRPlibrary"))
      report_crp_data_set(out, style, dynamic_cast<const CrpPotential&>(*potential), grids, analysis_only);
    std::optional<PotentialGrid> potential_grid;
    std::vector<double> asymptotic;
    int flux_slot_count = 0;  // the flux analysis's slots in the memory estimate; 0 without it
    if (analysis_only) {
      asymptotic = compute_asymptotic_potential(*potential, config, grids);
    } else {
      // Surface symmetry (UseSymmetry): check the PES before the reduced potential is computed.
      GridSymmetry symmetry = make_grid_symmetry(config, angular);
      if (config.input.use_symmetry) {
        const SymmetryCheck check = check_potential_symmetry(*potential, symmetry, config, grids, angular);
        print_surface_symmetry(out, style, symmetry, check);
        require_symmetric_potential(check, config, symmetry);
      }
      const MemoryEstimate estimate = estimate_memory(config, grids, angular, symmetry, omp_get_max_threads());
      print_memory_estimate(out, style, estimate, physical_memory_bytes());
      flux_slot_count = estimate.flux_slot_count;
      check_device_memory(out, style, propagates ? device : Device::Cpu, config, grids, angular,
                          symmetry.n_representatives());
      // Written per (Z, r) group of the propagator (first touch, core/loop_schedule.hpp).
      const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
      const ScatteringShape shape{basis.size(), config.input.nx, config.input.ny, config.input.nr, config.input.nz};
      const auto potential_start = std::chrono::steady_clock::now();
      potential_grid =
          compute_potential_grid(*potential, config, grids, angular, std::move(symmetry), xy_group_points(shape));
      const std::chrono::duration<double> potential_time = std::chrono::steady_clock::now() - potential_start;
      out << "\n  6D potential computed in " << fixed(potential_time.count(), 2) << " s\n";
      asymptotic = potential_grid->asymptotic;
    }
    print_asymptotic_potential(out, style, config, grids, asymptotic);

    const FftwFft fft;
    const LapackEigensolver eigensolver;
    const GasPhase gas = compute_gas_phase(config, grids, angular, asymptotic, fft, eigensolver);
    print_gas_phase(out, style, config, gas);

    std::filesystem::create_directories(command_line.output_dir);
    write_asymptotic_potential_file(command_line.output_dir / "asymptotic_potential.dat", grids, asymptotic);
    write_gas_phase_file(command_line.output_dir / "gas_phase_states.dat", grids, gas);
    out << "\n  Data files written to " << command_line.output_dir.string() << "\n";

    if (config.derived.run_mode == 6) {
      out << '\n' << style.emphasis("Potential-only run finished.") << '\n';
      return 0;
    }
    if (analysis_only) {
      run_fine(out, style, config, angular, gas, fft, command_line.output_dir, analysis_input(config));
      out << '\n' << style.emphasis("Analysis-only run finished.") << '\n';
      return 0;
    }

    run_propagation(out, style, config, grids, angular, std::move(*potential_grid), gas, fft, command_line.output_dir,
                    device, flux_slot_count);
    // alsoAnalysis: the fine analysis of the BK file this run has just written (the
    // propagation's memory is freed by now). Fortran: Main.f90 reads fort.21 back.
    if (config.input.also_analysis)
      run_fine(out, style, config, angular, gas, fft, command_line.output_dir,
               BkInput{command_line.output_dir / bk_file, false});
    out << '\n' << style.emphasis("Run finished.") << '\n';
    return 0;
  } catch (const InputError& e) {
    err << style.warning("Input error") << ": " << e.what() << '\n';
    return 2;
  } catch (const CheckpointError& e) {
    err << style.warning("Checkpoint error") << ": " << e.what() << '\n';
    return 4;
  } catch (const NotSupportedError& e) {
    err << style.warning("Not supported") << ": " << e.what() << '\n';
    return 3;
  } catch (const std::exception& e) {
    err << style.warning("Error") << ": " << e.what() << '\n';
    return 1;
  }
}

}  // namespace spodvr

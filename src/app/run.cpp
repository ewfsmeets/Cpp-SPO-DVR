// Program driver: input -> grids -> potential -> gas-phase states -> propagation -> output.
// Fortran: Main.f90
#include "app/run.hpp"

#include <omp.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <optional>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include "analysis/norm_budget.hpp"
#include "app/backend_factory.hpp"
#include "app/banner.hpp"
#include "app/fine_analysis_run.hpp"
#include "core/error.hpp"
#include "core/loop_schedule.hpp"
#include "core/numerical_libraries.hpp"
#include "core/real_compare.hpp"
#include "core/system_memory.hpp"
#include "fft/fftw_fft.hpp"
#include "frames/frame_source.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/checkpoint.hpp"
#include "io/fortran_files.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/pes_cut.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/kinetic_operators.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"
#include "report/analysis_files.hpp"
#include "report/analysis_report.hpp"
#include "report/densities_files.hpp"
#include "report/energy_resolved_files.hpp"
#include "report/fine_analysis_report.hpp"
#include "report/frame_files.hpp"
#include "report/norm_budget_files.hpp"
#include "report/pes_cut_file.hpp"
#include "report/plot_data.hpp"
#include "report/potential_files.hpp"
#include "report/potential_report.hpp"
#include "report/propagation_report.hpp"
#include "report/run_info_file.hpp"
#include "report/sections.hpp"

namespace spodvr {

namespace {

// Checkpoint file names in an output / restart directory: native, and as written by the
// Fortran code (fort.9, fort.21; a Fortran restart reads them as fort.8, fort.22).
constexpr std::string_view checkpoint_file = "checkpoint.bin";
constexpr std::string_view bk_file = "bk_coefficients.bin";
constexpr std::string_view energy_resolved_sums_file = "energy_resolved/sums.bin";  // EnergyResolvedSites
constexpr std::string_view fortran_wave_function_file = "fort.9";
constexpr std::string_view fortran_bk_file = "fort.21";
constexpr std::string_view fortran_wave_function_input = "fort.8";
constexpr std::string_view fortran_bk_input = "fort.22";

// BK coefficients of the earlier runs, in time order.
using BkRecords = std::vector<std::pair<double, BkCoefficients>>;

// The first of the files that exists in dir (the Fortran input names fort.8 / fort.22 or the
// output names fort.9 / fort.21).
std::filesystem::path first_existing(const std::filesystem::path& dir, std::string_view a, std::string_view b) {
  if (std::filesystem::exists(dir / a)) return dir / a;
  if (std::filesystem::exists(dir / b)) return dir / b;
  throw CheckpointError("neither " + std::string(a) + " nor " + std::string(b) + " found in '" + dir.string() + "'");
}

// A file a run reads from RestartDir: its names in the order first_existing tries them, and what
// it holds (for the error message).
struct RestartFile {
  std::vector<std::string_view> names;
  std::string_view holds;
};

// The files run mode 1 (restart) or 2 (analysis-only) reads from RestartDir, as RestartFormat
// selects, in the order restore_from_checkpoint and analysis_input open them.
std::vector<RestartFile> restart_files(const Config& config) {
  const bool native = config.input.restart_format == "native";
  const RestartFile bk = native ? RestartFile{{bk_file}, "the BK coefficients"}
                                : RestartFile{{fortran_bk_input, fortran_bk_file}, "the BK coefficients"};
  const RestartFile wave_packet =
      native ? RestartFile{{checkpoint_file}, "the wave packet"}
             : RestartFile{{fortran_wave_function_input, fortran_wave_function_file}, "the wave packet"};
  // The analysis-only run reads only the analysis surface of the checkpoint: the specular analysis line.
  if (config.derived.run_mode == 2) return {bk, RestartFile{wave_packet.names, "the specular analysis line"}};
  if (!config.input.energy_resolved_sites.empty())
    return {wave_packet, bk, RestartFile{{energy_resolved_sums_file}, "the energy-resolved sums"}};
  return {wave_packet, bk};
}

// The status of a path in RestartDir. An error other than "not found" (e.g. RestartDir without
// search permission) is an InputError naming RestartDir, not a filesystem_error (exit code 1).
std::filesystem::file_status restart_path_status(const std::filesystem::path& path) {
  std::error_code error;
  const std::filesystem::file_status status = std::filesystem::status(path, error);
  if (error && status.type() != std::filesystem::file_type::not_found)
    throw InputError("RestartDir", "cannot access '" + path.string() + "': " + error.message());
  return status;
}

// The first of the file's names that exists in dir (the one its reader opens) must be a readable
// regular file; if none exists, the error names them all.
void require_restart_file(const std::filesystem::path& dir, const RestartFile& file, const std::string& format) {
  for (const std::string_view name : file.names) {
    const std::filesystem::path path = dir / name;
    const std::filesystem::file_status status = restart_path_status(path);
    if (!std::filesystem::exists(status)) continue;
    if (!std::filesystem::is_regular_file(status) || !std::ifstream(path, std::ios::binary))
      throw InputError("RestartDir",
                       "'" + path.string() + "' (" + std::string(file.holds) + ") is not a readable regular file");
    return;
  }
  std::string missing = "no " + std::string(file.names[0]);
  if (file.names.size() == 2) missing = "neither " + std::string(file.names[0]) + " nor " + std::string(file.names[1]);
  throw InputError("RestartDir", missing + " in '" + dir.string() + "' (" + std::string(file.holds) +
                                     ", RestartFormat = \"" + format + "\")");
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
      state = read_fortran_wave_function(first_existing(dir, fortran_wave_function_input, fortran_wave_function_file),
                                         psi, fresh);
    });
    propagation.restore_analysis_surface(state);
    read_fortran_bk_file(first_existing(dir, fortran_bk_input, fortran_bk_file), config, propagation.rotational_basis(),
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

// Run modes 1 and 2, right after the input and before the (possibly long) set-up: RestartDir
// holds every file the run will read (the readers check their contents later), and a restart
// does not overwrite the files it continues from. A write2DPES run reads nothing from RestartDir
// (the Fortran stops after the cut, before it reads any restart file), but it writes run_info.json
// and asymptotic_potential.dat into the output directory: never over those of the run it names.
void check_restart_directory(const Config& config, const std::filesystem::path& output_dir) {
  if (config.derived.run_mode != 1 && config.derived.run_mode != 2) return;
  const std::filesystem::path dir = config.input.restart_dir;
  if (config.input.write_2d_pes) {
    std::error_code not_comparable;  // e.g. no RestartDir: then the two are not the same directory
    if (std::filesystem::equivalent(dir, output_dir, not_comparable))
      throw InputError("RestartDir",
                       "a write2DPES run writes into the output directory, which must differ from "
                       "RestartDir (its run_info.json would be replaced)");
    return;
  }
  if (!std::filesystem::is_directory(restart_path_status(dir)))
    throw InputError("RestartDir", "'" + dir.string() + "' is not a directory");
  for (const RestartFile& file : restart_files(config)) require_restart_file(dir, file, config.input.restart_format);
  if (config.derived.run_mode == 1 && !config.input.energy_resolved_sites.empty())
    check_energy_resolved_restart(dir, config.input);  // before the set-up, which can take hours
  // An error here (e.g. no output directory yet) means the two are not the same directory.
  std::error_code not_comparable;
  if (config.derived.run_mode == 1 && std::filesystem::equivalent(dir, output_dir, not_comparable))
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

// A propagation from tStart to tEnd (new run or restart) followed by the asymptotic analysis.
// Fortran: Main.f90 from SetupRDBasis to StateToStateProbs
void run_propagation(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                     const AngularGrid& angular, PotentialGrid potential, const GasPhase& gas, const Fft& fft,
                     const std::filesystem::path& output_dir, Device device,
                     const std::optional<BackendFactory>& backend) {
  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft,
                                    backend ? *backend : backend_factory(device));
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
  std::optional<DensitiesFiles> densities;  // writeDensities: the averages' own densities, every dtAve
  if (p.write_densities) {
    densities.emplace(output_dir, config, grids);
    observer.on_average_densities = [&](const Averages& a, std::span<const double> block_densities,
                                        const SpecularWave& specular, int surface_index) {
      if (!densities) return;
      try {
        densities->append(a, cell_probabilities(block_densities, specular, config), surface_index);
      } catch (const CheckpointError& e) {  // a full disk: optional output, so the propagation goes on
        table.print_note(std::string("writeDensities: ") + e.what() + "; densities/ stops here", a.time);
        densities.reset();
      }
    };
  }
  std::optional<FrameSource> frame_source;  // dtFrames: a frame of the wave packet every dtFrames
  std::optional<FrameFiles> frames;
  int frames_written = 0;
  if (p.dt_frames) {
    frame_source.emplace(config, grids, angular, propagation.rotational_basis(), gas, propagation.tables(), fft);
    frames.emplace(output_dir, config, grids, angular, propagation.rotational_basis(), frame_source->layout());
    observer.on_frame = [&](double time, const PropagationBackend& wave_packet, int surface_index) {
      if (!frames) return;
      const FrameRecord& record = frame_source->frame(wave_packet);
      try {
        frames->append(time, record, surface_index);
        ++frames_written;
      } catch (const CheckpointError& e) {  // a full disk: optional output, so the propagation goes on
        table.print_note(std::string("dtFrames: ") + e.what() + "; frames/ stops here", time);
        frames.reset();
      }
    };
  }
  std::optional<NormBudget> budget;  // NormBudget: what each absorber removes per analysis window
  std::optional<NormBudgetFiles> budget_files;
  int windows_written = 0;
  if (p.norm_budget) {
    budget.emplace(config, grids, propagation.tables());
    budget_files.emplace(output_dir, config, budget->absorbers());
    observer.on_analysis_densities = [&](double time, std::span<const double> block_densities,
                                         const SpecularWave& specular) {
      if (!budget_files) return;
      const std::optional<NormBudgetWindow> window = budget->add(cell_probabilities(block_densities, specular, config));
      if (!window) return;  // tStart: the first window starts here
      try {
        budget_files->append(time, *window);
        ++windows_written;
      } catch (const CheckpointError& e) {  // a full disk: optional output, so the propagation goes on
        table.print_note(std::string("NormBudget: ") + e.what() + "; norm_budget/ stops here", time);
        budget_files.reset();
      }
    };
  }
  std::optional<EnergyResolvedFiles> sums_files;  // EnergyResolvedSites/Energies: energy_resolved/
  if (!p.energy_resolved_sites.empty()) {
    const double t_start = p.t_end - p.n_steps * p.time_step;
    EnergyResolvedSums start = empty_energy_resolved_sums(config, gas, propagation.rotational_basis().size(),
                                                          static_cast<int>(grids.z_specular.size()), t_start);
    if (config.derived.run_mode == 1)
      start = read_energy_resolved_sums(p.restart_dir, config, start, static_cast<int>(earlier.size()));
    sums_files.emplace(output_dir, config, grids, start, propagation.tables().specular_to_scattering);
    propagation.start_energy_resolved(std::move(start));
    observer.on_energy_resolved_partial = [&](const EnergyResolvedSums& sums) {
      if (!sums_files) return;
      try {
        sums_files->write_partial(sums);
      } catch (const CheckpointError& e) {  // a full disk: optional output, so the propagation goes on
        table.print_note(std::string("EnergyResolvedSites: ") + e.what() + "; energy_resolved/ stops here",
                         sums.last_time);
        sums_files.reset();
      }
    };
  }
  const auto loop_start = std::chrono::steady_clock::now();
  propagation.run(observer);
  const std::chrono::duration<double> loop_time = std::chrono::steady_clock::now() - loop_start;
  table.finish();
  if (frame_source) {
    print_frame_cost(out, frames_written, frame_source->seconds(), loop_time.count());
    print_frame_extras_cost(out, frames_written, frame_source->extra_seconds());
  }
  if (budget) print_norm_budget_cost(out, windows_written, propagation.analysis_densities_seconds(), loop_time.count());
  if (propagation.sums_energy_resolved())
    print_energy_resolved_cost(out, p.energy_resolved_sites.size() * p.energy_resolved_energies.size(),
                               propagation.energy_resolved_seconds(), loop_time.count());

  // Checkpoint for a restart. Fortran: SaveWaveFunc
  const AnalysisSurfaceState surface = propagation.analysis_surface().state();
  if (write_native)
    write_checkpoint(output_dir / checkpoint_file, config, config.input.t_end, propagation.wave_packet_on_host(),
                     surface);
  if (write_fortran)
    write_fortran_wave_function(output_dir / fortran_wave_function_file, propagation.wave_packet_on_host(), surface);

  // The energy-resolved sums after the checkpoint, so that a full disk cannot cost the propagation's restart.
  if (sums_files) {
    try {
      sums_files->write_final(propagation.energy_resolved_sums());
    } catch (const CheckpointError& e) {
      out << "  EnergyResolvedSites: " << e.what() << "; energy_resolved/ is incomplete\n";
      sums_files.reset();
    }
  }

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
  if (densities) out << ", densities/";
  if (frames) out << ", frames/";
  if (budget_files) out << ", norm_budget/";
  if (sums_files) out << ", energy_resolved/";
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

// write2DPES: the cut of the PES, the gas-phase potential V_asy (the plots' energy reference; the Fortran
// writes none) and run_info.json (without the gas phase, which is not computed). Fortran: Main.f90
void write_pes_cut_run(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                       const AngularGrid& angular, const Potential& potential, const ProgramInfo& program,
                       const CommandLine& command_line) {
  const PesCut cut = compute_pes_cut(potential, config.input, grids);
  const std::vector<double> asymptotic = compute_asymptotic_potential(potential, config, grids);
  const std::filesystem::path& dir = command_line.output_dir;
  std::filesystem::create_directories(dir);
  write_pes_cut_file(dir / pes_cut_file_name, cut);
  write_asymptotic_potential_file(dir / "asymptotic_potential.dat", grids, asymptotic);
  write_run_info(dir / "run_info.json", config, grids, angular, program, command_line.input_file, std::nullopt);
  print_pes_cut(out, style, cut, dir / pes_cut_file_name);
}

// The BK file and the checkpoint of the run to analyse (RestartDir, RestartFormat). Fortran: fort.22
BkInput analysis_input(const Config& config) {
  const std::filesystem::path dir = config.input.restart_dir;
  if (config.input.restart_format == "native")
    return {dir / bk_file, false, CheckpointInput{dir / checkpoint_file, false}};
  return {first_existing(dir, fortran_bk_input, fortran_bk_file), true,
          CheckpointInput{first_existing(dir, fortran_wave_function_input, fortran_wave_function_file), true}};
}

// The BK file and the checkpoint that a propagation has just written into dir (alsoAnalysis).
BkInput written_analysis_input(const Config& config, const std::filesystem::path& dir) {
  const bool native = config.input.checkpoint_format != "fortran";  // as run_propagation's write_native
  const CheckpointInput checkpoint =
      native ? CheckpointInput{dir / checkpoint_file, false} : CheckpointInput{dir / fortran_wave_function_file, true};
  return {dir / bk_file, false, checkpoint};
}

}  // namespace

// Fortran: Main.f90 (new runs, restarts, analysis-only and potential-only runs).
int run_program(const CommandLine& command_line, std::ostream& out, std::ostream& err, bool color,
                const ProgramInfo& program, const std::optional<BackendFactory>& backend) {
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

    // The PES is built once, before the input checks that need its system parameters.
    const InputParameters input = read_input_file(command_line.input_file);
    const std::unique_ptr<Potential> potential = make_potential(input);
    const Config config = build_config(input, *potential);
    check_restart_directory(config, command_line.output_dir);
    // Checked before the (possibly long) set-up: the GPU must be available in this build.
    const Device device = parse_device(command_line.device);
    const bool propagates = config.derived.run_mode != 2 && config.derived.run_mode != 6 && !config.input.write_2d_pes;
    if (device == Device::Cuda && !propagates)
      out << "  Note: --device cuda has no effect in this run mode (no propagation).\n";
    if (device == Device::Cuda && propagates) out << "  " << prepare_device(device) << '\n';
    print_input_summary(out, style, config);

    const AngularGrid angular = make_angular_grid(config);
    const SpatialGrids grids = make_spatial_grids(config);
    // The flux analysis is not ported: a propagating run refuses it once its input checks have
    // passed. The developer tools and the tests' run_cpp_propagation bypass this on purpose.
    if (config.input.do_flux_analysis && propagates)
      throw NotSupportedError("doFluxAnalysis = true: the flux analysis is not implemented yet");
    print_grids(out, style, config, grids, angular);
    // write2DPES: the Fortran writes the cut right after the PES is set up and stops, in every run mode,
    // before the 6D potential and the gas phase (Main.f90).
    if (config.input.write_2d_pes) {
      write_pes_cut_run(out, style, config, grids, angular, *potential, program, command_line);
      out << '\n' << style.emphasis("2D cut of the PES written; the run stops here (write2DPES).") << '\n';
      return 0;
    }
    // Plot data (writeDensities, ...): its projected size, refused before the long set-up when the output file
    // system lacks room. Nothing is queried without a plot key.
    const std::vector<PlotDataSize> plot_data = plot_data_sizes(config, grids, angular);
    if (!plot_data.empty()) {
      const std::filesystem::path& dir = command_line.output_dir;
      check_replaceable(dir, plot_data);
      check_free_space(plot_data, free_space(dir) + replaced_plot_data_bytes(dir, plot_data));
      print_plot_data(out, style, plot_data);
    }

    // Analysis-only runs (onlyAnalysis) need no 6D potential: only V(r) for the gas phase.
    const bool analysis_only = config.derived.run_mode == 2;
    // The PES's own report (CRP data set, SCM diagnostics); SCM D3 can stop the run here.
    report_potential(out, style, *potential, config, grids, angular);
    std::optional<PotentialGrid> potential_grid;
    std::vector<double> asymptotic;
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
      print_memory_estimate(out, style, estimate_memory(config, grids, angular, symmetry, omp_get_max_threads()),
                            physical_memory_bytes());
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
    // An analysis-only run may write into the directory of the propagation it analyses: its metadata go
    // next to its own results, so that the propagation's run_info.json stays.
    const std::filesystem::path run_info_dir =
        analysis_only ? command_line.output_dir / "fine_analysis" : command_line.output_dir;
    std::filesystem::create_directories(run_info_dir);
    write_run_info(run_info_dir / "run_info.json", config, grids, angular, program, command_line.input_file,
                   std::cref(gas));
    out << "\n  Data files written to " << command_line.output_dir.string() << "\n";
    if (config.input.write_potential && potential_grid) {  // writePotential: V on the grid, as the run uses it
      write_potential_files(command_line.output_dir, config, grids, angular, *potential_grid);
      out << "  Plot data written to " << (command_line.output_dir / "potential").string() << "\n";
    }

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
                    device, backend);
    // alsoAnalysis: the fine analysis of the BK file this run has just written (the
    // propagation's memory is freed by now). Fortran: Main.f90 reads fort.21 back.
    if (config.input.also_analysis)
      run_fine(out, style, config, angular, gas, fft, command_line.output_dir,
               written_analysis_input(config, command_line.output_dir));
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

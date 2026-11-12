// One-off validation of the flux analysis (docs/validation_flux_replay/README.md): recomputes
// c_n psi_c and c_n d psi_c / dr of one Fortran analysis step from the Fortran's own wave function
// and compares them with its flux coefficient files. Fortran: Main.f90, FluxAnalyzor.f90 storeFluxCoeffs
//
// Usage: replay_flux <hardwall_small input.json> <Fortran part-2 run directory> [C++ output directory]
//   The run directory holds the files of the restarted one-step Fortran run (run_fortran.sh part2):
//   fort.9 (psi at the end of the step), fort.21 (eviroi), flux_coeffs.dat, flux_coeffs_deriv.dat
//   (phi and phi' after the step, summed onto zero-filled input files: c_n psi_c and c_n d psi_c / dr).
//   With an output directory, the C++ coefficients are also written there in the Fortran layout
//   and their bytes compared with the Fortran's files.
// Exit status: 0 every value equal, 1 a difference, 2 bad arguments. Run time: about a second.
#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <span>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/fortran_files.hpp"
#include "io/fortran_unformatted.hpp"
#include "propagation/analysis_surface.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/wave_function.hpp"

using namespace spodvr;
namespace fs = std::filesystem;
using cplx = std::complex<double>;

namespace {

// The part-2 input of run_fortran.sh: hardwall_small restarted (from run_dir) at t = 7000 for one
// analysis step (5 steps of 2), with the flux keys of part 1. Keep the two in step.
Config part2_config(const fs::path& case_input, const fs::path& run_dir) {
  std::ifstream file(case_input);
  nlohmann::json json = nlohmann::json::parse(file);
  json["NewCalc"] = false;
  json["RestartDir"] = run_dir.string();
  json["tStart"] = 7000.0;
  json["tEnd"] = 7010.0;
  json["nStep"] = 5;
  json["doFluxAnalysis"] = true;
  json["rFluxAnalysis"] = 2.0;
  json["EMinFluxAna"] = 0.2;
  json["EMaxFluxAna"] = 0.6;
  json["DeltaEFluxAna"] = 0.15;
  return build_config(parse_input(json));
}

// eviroi, the energy of the initial rovibrational state (hartree): record 2 of the Fortran BK file
// holds deltat, rkz, z0, alfwv, ZAnalysis, eviroi. Fortran: AsympInfo.f90 ReadOrWriteAsympInfo
double fortran_initial_state_energy(const fs::path& bk_file) {
  FortranRecordReader file(bk_file);
  (void)file.read("record 1 (sizes)");
  FortranRecordParser record = file.read("record 2 (deltat, rkz, z0, alfwv, ZAnalysis, eviroi)");
  for (int k = 0; k < 5; ++k) (void)record.read<double>();
  const double eviroi = record.read<double>();
  record.expect_end();
  return eviroi;
}

// phi and phi' of one analysis step ending at `time`, summed onto zero coefficients, from the wave
// function psi: the C++ accumulator with the reference r-derivative (flux_r_derivative).
FluxCoefficients replay_step(const ScatteringWave& psi, const FluxEnergies& energies, std::span<const double> weights,
                             int ir_flux, int n_planes, double time, double time_step) {
  const int n_blocks = psi.block_size();
  const int nr = static_cast<int>(weights.size());
  FluxAccumulator accumulator(n_blocks, n_planes, energies, time_step, 1);
  accumulator.add_step(time, [&](FluxSurfaceWave& slot) {
    std::vector<cplx> column(nr);
    for (int iz = 0; iz < n_planes; ++iz)
      for (int block = 0; block < n_blocks; ++block) {
        for (int ir = 0; ir < nr; ++ir) column[ir] = psi.block(ir + nr * iz)[block];
        slot.psi[slot.index(block, iz)] = column[ir_flux];
        slot.derivative[slot.index(block, iz)] = flux_r_derivative(weights, column);
      }
  });
  return accumulator.coefficients();
}

struct Agreement {
  long values = 0;   // complex values compared
  long nonzero = 0;  // of them, Fortran values other than (0, 0)
  long equal = 0;    // equal by value: real == real and imag == imag
  long bitwise = 0;  // equal bit for bit (the sign of zero included)
  double worst = 0.0;  // largest |C++ - Fortran| / max |Fortran| of one array
};

Agreement compare(std::span<const cplx> cpp, std::span<const cplx> fortran) {
  Agreement a;
  double scale = 0.0;
  for (const cplx v : fortran) scale = std::max(scale, std::abs(v));
  for (std::size_t k = 0; k < cpp.size(); ++k) {
    const cplx c = cpp[k];
    const cplx f = fortran[k];
    a.values = a.values + 1;
    if (f != cplx(0.0, 0.0)) a.nonzero = a.nonzero + 1;
    if (c.real() == f.real() && c.imag() == f.imag()) a.equal = a.equal + 1;
    if (std::bit_cast<std::uint64_t>(c.real()) == std::bit_cast<std::uint64_t>(f.real()) &&
        std::bit_cast<std::uint64_t>(c.imag()) == std::bit_cast<std::uint64_t>(f.imag()))
      a.bitwise = a.bitwise + 1;
    // NaN-aware: a NaN difference makes worst NaN (never smaller than a finite one).
    const double d = std::abs(c - f) / scale;
    if (!(d <= a.worst)) a.worst = d;
  }
  return a;
}

void print(const char* what, const Agreement& a) {
  std::printf(
      "%s %ld values (%ld nonzero): %ld equal by value, %ld bitwise; worst |C++ - Fortran| = %.3g of max "
      "|Fortran|\n",
      what, a.values, a.nonzero, a.equal, a.bitwise, a.worst);
}

std::string file_bytes(const fs::path& path) {
  std::ifstream file(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

int replay(const fs::path& case_input, const fs::path& run_dir, const fs::path& out_dir) {
  const Config config = part2_config(case_input, run_dir);
  const InputParameters& p = config.input;
  if (p.n_steps != config.derived.steps_per_analysis) {
    std::fprintf(stderr, "error: the replay covers one analysis step, the input has %d steps; fix part2_config\n",
                 p.n_steps);
    return 1;
  }
  const SpatialGrids grids = make_spatial_grids(config);
  const RotationalBasis rotations(make_angular_grid(config), p.j_in, p.mj_in);
  const int n_planes = grids.nz_potential;
  const int ir_flux = *grids.ir_flux_analysis;

  // psi at the end of the step: the Fortran saves it after the time loop (Main.f90 SaveWaveFunc),
  // and storeFluxCoeffs used it after Propagate and BalintKurtiCoeffs of the same step.
  WavePacket psi{ScatteringWave(rotations.size(), p.nx * p.ny, p.nz * p.nr), SpecularWave(p.nr, p.nz_specular)};
  (void)read_fortran_wave_function(run_dir / "fort.9", psi, AnalysisSurfaceState{0, 0, false, 0.0});

  // The time of the step as the C++ run computes it (propagation.cpp), = the Fortran's the_time =
  // tottym + dfloat(propStep - nstep) * dtlan = tEnd for the last step.
  const double time = p.t_end + (config.derived.steps_per_analysis - p.n_steps) * p.time_step;
  const double parallel_energy = parallel_kinetic_energies(config)[0];  // Fortran: edifin
  const double eviroi = fortran_initial_state_energy(run_dir / "fort.21");
  const FluxEnergies energies = make_flux_energies(config, parallel_energy, eviroi);
  const std::vector<double> weights = flux_derivative_weights(p.nr, grids.r.spacing(), ir_flux);
  const double time_step = config.derived.analysis_time_step;

  std::printf("time %.17g, analysis step %.17g, edifin %.17g, eviroi (Fortran) %.17g\n", time, time_step,
              parallel_energy, eviroi);
  std::printf(
      "%d flux energies, %d planes Z < zasymp, %d blocks (rotational functions x nX x nY), r_c = %.17g (ir %d)\n",
      energies.size(), n_planes, psi.scattering.block_size(), grids.r[ir_flux], ir_flux + 1);

  const FluxCoefficients cpp = replay_step(psi.scattering, energies, weights, ir_flux, n_planes, time, time_step);
  FluxCoefficients fortran(psi.scattering.block_size(), n_planes, energies.size());
  read_fortran_flux_coefficients(run_dir / "flux_coeffs.dat", run_dir / "flux_coeffs_deriv.dat", fortran);

  const Agreement phi = compare(cpp.phi, fortran.phi);
  const Agreement derivative = compare(cpp.phi_derivative, fortran.phi_derivative);
  print("phi = c_n psi_c:", phi);
  print("phi' = c_n d psi_c / dr:", derivative);

  // Control: the time of the step's start instead of its end must not match.
  const FluxCoefficients start =
      replay_step(psi.scattering, energies, weights, ir_flux, n_planes, time - time_step, time_step);
  print("control, t = start of the step, phi:", compare(start.phi, fortran.phi));

  bool bytes_identical = true;
  if (!out_dir.empty()) {
    fs::create_directories(out_dir);
    write_fortran_flux_coefficients(out_dir / "flux_coeffs.dat", out_dir / "flux_coeffs_deriv.dat", cpp);
    for (const char* name : {"flux_coeffs.dat", "flux_coeffs_deriv.dat"}) {
      const bool same = file_bytes(out_dir / name) == file_bytes(run_dir / name);
      std::printf("C++ %s byte-identical to the Fortran's: %s\n", name, same ? "yes" : "no");
      bytes_identical = bytes_identical && same;
    }
  }
  const bool all_equal = phi.equal == phi.values && derivative.equal == derivative.values && phi.values > 0;
  std::printf("result: %s\n", all_equal ? "every value equal" : "DIFFERENCES");
  return (all_equal && bytes_identical) ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3 || argc > 4) {
    std::fprintf(stderr,
                 "usage: replay_flux <hardwall_small input.json> <Fortran part-2 run directory> [C++ output "
                 "directory]\n");
    return 2;
  }
  try {
    return replay(argv[1], argv[2], argc == 4 ? fs::path(argv[3]) : fs::path());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s; check the arguments and the run directory\n", e.what());
    return 1;
  }
}

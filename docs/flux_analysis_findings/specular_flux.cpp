// Scratch measurement (not for merging): the flux through r = r_c of the specular wave (the initial
// channel), which the flux analysis leaves out. Accumulates phi_s and phi_s' of the specular wave with
// the run's c_n(t) and D(ir), embeds initXY phi_s into the (ichn, ixy) blocks of a copy of the run's flux
// coefficients and evaluates P_scat, P_spec and P_total with flux_results.
// Kept as a record: it builds only against the reverted flux port, see docs/flux_analysis_findings.md
// (comments edited, code unchanged).
// Usage: spodvr_q7_specular_flux <input.json> <output.tsv>   (threads: OMP_NUM_THREADS)
#include <chrono>
#include <cmath>
#include <complex>
#include <cstring>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <span>
#include <string>
#include <vector>

#include "analysis/flux_analysis.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"

using namespace spodvr;
using cplx = std::complex<double>;

namespace {

bool bitwise_equal(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

int run(const std::filesystem::path& input_path, const std::filesystem::path& output_path) {
  const auto start = std::chrono::steady_clock::now();
  const Config config = build_config(read_input_file(input_path));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  PotentialGrid potential = compute_potential_grid(*make_potential(config.input), config, grids, angular);
  const FftwFft fft;
  const GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver());
  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft, make_cpu_backend,
                                    estimate_memory(config, grids, angular, 1).flux_slot_count);
  if (!propagation.flux()) throw std::runtime_error("doFluxAnalysis is off");

  const int nr = config.input.nr;
  const int ir_c = *grids.ir_flux_analysis;
  const int n_planes = grids.nz_potential;
  const int n_rot = propagation.rotational_basis().size();
  const int n_xy = config.input.nx * config.input.ny;
  const int ichn = propagation.rotational_basis().initial_index();
  const FluxEnergies energies = propagation.flux()->energies();  // copy: the run's E_n, A_n
  const int n_energies = energies.size();
  const double time_step = config.derived.analysis_time_step;  // the FluxAccumulator's dt
  const std::vector<double>& weights = propagation.tables().flux_derivative;
  const std::vector<cplx>& init_xy = propagation.tables().specular_to_scattering;    // initXY
  const std::vector<cplx>& xy_proj = propagation.tables().scattering_to_specular;   // XYProj

  std::cout << std::format(
      "r_c = {:.17g} bohr (ir = {}), nz_potential = {}, z(nzp-1) = {:.6g}, ichn = {} (J {} mJ {}),"
      " n_rot = {}, nxy = {}, flux energies = {}, analysis dt = {}\n",
      grids.r[ir_c], ir_c, n_planes, grids.z_scattering[n_planes - 1], ichn, propagation.rotational_basis().j(ichn),
      propagation.rotational_basis().mj(ichn), n_rot, n_xy, n_energies, time_step);
  cplx sum_proj_init(0.0, 0.0);
  for (int ixy = 0; ixy < n_xy; ++ixy) sum_proj_init += xy_proj[ixy] * init_xy[ixy];
  std::cout << std::format("sum XYProj initXY = ({:.17g}, {:.3g})\n", sum_proj_init.real(), sum_proj_init.imag());

  // phi_s(iz, n), phi_s'(iz, n): layout iz + n_planes n.
  std::vector<cplx> phi_s(std::size_t(n_planes) * n_energies, cplx(0.0, 0.0));
  std::vector<cplx> dphi_s(phi_s.size(), cplx(0.0, 0.0));
  int steps = 0;
  bool finite = true;
  double max_spec_rc = 0.0;
  PropagationObserver observer;
  observer.on_averages = [&](const Averages& a) {
    std::cout << std::format("t = {:9.1f}  norm = {:.10e}  <Z> = {:.6f}  <r> = {:.6f}\n", a.time, a.norm, a.z_mean,
                             a.r_mean)
              << std::flush;
  };
  observer.on_flux_surface = [&](double time, const FluxSurfaceWave&) {
    const SpecularWave& spec = propagation.specular();
    std::vector<cplx> psi(n_planes), dpsi(n_planes);
    for (int iz = 0; iz < n_planes; ++iz) {
      const std::span<const cplx> column = std::span<const cplx>(spec.values).subspan(spec.index(0, iz), nr);
      psi[iz] = spec(ir_c, iz);
      dpsi[iz] = flux_r_derivative(weights, column);
      if (!std::isfinite(std::abs(psi[iz])) || !std::isfinite(std::abs(dpsi[iz]))) finite = false;
      max_spec_rc = std::max(max_spec_rc, std::abs(psi[iz]));
    }
    // phi = phi + c psi in time order, as FluxAccumulator::add_tile.
    for (int n = 0; n < n_energies; ++n) {
      const cplx c = flux_step_coefficient(energies, n, time, time_step);
      for (int iz = 0; iz < n_planes; ++iz) {
        const std::size_t k = iz + std::size_t(n_planes) * n;
        phi_s[k] = phi_s[k] + c * psi[iz];
        dphi_s[k] = dphi_s[k] + c * dpsi[iz];
      }
    }
    steps = steps + 1;
  };
  propagation.run(observer);
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

  const FluxCoefficients& run_coefficients = propagation.flux()->coefficients();
  const FluxGrid grid = make_flux_grid(config, grids, n_rot);
  const FluxResults run_results = propagation.flux_results();
  const FluxResults scat = flux_results(run_coefficients, energies, grid);

  // Embedding: block (ichn, ixy) += initXY(ixy) phi_s. spec_only holds the embedded part alone;
  // scat_channel the scattering wave's own component along initXY in channel ichn (diagnostic).
  FluxCoefficients total = run_coefficients;
  FluxCoefficients spec_only(run_coefficients.n_blocks(), n_planes, n_energies);
  FluxCoefficients scat_channel(run_coefficients.n_blocks(), n_planes, n_energies);
  for (int n = 0; n < n_energies; ++n)
    for (int iz = 0; iz < n_planes; ++iz) {
      const std::size_t k = iz + std::size_t(n_planes) * n;
      cplx proj(0.0, 0.0), dproj(0.0, 0.0);
      for (int ixy = 0; ixy < n_xy; ++ixy) {
        const std::size_t idx = run_coefficients.index(ichn + n_rot * ixy, iz, n);
        proj += xy_proj[ixy] * run_coefficients.phi[idx];
        dproj += xy_proj[ixy] * run_coefficients.phi_derivative[idx];
      }
      for (int ixy = 0; ixy < n_xy; ++ixy) {
        const std::size_t idx = total.index(ichn + n_rot * ixy, iz, n);
        const cplx v = init_xy[ixy] * phi_s[k];
        const cplx dv = init_xy[ixy] * dphi_s[k];
        spec_only.phi[idx] = v;
        spec_only.phi_derivative[idx] = dv;
        total.phi[idx] = total.phi[idx] + v;
        total.phi_derivative[idx] = total.phi_derivative[idx] + dv;
        scat_channel.phi[idx] = init_xy[ixy] * proj;
        scat_channel.phi_derivative[idx] = init_xy[ixy] * dproj;
      }
    }
  const FluxResults spec = flux_results(spec_only, energies, grid);
  const FluxResults tot = flux_results(total, energies, grid);
  const FluxResults chan = flux_results(scat_channel, energies, grid);
  const ScatteringResults s = propagation.scattering_results();

  std::cout << std::format("analysis steps seen: {}; specular values finite: {}; max |spec(r_c)| = {:.3e}\n", steps,
                           finite, max_spec_rc);
  std::cout << std::format("check (i) own P_scat == run flux_results bitwise: {}\n",
                           bitwise_equal(scat.probability, run_results.probability));
  std::cout << std::format("run time {:.1f} s\n", seconds);

  std::ofstream out(output_path);
  out << "n\tE_flux_eV\tE_coll_smatrix_eV\tP_scat\tP_spec\tP_total\tcross\tP_scat_ichn_initXY\tone_minus_sumR\tsumR\n";
  for (int n = 0; n < n_energies; ++n) {
    const double e_s = n < static_cast<int>(s.energies.size()) ? s.energies[n].collision_energy : NAN;
    const double react = n < static_cast<int>(s.energies.size()) ? s.energies[n].p_reaction : NAN;
    const double sum_r = n < static_cast<int>(s.energies.size()) ? s.energies[n].p_total : NAN;
    const double cross = tot.probability[n] - scat.probability[n] - spec.probability[n];
    out << std::format("{}\t{:.10f}\t{:.10f}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\t{:.16e}\n", n,
                       scat.energy[n], e_s, scat.probability[n], spec.probability[n], tot.probability[n], cross,
                       chan.probability[n], react, sum_r);
  }
  std::cout << "S-matrix energies: " << s.energies.size() << "; written " << output_path << '\n';
  return finite ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: spodvr_q7_specular_flux <input.json> <output.tsv>\n";
    return 2;
  }
  try {
    return run(argv[1], argv[2]);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
}

#pragma once

#include <complex>
#include <functional>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "fft/fft.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

struct FluxSurfaceWave;  // analysis/flux_analysis.hpp

// Where the wave packet lives and is propagated: in host memory (OpenMP) or on a GPU (CUDA).
enum class Device { Cpu, Cuda };

// "cpu" / "cuda" (command line --device). Throws InputError for other values.
[[nodiscard]] Device parse_device(std::string_view name);
[[nodiscard]] std::string to_string(Device device);

// Owns the wave packet and the split-operator stages that act on it, on the host or on a GPU.
// The time loop (WavePacketPropagation::run) talks to it once per analysis step: propagate, then
// read the small pieces the analysis needs. Everything that needs all of psi (checkpoint,
// restart, dumps) goes through read_wave_packet() / wave_packet_on_host(), outside the time loop.
// All data handed out are in the host layout of wave_function.hpp.
class PropagationBackend {
 public:
  virtual ~PropagationBackend() = default;

  // Replaces the wave packet (restart): `read` fills a host wave packet of the right shape. The
  // CPU backend passes its own wave packet (read in place, no second copy of psi in memory); a
  // GPU backend passes a host buffer and uploads it afterwards.
  virtual void read_wave_packet(const std::function<void(WavePacket&)>& read) = 0;
  // The complete wave packet in host memory (checkpoint, dump). A GPU backend copies it into a
  // host buffer first, so this is expensive there: only call it outside the time loop.
  [[nodiscard]] virtual const WavePacket& wave_packet_on_host() = 0;

  // [T(dt/2) V(dt) T(dt/2)]^nanspo. Fortran: SplitOperator.f90 Propagate
  virtual void propagate_analysis_step() = 0;

  // The specular wave after the last propagate_analysis_step (analysis surface, BK, averages).
  [[nodiscard]] virtual const SpecularWave& specular() const = 0;
  // The nr (Z, r) blocks of Z index iz, layout (irot, ixy, ir): the input of the BK analysis.
  [[nodiscard]] virtual std::span<const std::complex<double>> z_plane(int iz) const = 0;
  // sum |psi|^2 over (irot, ixy) at every (Z, r) point izr, summed in memory order.
  // Fortran: Averages.f90 complexInnerProd(1, izr)
  [[nodiscard]] virtual std::vector<double> block_densities() const = 0;
  // psi and d psi / dr at the flux surface r = r_c after the last propagate_analysis_step, on the
  // planes iz < nz_potential, into every element of `out` (caller-owned, a slot of the
  // FluxAccumulator; host layout (block, iz)). Throws std::logic_error without flux analysis
  // (doFluxAnalysis) or unless out has nz_potential planes of n_rot nx ny blocks.
  // Fortran: WaveFunction.f90 getWaveFunctionAtZr, getWaveFuncDerivAtZr
  virtual void flux_surface(FluxSurfaceWave& out) const = 0;
};

// Everything a backend is built from. The references must outlive the backend; the potential
// (premultiplied by -dt inside the backend) and the initial specular wave are taken over.
struct BackendInput {
  const Config& config;
  const SpatialGrids& grids;
  const AngularGrid& angular;
  const RotationalBasis& basis;
  const PropagationTables& tables;
  const Fft& fft;  // host FFT library (CPU backend)
};

// Creates the backend of a propagation (see app/backend_factory.cpp for the device choice).
using BackendFactory =
    std::function<std::unique_ptr<PropagationBackend>(const BackendInput&, PotentialGrid, SpecularWave)>;

}  // namespace spodvr

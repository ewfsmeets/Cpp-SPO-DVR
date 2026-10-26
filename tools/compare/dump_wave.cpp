// Runs the propagation of an input file and writes the final wave packet (scattering grid, then
// specular grid) as raw complex<double> values, for comparisons between builds or versions
// with spodvr_compare_wave.
// Usage: spodvr_dump_wave <input.json> <output file>   (threads: OMP_NUM_THREADS)
#include <cstdio>
#include <fstream>

#include <nlohmann/json.hpp>

#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "propagation/propagation.hpp"

using namespace spodvr;

namespace {

void write_values(std::ofstream& out, const ComplexVector& values) {
  for (const auto& c : values) {
    const double parts[2] = {c.real(), c.imag()};
    out.write(reinterpret_cast<const char*>(parts), sizeof(parts));
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: %s <input.json> <output file>\n", argv[0]);
    return 1;
  }
  std::ifstream file(argv[1]);
  const Config config = build_config(parse_input(nlohmann::json::parse(file)));
  const InputParameters& p = config.input;
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular(p.j_max, config.derived.mj_max_basis, p.homonuclear, p.j_in, p.mj_in);
  PotentialGrid potential = compute_potential_grid(*make_potential(p.system_id), config, grids, angular);
  const FftwFft fft;
  const GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver());

  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft);
  propagation.run([](const Averages&) {}, [](double) {});

  std::ofstream out(argv[2], std::ios::binary);
  write_values(out, propagation.wave_packet().scattering.values);
  write_values(out, propagation.wave_packet().specular.values);
  return out ? 0 : 1;
}

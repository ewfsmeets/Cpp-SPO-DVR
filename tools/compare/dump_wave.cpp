// Runs the propagation of an input file on the CPU or a GPU and writes the final wave packet
// (scattering grid, then specular grid) as raw complex<double> values (re, im), for comparisons
// between builds, versions or devices with spodvr_compare_wave (docs/performance.md).
// Usage: spodvr_dump_wave <input.json> <output file> [--device cpu|cuda]
//   --device: default cpu; cuda needs a build with -DSPODVR_USE_CUDA=ON. Threads: OMP_NUM_THREADS.
//   Exit 0 written, 1 error (input, device, file), 2 bad arguments.
#include <complex>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "app/backend_factory.hpp"
#include "core/error.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/binary_file.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "propagation/propagation.hpp"

using namespace spodvr;

namespace {

// Printed on stderr for every bad command line.
constexpr std::string_view usage_line = "usage: spodvr_dump_wave <input.json> <output file> [--device cpu|cuda]\n";

// The command line of the tool.
struct Arguments {
  std::filesystem::path input;
  std::filesystem::path output;
  Device device = Device::Cpu;
};

// Two file names and an optional "--device cpu|cuda", in any order; nullopt if they do not parse.
// An unknown device name throws InputError (parse_device).
std::optional<Arguments> parse_arguments(const std::vector<std::string>& args) {
  Arguments parsed;
  std::vector<std::string> files;
  for (std::size_t i = 0; i < args.size(); ++i) {
    if (args[i] == "--device") {
      if (i + 1 == args.size()) return std::nullopt;
      parsed.device = parse_device(args[++i]);
    } else if (args[i].starts_with("-")) {
      return std::nullopt;
    } else {
      files.push_back(args[i]);
    }
  }
  if (files.size() != 2) return std::nullopt;
  parsed.input = files[0];
  parsed.output = files[1];
  return parsed;
}

// std::complex<double> is laid out as (re, im), so the bytes are the values' real and imaginary
// parts in memory order, scattering grid first.
void write_wave_packet(const std::filesystem::path& path, const WavePacket& wave_packet) {
  BinaryWriter out(path);
  out.write_array<std::complex<double>>(wave_packet.scattering.values);
  out.write_array<std::complex<double>>(wave_packet.specular.values);
  out.close();
}

// The propagation of the input with the backend of `device` (as the program's run mode 0, without
// analysis output), then the dump of the final wave packet.
void dump_final_wave_packet(const Arguments& arguments) {
  // Checked before the (possibly long) set-up: the GPU must be available in this build.
  const std::string device_summary = prepare_device(arguments.device);
  if (!device_summary.empty()) std::cout << device_summary << '\n';

  const Config config = build_config(read_input_file(arguments.input));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  PotentialGrid potential = compute_potential_grid(*make_potential(config.input), config, grids, angular);
  const FftwFft fft;
  const GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver());

  WavePacketPropagation propagation(config, grids, angular, std::move(potential), gas, fft,
                                    backend_factory(arguments.device));
  propagation.run({});
  write_wave_packet(arguments.output, propagation.wave_packet_on_host());
}

}  // namespace

int main(int argc, char** argv) {
  const std::vector<std::string> args(argv + 1, argv + argc);
  std::optional<Arguments> arguments;
  try {
    arguments = parse_arguments(args);
  } catch (const InputError& e) {
    std::cerr << "error: " << e.what() << '\n' << usage_line;
    return 2;
  }
  if (!arguments) {
    std::cerr << usage_line;
    return 2;
  }
  try {
    dump_final_wave_packet(*arguments);
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << '\n';
    return 1;
  }
  return 0;
}

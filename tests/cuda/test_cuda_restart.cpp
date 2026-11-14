// Tests: checkpoint / restart and the command line with the GPU backend (--device cuda). A GPU
// run split by a restart (native or Fortran files) is bitwise identical to the unsplit GPU run,
// also with the flux analysis; a GPU run can continue a CPU run (and the results then agree to rounding).
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <span>
#include <sstream>
#include <string>

#include <unistd.h>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "app/cli.hpp"
#include "app/run.hpp"
#include "cuda_test_setup.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/checkpoint.hpp"

using namespace spodvr;
using namespace spodvr::testing;
namespace fs = std::filesystem;

namespace {

nlohmann::json short_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/cu111_dailight_short_sym/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["CheckpointFormat"] = "both";
  return json;
}

std::string read_file(const fs::path& path) {
  std::ifstream file(path, std::ios::binary);
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

int run(const nlohmann::json& json, const fs::path& dir, const std::string& device, std::string& err) {
  fs::create_directories(dir);
  const fs::path input = dir / "input.json";
  std::ofstream(input) << json.dump(1);
  CommandLine cl;
  cl.input_file = input;
  cl.output_dir = dir;
  cl.threads = 2;
  cl.device = device;
  std::ostringstream out;
  std::ostringstream errors;
  const int code = run_program(cl, out, errors, false);
  err = errors.str();
  return code;
}

nlohmann::json with_times(nlohmann::json json, double t_start, double t_end, int n_steps) {
  json["tStart"] = t_start;
  json["tEnd"] = t_end;
  json["nStep"] = n_steps;
  return json;
}

nlohmann::json restart_from(nlohmann::json json, const fs::path& dir, const char* format = "native") {
  json["NewCalc"] = false;
  json["RestartDir"] = dir.string();
  json["RestartFormat"] = format;
  return json;
}

// short_json() with the flux analysis of the gate case cu111_dailight_short: r_c = 1.7 bohr,
// 0.2 .. 0.6 eV by 0.15 eV (3 energies).
nlohmann::json flux_json() {
  nlohmann::json json = short_json();
  json.update({{"doFluxAnalysis", true},
               {"rFluxAnalysis", 1.7},
               {"EMinFluxAna", 0.2},
               {"EMaxFluxAna", 0.6},
               {"DeltaEFluxAna", 0.15}});
  return json;
}

// The unsplit run (0 .. 80 dt) in root/full, and the same run split at 40 dt in root/first and
// root/second (native files) or root/second_fortran; the first part on `first_device`, the rest on
// the GPU. Every run must succeed.
void run_split(const nlohmann::json& json, const fs::path& root, const std::string& first_device) {
  const double dt = json["dtStep"].get<double>();
  std::string err;
  ASSERT_EQ(run(with_times(json, 0.0, 80 * dt, 80), root / "full", "cuda", err), 0) << err;
  ASSERT_EQ(run(with_times(json, 0.0, 40 * dt, 40), root / "first", first_device, err), 0) << err;
  for (const char* format : {"native", "fortran"}) {
    const fs::path second = root / (std::string(format) == "native" ? "second" : "second_fortran");
    ASSERT_EQ(run(restart_from(with_times(json, 40 * dt, 80 * dt, 40), root / "first", format), second, "cuda", err), 0)
        << format << ": " << err;
  }
}

// phi and phi' of the native flux coefficient file of the run in dir (sizes from its input.json).
FluxCoefficients read_flux_file(const fs::path& dir) {
  const Config config = build_config(read_input_file(dir / "input.json"));
  const SpatialGrids grids = make_spatial_grids(config);
  const int n_rot = RotationalBasis(make_angular_grid(config), config.input.j_in, config.input.mj_in).size();
  FluxCoefficients c(n_rot * config.input.nx * config.input.ny, grids.nz_potential, flux_energy_count(config.input));
  (void)read_flux_coefficients(dir / "flux_coefficients.bin", config, c);
  return c;
}

// The number of values that are zero or not finite.
int zero_or_non_finite(std::span<const std::complex<double>> values) {
  int count = 0;
  for (const std::complex<double>& v : values)
    if (v == 0.0 || !std::isfinite(v.real()) || !std::isfinite(v.imag())) ++count;
  return count;
}

}  // namespace

TEST(CudaRestart, SplitGpuRunIsBitwiseIdenticalToTheUnsplitGpuRun) {
  SPODVR_REQUIRE_GPU();
  const fs::path root = fs::temp_directory_path() / ("spodvr_cuda_restart_" + std::to_string(::getpid()));
  fs::remove_all(root);
  run_split(short_json(), root, "cuda");
  for (const char* name : {"checkpoint.bin", "bk_coefficients.bin", "s_matrix.tsv", "probabilities_rovib.tsv",
                           "probabilities_diffraction.tsv", "fort.9", "fort.21"}) {
    EXPECT_EQ(read_file(root / "full" / name), read_file(root / "second" / name)) << name;
    EXPECT_EQ(read_file(root / "full" / name), read_file(root / "second_fortran" / name)) << name;
  }
  fs::remove_all(root);
}

// With the flux analysis the GPU restart continues phi and phi' from the flux coefficient files of
// either format: the flux files of the split GPU run are bitwise those of the unsplit GPU run.
TEST(CudaRestart, SplitGpuFluxRunIsBitwiseIdenticalToTheUnsplitGpuRun) {
  SPODVR_REQUIRE_GPU();
  const fs::path root = fs::temp_directory_path() / ("spodvr_cuda_restart_flux_" + std::to_string(::getpid()));
  fs::remove_all(root);
  run_split(flux_json(), root, "cuda");
  // The comparison sees the restored part: every coefficient of the first part is non-zero and finite.
  const FluxCoefficients first = read_flux_file(root / "first");
  EXPECT_EQ(zero_or_non_finite(first.phi), 0);
  EXPECT_EQ(zero_or_non_finite(first.phi_derivative), 0);
  for (const char* name : {"flux_coefficients.bin", "flux_coeffs.dat", "flux_coeffs_deriv.dat", "flux_reaction.tsv",
                           "flux_density.tsv", "checkpoint.bin", "bk_coefficients.bin", "fort.9", "fort.21"}) {
    const std::string full = read_file(root / "full" / name);
    EXPECT_FALSE(full.empty()) << name;
    EXPECT_EQ(full, read_file(root / "second" / name)) << name;
    EXPECT_EQ(full, read_file(root / "second_fortran" / name)) << name;
  }
  fs::remove_all(root);
}

// A GPU run continues a CPU run, its flux analysis included (native and Fortran files): the flux
// coefficients agree with those of the unsplit GPU run to the rounding of CPU against GPU. Steps of
// 10 a.u. (0 .. 800) bring the packet to the coupling region, max |phi| 1.2e-6. Measured 2026-11-13:
// 2.1e-10 of max |phi| (phi and phi'); with the restore skipped 7.2e-5 (phi) and 7.4e-5 (phi'), the
// first part's share. The tolerance 1e-8 lies between. (With steps of 2 a.u. max |phi| is 2e-14, and
// psi's rounding, relative to max |psi| = 1.6, makes the CPU-GPU difference 3e-3 of it.)
TEST(CudaRestart, GpuRunContinuesTheFluxAnalysisOfACpuRun) {
  SPODVR_REQUIRE_GPU();
  const fs::path root = fs::temp_directory_path() / ("spodvr_cuda_restart_cpu_" + std::to_string(::getpid()));
  fs::remove_all(root);
  nlohmann::json json = flux_json();
  json["dtStep"] = 10.0;
  json["dtAna"] = 10.0;
  run_split(json, root, "cpu");
  const FluxCoefficients gpu = read_flux_file(root / "full");
  // relative_difference falls back to the absolute difference for all-zero coefficients: exclude them.
  ASSERT_EQ(zero_or_non_finite(gpu.phi), 0);
  ASSERT_EQ(zero_or_non_finite(gpu.phi_derivative), 0);
  for (const char* second : {"second", "second_fortran"}) {
    const FluxCoefficients cpu_then_gpu = read_flux_file(root / second);
    EXPECT_LE(relative_difference(gpu.phi, cpu_then_gpu.phi), 1e-8) << second;
    EXPECT_LE(relative_difference(gpu.phi_derivative, cpu_then_gpu.phi_derivative), 1e-8) << second;
  }
  fs::remove_all(root);
}

TEST(CudaRestart, UnknownDeviceIsAnError) {
  EXPECT_THROW((void)parse_command_line({"input.json", "--device", "tpu"}), std::invalid_argument);
  EXPECT_EQ(parse_command_line({"input.json", "--device", "cuda"}).device, "cuda");
  EXPECT_EQ(parse_command_line({"input.json"}).device, "cpu");
}

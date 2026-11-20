// Tests: the energy-resolved sums (EnergyResolvedSites, EnergyResolvedEnergies): the accumulator against the sum
// written out by hand (layout (izr, irot) -> [Z, r, rot], phase factors, specular sums), independence of the thread
// count, the files' round trip through a restart's reader and its checks, and the memory estimate's item against the
// sums the program allocates.
#include <gtest/gtest.h>

#include <complex>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <omp.h>

#include "analysis/energy_resolved.hpp"
#include "core/constants.hpp"
#include "core/error.hpp"
#include "counting_backend.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "input/validate_outputs.hpp"
#include "propagation/cpu_backend.hpp"
#include "propagation/memory_estimate.hpp"
#include "report/energy_resolved_files.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using namespace spodvr::test_support;
using ::test_support::TemporaryDirectory;

namespace {

using cplx = std::complex<double>;

// Hand-made sums of 2 sites, 3 energies, nZ 3, nR 4, 2 rotational states, 5 specular Z points.
EnergyResolvedSums small_sums() {
  EnergyResolvedSums s;
  s.site_names = {"a", "b"};
  s.site_ixy = {0, 5};
  s.collision_ev = {0.3, 0.4, 0.5};
  s.total_hartree = {0.011, 0.0147, 0.0184};
  s.n_z = 3;
  s.n_r = 4;
  s.n_rot = 2;
  s.n_z_specular = 5;
  s.analysis_time_step = 10.0;
  s.sums.assign(2 * 3 * s.site_block(), cplx(0.0, 0.0));
  s.specular_sums.assign(3 * s.specular_block(), cplx(0.0, 0.0));
  return s;
}

// A made-up wave packet value at (ixy, izr, irot) and time t, different everywhere.
cplx psi_value(int ixy, std::size_t izr, std::size_t irot, double t) {
  const double z = double(izr), rot = double(irot);
  return {std::sin(0.1 * ixy + 0.01 * z + 0.3 * rot + 1e-3 * t), std::cos(0.2 * ixy - 0.02 * z + 1e-3 * t)};
}

// The column reader of a made-up wave packet at time t: (izr, irot) with izr fastest, n_zr values per rotation.
EnergyResolvedAccumulator::ColumnReader reader_at(double t, std::size_t n_zr, std::size_t n_rot) {
  return [=](int ixy, std::span<cplx> out) {
    for (std::size_t irot = 0; irot < n_rot; ++irot)
      for (std::size_t izr = 0; izr < n_zr; ++izr) out[izr + n_zr * irot] = psi_value(ixy, izr, irot, t);
  };
}

SpecularWave specular_at(double t, int n_r, int n_z) {
  SpecularWave chi(n_r, n_z);
  for (std::size_t k = 0; k < chi.values.size(); ++k) {
    const double x = double(k);
    chi.values[k] = cplx(1e-2 * x, std::sin(1e-3 * t + x));
  }
  return chi;
}

// Adds the steps t = 10, 20, 30 to the sums with the given thread count.
EnergyResolvedSums summed(int threads) {
  const int before = omp_get_max_threads();
  omp_set_num_threads(threads);
  EnergyResolvedAccumulator accumulator(small_sums());
  for (const double t : {10.0, 20.0, 30.0}) accumulator.add(t, reader_at(t, 12, 2), specular_at(t, 4, 5));
  omp_set_num_threads(before);
  return accumulator.sums();
}

}  // namespace

// Every element equals the sum written out by hand, in the same operation order (so bitwise): site isite, energy ien,
// (iz, ir, irot) in [Z, r, rot] order, from column element izr = ir + nR iz of rotation irot; the specular sums in
// the SpecularWave's (ir, iz) order. The time and the step count are those of the last step.
TEST(EnergyResolved, AccumulatorEqualsTheSumByHand) {
  const EnergyResolvedSums s = summed(1);
  const EnergyResolvedSums layout = small_sums();
  int checked = 0;
  for (std::size_t isite = 0; isite < 2; ++isite)
    for (std::size_t ien = 0; ien < 3; ++ien)
      for (std::size_t iz = 0; iz < 3; ++iz)
        for (std::size_t ir = 0; ir < 4; ++ir)
          for (std::size_t irot = 0; irot < 2; ++irot) {
            cplx expected(0.0, 0.0);
            for (const double t : {10.0, 20.0, 30.0}) {
              const cplx phase = std::exp(Constants::imaginary_unit * layout.total_hartree[ien] * t);
              expected = expected + psi_value(layout.site_ixy[isite], ir + 4 * iz, irot, t) * phase * 10.0;
            }
            const std::size_t k = (isite * 3 + ien) * layout.site_block() + irot + 2 * (ir + 4 * iz);
            EXPECT_EQ(s.sums[k], expected) << isite << " " << ien << " " << iz << " " << ir << " " << irot;
            ++checked;
          }
  EXPECT_EQ(checked, 144);
  for (std::size_t ien = 0; ien < 3; ++ien)
    for (std::size_t k = 0; k < 20; ++k) {
      cplx expected(0.0, 0.0);
      for (const double t : {10.0, 20.0, 30.0}) {
        const cplx phase = std::exp(Constants::imaginary_unit * layout.total_hartree[ien] * t);
        expected = expected + specular_at(t, 4, 5).values[k] * phase * 10.0;
      }
      EXPECT_EQ(s.specular_sums[ien * 20 + k], expected) << ien << " " << k;
    }
  EXPECT_EQ(s.last_time, 30.0);
  EXPECT_EQ(s.steps_summed, 3);
}

// One owner per element: 1 and 4 threads give bitwise the same sums.
TEST(EnergyResolved, ThreadCountDoesNotChangeTheSums) {
  const EnergyResolvedSums one = summed(1);
  const EnergyResolvedSums four = summed(4);
  EXPECT_EQ(one.sums, four.sums);
  EXPECT_EQ(one.specular_sums, four.specular_sums);
}

// The files hold the sums a restart reads back bitwise; other energies, other sites or sums that end at another time
// than tStart are refused.
TEST(EnergyResolved, FilesRoundTripThroughTheRestartReader) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_short/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["PlotSites"] = {{{"name", "a"}, {"x", 0.0}, {"y", 0.0}}, {{"name", "b"}, {"x", 0.5}, {"y", 0.0}}};
  json["EnergyResolvedSites"] = {"a", "b"};
  json["EnergyResolvedEnergies"] = {0.3, 0.4, 0.5};
  const Config config = build_config(parse_input(json));
  const SpatialGrids grids = make_spatial_grids(config);
  EnergyResolvedSums sums = summed(1);
  sums.site_ixy = {plot_site_ixy(config.input.plot_sites[0], config.input.nx, config.input.ny),
                   plot_site_ixy(config.input.plot_sites[1], config.input.nx, config.input.ny)};
  sums.n_z = config.input.nz;  // the shape of this run, filled with a pattern
  sums.n_r = config.input.nr;
  sums.n_z_specular = static_cast<int>(grids.z_specular.size());
  sums.sums.assign(2 * 3 * sums.site_block(), cplx(0.0, 0.0));
  sums.specular_sums.assign(3 * sums.specular_block(), cplx(0.0, 0.0));
  for (std::size_t k = 0; k < sums.sums.size(); ++k) sums.sums[k] = cplx(double(k), -0.5 * double(k));
  for (std::size_t k = 0; k < sums.specular_sums.size(); ++k) sums.specular_sums[k] = cplx(-double(k), 2.0);
  sums.last_time = 600.0;
  const std::vector<cplx> s2s(std::size_t(config.input.nx * config.input.ny), cplx(0.5, 0.25));

  TemporaryDirectory dir("energy_resolved_files");
  EnergyResolvedFiles files(dir.path(), config, grids, sums, s2s);
  files.write_partial(sums);
  files.write_final(sums);

  EnergyResolvedSums expected = sums;
  std::fill(expected.sums.begin(), expected.sums.end(), cplx(0.0, 0.0));
  std::fill(expected.specular_sums.begin(), expected.specular_sums.end(), cplx(0.0, 0.0));
  Config restart = config;
  restart.input.t_start = 600.0;
  const EnergyResolvedSums read = read_energy_resolved_sums(dir.path(), restart, expected, 3);
  EXPECT_EQ(read.sums, sums.sums);
  EXPECT_EQ(read.specular_sums, sums.specular_sums);
  EXPECT_EQ(read.steps_summed, 3);
  EXPECT_EQ(read.last_time, 600.0);
  EXPECT_NO_THROW(check_energy_resolved_restart(dir.path(), restart.input));

  Config earlier_start = restart;
  earlier_start.input.t_start = 500.0;  // the sums end at 600
  EXPECT_THROW((void)read_energy_resolved_sums(dir.path(), earlier_start, expected, 3), CheckpointError);
  EXPECT_THROW((void)read_energy_resolved_sums(dir.path(), restart, expected, 4), CheckpointError);  // another run's
  Config other_energy = restart;
  other_energy.input.energy_resolved_energies[2] = 0.51;
  EXPECT_THROW(check_energy_resolved_restart(dir.path(), other_energy.input), InputError);
  EXPECT_THROW((void)read_energy_resolved_sums(dir.path(), other_energy, expected, 3), InputError);
  Config moved_site = restart;
  moved_site.input.plot_sites[1].x = 0.0;  // "b" at another grid point, the same name
  moved_site.input.plot_sites[1].y = 0.5;
  EXPECT_THROW(check_energy_resolved_restart(dir.path(), moved_site.input), InputError);
  Config renamed = restart;
  renamed.input.plot_sites[1].name = "c";
  renamed.input.energy_resolved_sites[1] = "c";
  EXPECT_THROW(check_energy_resolved_restart(dir.path(), renamed.input), InputError);
}

// A run without analysis steps (nStep = 0) writes its (empty) sums as both the partial and the final ones, and its
// checkpoint.
TEST(EnergyResolved, RunWithoutAnalysisStepsWritesItsSums) {
  TemporaryDirectory dir("energy_resolved_no_steps");
  nlohmann::json json = small_hardwall_input();
  json["PlotSites"] = {{{"name", "a"}, {"x", 0.0}, {"y", 0.0}}};
  json["EnergyResolvedSites"] = {"a"};
  json["EnergyResolvedEnergies"] = {0.3};
  json["tEnd"] = 0.0;
  json["nStep"] = 0;
  BackendCallCounts counts;
  ASSERT_EQ(run_counted(json, dir.path(), make_cpu_backend, counts), 0);
  EXPECT_EQ(counts.site_column, 0);
  EXPECT_TRUE(std::filesystem::exists(dir.path() / "checkpoint.bin"));
  EXPECT_EQ(std::filesystem::file_size(dir.path() / "energy_resolved" / "partial_time.bin"), 8u);
  EXPECT_EQ(std::filesystem::file_size(dir.path() / "energy_resolved" / "steps_summed.bin"), 4u);
}

// The memory estimate's item equals what the program allocates: the sums of a new run plus one site's column.
TEST(EnergyResolved, MemoryEstimateItemMatchesTheAllocation) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_short/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["PlotSites"] = {{{"name", "a"}, {"x", 0.0}, {"y", 0.0}}};
  json["EnergyResolvedSites"] = {"a"};
  json["EnergyResolvedEnergies"] = {0.3, 0.5};
  const Config config = build_config(parse_input(json));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  const MemoryEstimate estimate = estimate_memory(config, grids, angular, 1);
  ASSERT_EQ(estimate.items.back().name, "Energy-resolved sums (EnergyResolvedSites)");

  const GasPhase gas_stub{RovibrationalStates(1, 0)};  // E(vin, jin) does not change the sizes
  const EnergyResolvedSums start =
      empty_energy_resolved_sums(config, gas_stub, basis.size(), static_cast<int>(grids.z_specular.size()), 0.0);
  const std::size_t column = start.site_block();
  EXPECT_EQ(estimate.items.back().values, start.sums.size() + start.specular_sums.size() + column);
  EXPECT_EQ(estimate.items.back().bytes, estimate.items.back().values * sizeof(cplx));
}

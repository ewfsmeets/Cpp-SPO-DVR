// Tests: frames/ (dtFrames; frames/frame_reducers.*, frames/frame_source.*, report/frame_files.*): the reductions
// of one plane against a hand computation, with a negative test; the same bytes with 1 and 8 threads; the
// projected size; the memory-estimate item. The physics checks (sums, orthogonality, numpy against the
// checkpoints, restart joins) run on the plot fixtures (tools/plots/tests/test_runs_frames.py).
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include <omp.h>

#include "comparison.hpp"
#include "counting_backend.hpp"
#include "fft/fftw_fft.hpp"
#include "frames/frame_extras.hpp"
#include "frames/frame_source.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/associated_legendre.hpp"
#include "grids/gauss_legendre.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/binary_file.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation_tables.hpp"
#include "report/frame_files.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using namespace spodvr::test_support;
using ::test_support::TemporaryDirectory;
namespace fs = std::filesystem;
using Complex = std::complex<double>;

namespace {

std::string file_bytes(const fs::path& path) {
  std::ifstream file(path, std::ios::binary);
  std::ostringstream content;
  content << file.rdbuf();
  return content.str();
}

// The small hard-wall input of the call tests with frames every dtAve (100 a.u.).
nlohmann::json frames_input() {
  nlohmann::json json = small_hardwall_input();
  json["dtFrames"] = 100.0;
  return json;
}

// Runs the program on `json` with `threads` OpenMP threads into `dir`.
std::string run_into(const nlohmann::json& json, const fs::path& dir, int threads) {
  fs::create_directories(dir);
  std::ofstream(dir / "input.json") << json.dump(1);
  CommandLine cl;
  cl.input_file = dir / "input.json";
  cl.output_dir = dir;
  cl.threads = threads;
  std::ostringstream out;
  std::ostringstream err;
  const int previous = omp_get_max_threads();
  const int status = run_program(cl, out, err, false);
  omp_set_num_threads(previous);
  EXPECT_EQ(status, 0) << err.str();
  return out.str();
}

// A small layout: 3 rotational states (the initial one is 1), 2 x 2 (x, y) points, 2 r points, 3 + 2 Z points.
struct SyntheticFrame {
  std::vector<Complex> s2s = {{0.5, 0.1}, {-0.2, 0.4}, {0.3, -0.3}, {0.1, 0.2}};
  FrameLayout layout{3, 5, 2, 4, 3, 1, 0.25, 0.5, s2s};
  SpecularWave specular{2, 5};
  std::vector<Complex> plane = std::vector<Complex>(layout.plane().size());

  SyntheticFrame() {
    for (std::size_t k = 0; k < plane.size(); ++k) plane[k] = {0.1 * double(k % 7) - 0.3, 0.05 * double(k % 5)};
    for (int iz = 0; iz < 5; ++iz)
      for (int ir = 0; ir < 2; ++ir) specular(ir, iz) = {0.2 + 0.1 * iz, -0.1 * ir + 0.05 * iz};
    specular(1, 2) = {0.0, 0.0};  // one point without chi
  }
};

}  // namespace

// One plane by hand: prob_zrxy = |psi_scat + s2s chi (initial state only)|^2 dZ dr dx dy sin, pop_jmj_z its sum
// over r and (x, y) per state, the tail |chi|^2 dZ dr beyond the scattering grid (|chi|^2 as abs(chi)^2, as the
// densities' specular cells, bitwise). Negative test: chi left out of
// the initial state (a reducer that forgets psi_full) changes prob_zrxy at exactly the (r, x, y) points of the
// plane where chi is not zero.
TEST(FrameReducers, OnePlaneByHand) {
  SyntheticFrame f;
  const int iz = 2;
  FrameRecord record = make_frame_record(f.layout);
  std::vector<double> scratch(std::size_t(f.layout.nr) * f.layout.n_rot);
  reduce_scattering_plane(f.layout, iz, f.plane, f.specular, scratch, record);
  reduce_specular_tail(f.layout, f.specular, record);
  std::vector<double> pop(3, 0.0);
  for (int ir = 0; ir < 2; ++ir)
    for (int ixy = 0; ixy < 4; ++ixy) {
      double cell = 0.0;
      for (int irot = 0; irot < 3; ++irot) {
        Complex psi = f.plane[f.layout.plane()(irot, ixy, ir)];
        if (irot == 1) psi += f.s2s[ixy] * f.specular(ir, iz);
        cell += std::norm(psi);
        pop[irot] += std::norm(psi) * 0.25;
      }
      EXPECT_NEAR(record.prob_zrxy[f.layout.zrxy()(ixy, ir, iz)], cell * 0.25, 1e-15) << ir << " " << ixy;
    }
  for (int irot = 0; irot < 3; ++irot) EXPECT_NEAR(record.pop_jmj_z[f.layout.jmj_z()(irot, iz)], pop[irot], 1e-15);
  for (int izt = 3; izt < 5; ++izt)
    for (int ir = 0; ir < 2; ++ir)
      EXPECT_EQ(record.prob_tail_zr[f.layout.tail_zr()(ir, izt - 3)],
                std::abs(f.specular(ir, izt)) * std::abs(f.specular(ir, izt)) * 0.5);  // as densities/ squares

  FrameLayout slip = f.layout;
  slip.initial_channel = -1;  // chi never added
  FrameRecord without_chi = make_frame_record(slip);
  reduce_scattering_plane(slip, iz, f.plane, f.specular, scratch, without_chi);
  int changed = 0;
  for (int ir = 0; ir < 2; ++ir)
    for (int ixy = 0; ixy < 4; ++ixy) {
      const std::size_t k = f.layout.zrxy()(ixy, ir, iz);
      changed += (std::abs(without_chi.prob_zrxy[k] - record.prob_zrxy[k]) <= 1e-15) ? 0 : 1;
    }
  EXPECT_EQ(changed, 4);  // r index 0 only: chi(1, 2) = 0
}

// The reductions do not depend on the thread count: a run with 1 and with 8 threads writes the same bytes.
TEST(FrameFiles, OneAndEightThreadsWriteTheSameBytes) {
  TemporaryDirectory dir("frames_threads");
  run_into(frames_input(), dir.path() / "one", 1);
  const std::string out = run_into(frames_input(), dir.path() / "eight", 8);
  EXPECT_NE(out.find("Frames (dtFrames): 11 written"), std::string::npos) << out;
  for (const std::string name :
       {"manifest.json", "times.bin", "prob_zrxy.bin", "prob_tail_zr.bin", "pop_jmj_z.bin", "surface_index.bin"}) {
    const std::string one = file_bytes(dir.path() / "one" / "frames" / name);
    EXPECT_FALSE(one.empty()) << name;
    EXPECT_EQ(one, file_bytes(dir.path() / "eight" / "frames" / name)) << name;
  }
}

// The projected size is the written size: 11 records (tStart and every 100 a.u. to 1000) of the datasets and the
// time; the PLOT DATA line names the key.
TEST(FrameFiles, ProjectedSizeIsTheWrittenSize) {
  TemporaryDirectory dir("frames_size");
  const std::string out = run_into(frames_input(), dir.path(), 8);
  EXPECT_NE(out.find("dtFrames"), std::string::npos);
  const Config config = build_config(parse_input(frames_input()));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  const auto size = frames_data_size(config, grids, angular, basis);
  ASSERT_TRUE(size.has_value());
  std::uintmax_t written = 0;
  for (const auto& entry : fs::directory_iterator(dir.path() / "frames"))
    if (entry.path().extension() == ".bin") written += fs::file_size(entry.path());
  EXPECT_EQ(size->bytes, written);
}

// With dtFrames the memory estimate lists the frame record and its per-r sums: the sizes FrameSource allocates.
TEST(FrameFiles, MemoryEstimateItem) {
  const Config config = build_config(parse_input(frames_input()));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  const MemoryEstimate e = estimate_memory(config, grids, angular, make_grid_symmetry(config, angular), 1);
  ASSERT_EQ(e.items.back().name, "Frame record (dtFrames)");
  const FrameLayout layout = frame_layout(config, grids, basis, {});
  const FrameRecord record = make_frame_record(layout);
  EXPECT_EQ(e.items.back().values, record.prob_zrxy.size() + record.prob_tail_zr.size() + record.pop_jmj_z.size() +
                                       std::size_t(layout.nr) * layout.n_rot);
  const Config without = build_config(parse_input(small_hardwall_input()));
  const MemoryEstimate base = estimate_memory(without, grids, angular, make_grid_symmetry(without, angular), 1);
  EXPECT_EQ(base.items.size() + 1, e.items.size());
}

namespace {

// The small hard-wall input with every FrameExtras extra and densities/ at the frames' times (100 a.u.).
nlohmann::json extras_input() {
  nlohmann::json json = frames_input();
  json["FrameExtras"] = nlohmann::json::array({"angular", "diffraction", "vibration", "rotation_moments", "momentum"});
  json["writeDensities"] = true;
  return json;
}

// Every value of a little-endian binary file of doubles.
std::vector<double> doubles(const fs::path& path) {
  std::vector<double> values(fs::file_size(path) / sizeof(double));
  BinaryReader reader(path);
  reader.read_array(std::span(values));
  return values;
}

// The largest |sum of `extra` over a record's Z row - the same row of `reference`|, relative to the record's largest
// row sum: records x nz rows of `extra_per_z` and `reference_per_z` values each.
double worst_row_difference(const std::vector<double>& extra, std::size_t extra_per_z,
                            const std::vector<double>& reference, std::size_t reference_per_z, std::size_t rows) {
  std::vector<double> a(rows, 0.0);
  std::vector<double> b(rows, 0.0);
  for (std::size_t i = 0; i < extra.size(); ++i) a[i / extra_per_z] += extra[i];
  for (std::size_t i = 0; i < reference.size(); ++i) b[i / reference_per_z] += reference[i];
  double scale = 0.0;
  for (const double v : b) scale = spodvr::testing::nan_aware_max(scale, std::abs(v));
  double worst = 0.0;
  for (std::size_t k = 0; k < rows; ++k) worst = spodvr::testing::nan_aware_max(worst, std::abs(a[k] - b[k]));
  return worst / scale;
}

}  // namespace

// FrameExtras on the small hard-wall run, every record: at each Z the probability over the DVR nodes (angular), over
// the parallel momenta (diffraction: Parseval) and over all gas-phase states (vibration: completeness) is the frame's
// probability at that Z (prob_zrxy summed over r, x, y) to 1e-13 of the largest; prob_kr sums to the frame's
// scattering-grid probability, prob_kz_spec to the densities' specular cells (both Parseval). The hard wall keeps
// J = 0: the rotational moments vanish. Negative: one angular value off by 1e-10 of the largest row fails.
TEST(FrameExtras, EachExtraSumsToTheFramesProbability) {
  TemporaryDirectory dir("frame_extras");
  const std::string out = run_into(extras_input(), dir.path(), 8);
  EXPECT_NE(out.find("Frame extras (FrameExtras), seconds per frame: angular "), std::string::npos) << out;
  const Config config = build_config(parse_input(extras_input()));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  const std::size_t records = 11;
  const std::size_t nz = grids.z_scattering.size();
  const std::size_t nr = config.input.nr;
  const std::size_t nxy = std::size_t(config.input.nx) * config.input.ny;
  const fs::path frames = dir.path() / "frames";
  const std::vector<double> cells = doubles(frames / "prob_zrxy.bin");
  ASSERT_EQ(cells.size(), records * nz * nr * nxy);
  const std::vector<double> angle = doubles(frames / "prob_angle_z.bin");
  const std::size_t nodes = std::size_t(angular.n_theta()) * angular.n_phi();
  EXPECT_LE(worst_row_difference(angle, nodes, cells, nr * nxy, records * nz), 1e-13);
  EXPECT_LE(worst_row_difference(doubles(frames / "pop_nm_z.bin"), nxy, cells, nr * nxy, records * nz), 1e-13);
  const std::size_t n_j = basis_j_values(basis).size();
  EXPECT_LE(worst_row_difference(doubles(frames / "pop_vj_z.bin"), nr * n_j, cells, nr * nxy, records * nz), 1e-13);
  EXPECT_LE(worst_row_difference(doubles(frames / "prob_kr.bin"), nr, cells, nz * nr * nxy, records), 1e-13);
  const std::vector<double> specular = doubles(dir.path() / "densities" / "prob_spec_zr.bin");
  const std::size_t nzsp = grids.z_specular.size();
  EXPECT_LE(worst_row_difference(doubles(frames / "prob_kz_spec.bin"), nzsp, specular, nzsp * nr, records), 1e-13);
  for (const std::string name : {"p2_numerator_z.bin", "rot_numerator_z.bin"}) {
    const std::vector<double> moments = doubles(frames / name);
    ASSERT_EQ(moments.size(), records * nz) << name;
    for (const double v : moments) EXPECT_LE(std::abs(v), 1e-14) << name;  // J = 0 throughout
  }
  std::vector<double> wrong = angle;
  wrong[5 * nz * nodes + 3] += 1e-10;
  EXPECT_GT(worst_row_difference(wrong, nodes, cells, nr * nxy, records * nz), 1e-13);
}

// The extras do not depend on the thread count, and their projected size and memory item are what the run writes
// and allocates.
TEST(FrameExtras, ThreadsSizeAndMemory) {
  TemporaryDirectory dir("frame_extras_threads");
  run_into(extras_input(), dir.path() / "one", 1);
  run_into(extras_input(), dir.path() / "eight", 8);
  for (const std::string name : {"prob_angle_z.bin", "pop_nm_z.bin", "pop_vj_z.bin", "p2_numerator_z.bin",
                                 "rot_numerator_z.bin", "prob_kr.bin", "prob_kz_spec.bin", "manifest.json"}) {
    const std::string one = file_bytes(dir.path() / "one" / "frames" / name);
    EXPECT_FALSE(one.empty()) << name;
    EXPECT_EQ(one, file_bytes(dir.path() / "eight" / "frames" / name)) << name;
  }
  const Config config = build_config(parse_input(extras_input()));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  std::uintmax_t written = 0;
  for (const auto& entry : fs::directory_iterator(dir.path() / "eight" / "frames"))
    if (entry.path().extension() == ".bin") written += fs::file_size(entry.path());
  EXPECT_EQ(frames_data_size(config, grids, angular, basis)->bytes, written);
  const MemoryEstimate e = estimate_memory(config, grids, angular, make_grid_symmetry(config, angular), 4);
  const FrameLayout layout = frame_layout(config, grids, basis, {});
  const std::size_t expected = FrameExtras::memory_values(layout, angular, FrameExtraSwitches::of(config.input),
                                                          static_cast<int>(basis_j_values(basis).size()), 4);
  int found = 0;
  for (const MemoryItem& item : e.items)
    if (item.name == "Frame extras (FrameExtras)") {
      ++found;
      EXPECT_EQ(item.values, expected);
    }
  EXPECT_EQ(found, 1);
}

// The angular extra transforms the nX nY columns of a (Z, r) point in passes of columns_per_pass, as the potential
// step: on 4 x 4 (x, y) points with JMax 28 (n_theta 15, n_phi 60: 4 passes of 4 columns) its record equals, bitwise,
// the reduction of all 16 columns in one pass written out below (the additions in the same order). Negative: one
// plane value changed by 1e-12 of its modulus changes the record.
TEST(FrameExtras, AngularPassesGiveTheBitsOfOnePass) {
  nlohmann::json json = frames_input();
  json["nX"] = 4;
  json["nY"] = 4;
  json["JMax"] = 28;
  json["MJMax"] = 28;
  json["FrameExtras"] = nlohmann::json::array({"angular"});
  const InputParameters input = parse_input(json);
  const std::unique_ptr<Potential> pes = make_potential(input);
  const Config config = build_config(input, *pes);
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  ASSERT_EQ(AngularTransform::columns_per_pass(16, angular.n_theta(), angular.n_phi()), 4);
  const std::vector<double> asymptotic = compute_asymptotic_potential(*pes, config, grids);
  const PropagationTables tables = make_propagation_tables(config, grids, basis, asymptotic);
  const FrameLayout layout = frame_layout(config, grids, basis, tables.specular_to_scattering);
  const FftwFft fft;
  const GasPhase gas = compute_gas_phase(config, grids, angular, asymptotic, fft, LapackEigensolver());
  FrameExtras extras(config, grids, angular, basis, gas, tables, fft, layout);

  // A plane of psi_full (layout (irot, ixy, ir)): cos and sin of the index.
  const std::size_t column = std::size_t(layout.n_rot) * layout.nxy;
  std::vector<std::complex<double>> plane(column * layout.nr);
  for (std::size_t i = 0; i < plane.size(); ++i) plane[i] = {std::cos(0.1 * double(i)), std::sin(0.3 * double(i))};
  const std::size_t nodes = std::size_t(angular.n_theta()) * angular.n_phi();
  // The estimate counts what the extra allocates: the plane of psi_full, the dataset and the per-r sums, and per
  // thread the DVR blocks of one column pass of the transform.
  ASSERT_EQ(extras.angular_columns(), 4);
  const std::size_t dvr_block = AngularTransform::dvr_block_size(angular.n_theta(), angular.n_phi());
  EXPECT_EQ(FrameExtras::memory_values(layout, angular, FrameExtraSwitches::of(config.input),
                                       static_cast<int>(basis_j_values(basis).size()), 1),
            2 * plane.size() + std::size_t(layout.nz) * nodes + std::size_t(layout.nr) * nodes +
                2 * dvr_block * std::size_t(extras.angular_columns()));
  auto reduced = [&](const std::vector<std::complex<double>>& values) {
    FrameRecord record = make_frame_record(layout);
    record.prob_angle_z.assign(std::size_t(layout.nz) * nodes, 0.0);
    extras.reduce_plane(0, values, record);
    extras.reduce_plane(1, values, record);
    return record.prob_angle_z;
  };
  const std::vector<double> record = reduced(plane);

  // The reduction before the column passes: all nX nY columns in one transform, sums over ixy, then r, ascending.
  const AngularTransform one_pass(angular, basis, fft, layout.nxy);
  const int block = one_pass.dvr_block_size();
  const int ld = one_pass.leading_dimension();
  std::vector<std::complex<double>> dvr(std::size_t(block) * layout.nxy);
  std::vector<double> per_r(std::size_t(layout.nr) * nodes, 0.0);
  for (int ir = 0; ir < layout.nr; ++ir) {
    one_pass.fbr_to_dvr(std::span<const std::complex<double>>(plane).subspan(column * ir, column), dvr);
    for (int ixy = 0; ixy < layout.nxy; ++ixy)
      for (int it = 0; it < angular.n_theta(); ++it)
        for (int ip = 0; ip < angular.n_phi(); ++ip)
          per_r[std::size_t(ir) * nodes + ip + angular.n_phi() * it] +=
              std::norm(dvr[std::size_t(block) * ixy + ip + ld * it]);
  }
  std::vector<double> expected(nodes);
  for (std::size_t node = 0; node < nodes; ++node) {
    double sum = 0.0;
    for (int ir = 0; ir < layout.nr; ++ir) sum += per_r[std::size_t(ir) * nodes + node];
    expected[node] = sum * layout.cell_volume;
  }
  for (int iz = 0; iz < 2; ++iz)
    EXPECT_TRUE(std::equal(expected.begin(), expected.end(), record.begin() + std::ptrdiff_t(iz * nodes)))
        << "Z plane " << iz;

  std::vector<std::complex<double>> changed = plane;
  changed[5] *= 1.0 + 1e-12;
  EXPECT_NE(reduced(changed), record);
}

// FrameExtras on a coupled case (cu111_dailight_small on 4 x 4 (x, y) points, JMax 2, to t = 3000 a.u., when the
// packet is at the surface and J = 2 is populated): at every frame and Z the vibrational populations summed over v
// are the rotational populations summed over the states of each J (completeness per J: catches a mixed-up v/J
// layout, which the J = 0 hard wall cannot). Negative: one value off by 10 times the tolerance fails exactly one
// (frame, Z, J).
TEST(FrameExtras, VibrationPerJIsTheRotationalPopulation) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/cu111_dailight_small/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json.update({{"nX", 4},
               {"nY", 4},
               {"JMax", 2},
               {"MJMax", 2},
               {"tEnd", 3000.0},
               {"nStep", 1500},
               {"dtFrames", 1000.0},
               {"FrameExtras", nlohmann::json::array({"vibration"})}});
  TemporaryDirectory dir("frame_extras_coupled");
  run_into(json, dir.path(), 8);
  const Config config = build_config(parse_input(json));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  const std::vector<int> j_values = basis_j_values(basis);
  const std::size_t nz = grids.z_scattering.size(), nr = config.input.nr, nj = j_values.size(), nrot = basis.size();
  const std::vector<double> vib = doubles(dir.path() / "frames" / "pop_vj_z.bin");
  const std::vector<double> pop = doubles(dir.path() / "frames" / "pop_jmj_z.bin");
  const std::size_t rows = pop.size() / nrot;  // records x nz
  ASSERT_EQ(rows, 4 * nz);
  std::vector<double> from_vib(rows * nj, 0.0), from_pop(rows * nj, 0.0);
  for (std::size_t row = 0; row < rows; ++row) {
    for (std::size_t v = 0; v < nr; ++v)
      for (std::size_t k = 0; k < nj; ++k) from_vib[row * nj + k] += vib[(row * nr + v) * nj + k];
    for (std::size_t s = 0; s < nrot; ++s) {
      const std::size_t k = std::ranges::find(j_values, basis.j(int(s))) - j_values.begin();
      from_pop[row * nj + k] += pop[row * nrot + s];
    }
  }
  double scale = 0.0, top_j = 0.0;
  for (std::size_t i = 0; i < from_pop.size(); ++i) {
    scale = spodvr::testing::nan_aware_max(scale, from_pop[i]);
    if (i % nj == nj - 1) top_j = spodvr::testing::nan_aware_max(top_j, from_pop[i]);
  }
  EXPECT_GT(top_j, 1e-6 * scale);  // J = 2 is populated: the check sees the J axis
  const auto failing = [&](const std::vector<double>& values) {
    int n = 0;
    for (std::size_t i = 0; i < values.size(); ++i) n += !(std::abs(values[i] - from_pop[i]) <= 1e-13 * scale);
    return n;
  };
  EXPECT_EQ(failing(from_vib), 0);
  std::vector<double> wrong = from_vib;
  wrong[3 * nz * nj + nj - 1] += 1e-12 * scale;
  EXPECT_EQ(failing(wrong), 1);
}

// The analytic P2 elements (rotation_moments) against a 32-point Gauss-Legendre quadrature of the normalised
// associated Legendre functions (exact for these degrees), for every J, J' <= 6 and every m: the m dependence of the
// J <-> J + 2 couplings and the heteronuclear J <-> J + 1 zeros, which the fixtures (J = 0, 2, m = 0 populated) do
// not reach. Negative: a slip in the diagonal (1 instead of 3 in front of <cos^2>) fails at every (J, m).
TEST(FrameExtras, P2ElementsEqualTheQuadrature) {
  const int j_max = 6;
  const QuadratureRule rule = gauss_associated_legendre(1.0, -1.0, 32, 0);  // descending roots, negative weights
  int checked = 0;
  for (int m = -j_max; m <= j_max; ++m) {
    std::vector<std::vector<double>> legendre;  // [node][J - |m|]
    for (const double x : rule.roots) legendre.push_back(associated_legendre(j_max, m, std::acos(x)));
    for (int j = std::abs(m); j <= j_max; ++j)
      for (int jp = std::abs(m); jp <= j_max; ++jp) {
        double quadrature = 0.0;
        for (std::size_t t = 0; t < rule.roots.size(); ++t) {
          const double x = rule.roots[t];
          const double p2 = (3.0 * x * x - 1.0) / 2.0;
          quadrature += -rule.weights[t] * p2 * legendre[t][j - std::abs(m)] * legendre[t][jp - std::abs(m)];
        }
        EXPECT_NEAR(p2_element(j, jp, m), quadrature, 1e-13) << "J " << j << " J' " << jp << " m " << m;
        ++checked;
        if (j == jp) {  // the slip: (cos^2 - 1) / 2 instead of (3 cos^2 - 1) / 2
          const double slip = (p2_element(j, j, m) * 2.0 + 1.0) / 3.0 / 2.0 - 0.5;
          EXPECT_GT(std::abs(slip - quadrature), 1e-3) << "J " << j << " m " << m;
        }
      }
  }
  EXPECT_EQ(checked, 231);  // sum over m of (7 - |m|)^2: the (J, J') pairs with both J >= |m|
}

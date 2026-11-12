// Tests: surface symmetry of the potential in a run: the grid symmetry from the input, the PES
// self-check, the reduced potential grid, the gather in the potential operator and the
// propagation with and without symmetry.
#include <gtest/gtest.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "comparison.hpp"
#include "core/error.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/dai_light_cu111_potential.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/memory_estimate.hpp"
#include "propagation/propagation.hpp"

using namespace spodvr;
using spodvr::testing::nan_aware_max;
using Json = nlohmann::json;
using cplx = std::complex<double>;

namespace {

// The input of a data-free PES, selected by its SystemID.
InputParameters input_with_id(const std::string& id) {
  InputParameters p;
  p.system_id = id;
  return p;
}

Json case_json(const std::string& name) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + name + "/input.json");
  return Json::parse(file);
}

struct RunSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  explicit RunSetup(const Json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)) {}
};

// The Dai-Light C6v case (nphi = 12) and the same input without UseSymmetry but with the same
// angular grid: the full potential on the same grid points.
const RunSetup& dai_light_symmetric() {
  static const RunSetup s(case_json("cu111_dailight_short_sym"));
  return s;
}

// A non-symmetric test PES: the Dai-Light potential plus a term linear in x.
class TiltedPotential final : public Potential {
 public:
  [[nodiscard]] double value(double x, double y, double z, double r, double theta, double phi, int iz,
                             int ir) const override {
    return base_.value(x, y, z, r, theta, phi, iz, ir) + 1e-3 * x;
  }
  [[nodiscard]] const SystemParameters& system_parameters() const override { return base_.system_parameters(); }
  [[nodiscard]] std::string name() const override { return "Tilted"; }

 private:
  DaiLightCu111Potential base_;
};

// The Dai-Light potential, but NaN at one configuration (x, y, Z, r, theta, phi); it records every
// configuration it is evaluated at, in call order.
class DaiLightWithNaN final : public Potential {
 public:
  using Configuration = std::array<double, 6>;
  explicit DaiLightWithNaN(std::optional<Configuration> nan_at) : nan_at_(nan_at) {}
  [[nodiscard]] double value(double x, double y, double z, double r, double theta, double phi, int iz,
                             int ir) const override {
    const Configuration q = {x, y, z, r, theta, phi};
    calls_.push_back(q);
    if (q == nan_at_) return std::numeric_limits<double>::quiet_NaN();
    return base_.value(x, y, z, r, theta, phi, iz, ir);
  }
  [[nodiscard]] const SystemParameters& system_parameters() const override { return base_.system_parameters(); }
  [[nodiscard]] std::string name() const override { return "DaiLightWithNaN"; }
  [[nodiscard]] const std::vector<Configuration>& calls() const { return calls_; }

 private:
  DaiLightCu111Potential base_;
  std::optional<Configuration> nan_at_;
  mutable std::vector<Configuration> calls_;
};

double scaled_difference(double a, double b) { return std::abs(a - b) / std::max({std::abs(a), std::abs(b), 1.0}); }

void fill_pseudo_random(std::span<cplx> v) {
  for (std::size_t i = 0; i < v.size(); ++i) {
    const double x = static_cast<double>(i);
    v[i] = cplx(std::sin(0.37 * x + 0.1), std::cos(1.31 * x));
  }
}

}  // namespace

// --- grid symmetry of a run ------------------------------------------------------------------

TEST(RunSymmetry, FromTheInput) {
  Json json = case_json("hardwall_small");  // C2v, 4 x 4, nphi = 10
  const RunSetup plain(json);
  const GridSymmetry none = make_grid_symmetry(plain.config, plain.angular);
  EXPECT_TRUE(none.is_trivial());
  EXPECT_EQ(none.n_representatives(), 16);
  EXPECT_EQ(none.point_group(), "none");

  json["UseSymmetry"] = true;
  const RunSetup c2v(json);
  const GridSymmetry s = make_grid_symmetry(c2v.config, c2v.angular);
  EXPECT_EQ(s.point_group(), "C2v");
  EXPECT_EQ(s.operations().size(), 4U);
  EXPECT_EQ(s.n_representatives(), 9);  // Fortran rectangle A: (4/2 + 1)^2
  // (Frozen x or y is rejected by GridSymmetry, see test_symmetry; the input allows 6D only.)
}

// Without nPhisOptimal the Dai-Light case keeps nphi = 10: only the C2v subgroup of C6v is
// usable (the Fortran's triangle A); the other operations are reported as dropped.
TEST(RunSymmetry, SubgroupWithoutNPhisOptimal) {
  Json json = case_json("cu111_dailight_short_sym");
  json["nPhisOptimal"] = false;
  const RunSetup s(json);
  EXPECT_EQ(s.angular.n_phi(), 10);
  const GridSymmetry symmetry = make_grid_symmetry(s.config, s.angular);
  EXPECT_EQ(symmetry.operations().size(), 4U);
  EXPECT_EQ(symmetry.dropped_operations().size(), 8U);
  EXPECT_EQ(symmetry.n_representatives(), 21);
  EXPECT_EQ(dai_light_symmetric().angular.n_phi(), 12);
  EXPECT_EQ(make_grid_symmetry(dai_light_symmetric().config, dai_light_symmetric().angular).n_representatives(), 10);
}

// --- the PES self-check ----------------------------------------------------------------------

TEST(SymmetryCheck, DaiLightPasses) {
  const RunSetup& s = dai_light_symmetric();
  const DaiLightCu111Potential pes;
  const GridSymmetry symmetry = make_grid_symmetry(s.config, s.angular);
  const SymmetryCheck check = check_potential_symmetry(pes, symmetry, s.config, s.grids, s.angular);
  EXPECT_EQ(check.n_samples, 1000);
  EXPECT_EQ(check.n_comparisons, 11000);
  EXPECT_LT(check.max_deviation, 1e-8);  // rounding only, at the steep wall
  EXPECT_TRUE(check.passed());
  EXPECT_NO_THROW(require_symmetric_potential(check, s.config, symmetry));
}

TEST(SymmetryCheck, NonSymmetricPotentialFails) {
  const RunSetup& s = dai_light_symmetric();
  const GridSymmetry symmetry = make_grid_symmetry(s.config, s.angular);
  const SymmetryCheck check = check_potential_symmetry(TiltedPotential(), symmetry, s.config, s.grids, s.angular);
  EXPECT_GT(check.max_deviation, 1e-4);
  EXPECT_FALSE(check.passed());
  try {
    require_symmetric_potential(check, s.config, symmetry);
    FAIL() << "no InputError";
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), "SymmetryTolerance");
    EXPECT_NE(std::string(e.what()).find("C6v"), std::string::npos) << e.what();
  }
}

// One NaN of the PES, far from the first comparison, fails the check at that configuration: a
// maximum taken with `deviation > max` alone skipped it and passed. The error names the value.
TEST(SymmetryCheck, NonFinitePotentialFailsAtItsConfiguration) {
  const RunSetup& s = dai_light_symmetric();
  const GridSymmetry symmetry = make_grid_symmetry(s.config, s.angular);
  // The check evaluates V(q) and then V(g q) for the operations g != identity: the configuration
  // q of sample 500 is call 500 * (number of operations).
  const DaiLightWithNaN finite(std::nullopt);
  ASSERT_TRUE(check_potential_symmetry(finite, symmetry, s.config, s.grids, s.angular).passed());
  const DaiLightWithNaN::Configuration q = finite.calls().at(500 * symmetry.operations().size());

  const SymmetryCheck check = check_potential_symmetry(DaiLightWithNaN(q), symmetry, s.config, s.grids, s.angular);
  EXPECT_FALSE(check.passed());
  EXPECT_TRUE(std::isnan(check.max_deviation));
  EXPECT_TRUE(std::isnan(check.value));
  EXPECT_TRUE(std::isfinite(check.image_value));
  EXPECT_EQ((DaiLightWithNaN::Configuration{check.x, check.y, check.z, check.r, check.theta, check.phi}), q);
  try {
    require_symmetric_potential(check, s.config, symmetry);
    FAIL() << "no InputError";
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), "SystemID");
    EXPECT_NE(std::string(e.what()).find("is not finite"), std::string::npos) << e.what();
    EXPECT_NE(std::string(e.what()).find("V(q) = nan eV"), std::string::npos) << e.what();
  }
}

// The Persson PES is symmetric about the origin (C6v about the top site) but its period is
// cc/0.529167 instead of the grid's alat = cc/0.529177 (Fortran quirk): the self-check wraps the
// images into the cell, as the potential grid uses them, and so detects it. Deep in the repulsive
// wall the difference is hundreds of eV; with the potential capped at 10 eV still meV.
TEST(SymmetryCheck, PerssonFailsBecauseOfItsPeriod) {
  Json json = case_json("cu111_persson_small");
  json["UseSymmetry"] = true;
  json["nPhisOptimal"] = true;
  const RunSetup s(json);
  const auto pes = make_potential(input_with_id("H2Cu111Persson"));
  const GridSymmetry symmetry = make_grid_symmetry(s.config, s.angular);
  const SymmetryCheck check = check_potential_symmetry(*pes, symmetry, s.config, s.grids, s.angular);
  EXPECT_GT(check.max_deviation, 1.0);
  EXPECT_FALSE(check.passed());
  EXPECT_THROW((void)compute_potential_grid(*pes, s.config, s.grids, s.angular), InputError);

  json["ImposeVMax"] = true;
  json["VMax"] = 10.0;
  const RunSetup capped(json);
  const SymmetryCheck capped_check =
      check_potential_symmetry(*pes, symmetry, capped.config, capped.grids, capped.angular);
  EXPECT_GT(capped_check.max_deviation, 1e-4);
  EXPECT_LT(capped_check.max_deviation, 1e-2);

  json["SymmetryTolerance"] = 1e-2;  // the user accepts a few meV
  const RunSetup tolerant(json);
  EXPECT_NO_THROW((void)compute_potential_grid(*pes, tolerant.config, tolerant.grids, tolerant.angular));
}

// --- the reduced potential grid ----------------------------------------------------------------

// The hard wall does not depend on x, y, phi: the C2v-reduced grid gathers to exactly the full one.
TEST(SymmetricPotentialGrid, HardWallIsBitwiseTheFullGrid) {
  Json json = case_json("hardwall_small");
  const RunSetup plain(json);
  json["UseSymmetry"] = true;
  const RunSetup sym(json);
  const auto pes = make_potential(input_with_id("HardWallPES"));
  const PotentialGrid full = compute_potential_grid(*pes, plain.config, plain.grids, plain.angular);
  const PotentialGrid reduced = compute_potential_grid(*pes, sym.config, sym.grids, sym.angular);
  EXPECT_EQ(reduced.values.size() * 16, full.values.size() * 9);
  EXPECT_EQ(reduced.asymptotic, full.asymptotic);
  for (std::size_t izr = 0; izr < full.index.n2; ++izr)
    for (std::size_t ixy = 0; ixy < full.index.n1; ++ixy)
      for (std::size_t itp = 0; itp < full.index.n0; ++itp) ASSERT_EQ(reduced(itp, ixy, izr), full(itp, ixy, izr));
}

// Dai-Light C6v: the gathered potential equals the directly computed one on the same grid to
// rounding (the representatives' values stand for their symmetry-equivalent points).
TEST(SymmetricPotentialGrid, DaiLightGathersTheFullGrid) {
  const RunSetup& s = dai_light_symmetric();
  const DaiLightCu111Potential pes;
  const PotentialGrid reduced = compute_potential_grid(pes, s.config, s.grids, s.angular);
  const PotentialGrid full = compute_potential_grid(
      pes, s.config, s.grids, s.angular, GridSymmetry(s.config.input.nx, s.config.input.ny, s.angular.n_phi()));
  EXPECT_EQ(reduced.symmetry.n_representatives(), 10);
  EXPECT_EQ(reduced.values.size() * 64, full.values.size() * 10);
  double worst = 0.0;
  for (std::size_t izr = 0; izr < full.index.n2; ++izr)
    for (std::size_t ixy = 0; ixy < full.index.n1; ++ixy)
      for (std::size_t itp = 0; itp < full.index.n0; ++itp)
        worst = nan_aware_max(worst, scaled_difference(reduced(itp, ixy, izr), full(itp, ixy, izr)));
  EXPECT_LT(worst, 1e-12);
  std::cout << "  max scaled difference reduced vs full: " << worst << '\n';
}

TEST(SymmetricPotentialGrid, IndependentOfThreadCount) {
  const RunSetup& s = dai_light_symmetric();
  const DaiLightCu111Potential pes;
  omp_set_num_threads(1);
  const PotentialGrid serial = compute_potential_grid(pes, s.config, s.grids, s.angular);
  omp_set_num_threads(4);
  const PotentialGrid parallel = compute_potential_grid(pes, s.config, s.grids, s.angular);
  EXPECT_TRUE(std::equal(serial.values.begin(), serial.values.end(), parallel.values.begin()));
}

TEST(SymmetricPotentialGrid, MemoryEstimateMatches) {
  const RunSetup& s = dai_light_symmetric();
  const PotentialGrid potential = compute_potential_grid(DaiLightCu111Potential(), s.config, s.grids, s.angular);
  const MemoryEstimate e = estimate_memory(s.config, s.grids, s.angular, 1);
  EXPECT_EQ(e.items[0].values, potential.values.size());
  ASSERT_EQ(e.items[1].name, "Surface symmetry maps");
  EXPECT_EQ(e.items[1].values, 2U * 64 + 12U * 12);
  EXPECT_EQ(e.items[2].values, potential.values.size());  // the exp(-i V dt) table, reduced too
}

// --- the potential operator ---------------------------------------------------------------------

namespace {

// The operators of the Dai-Light C6v case, with the reduced potential and with the same values
// expanded to the full grid (no symmetry).
struct OperatorCase {
  const RunSetup& s = dai_light_symmetric();
  FftwFft fft;
  PotentialGrid reduced = compute_potential_grid(DaiLightCu111Potential(), s.config, s.grids, s.angular);
  GasPhase gas = compute_gas_phase(s.config, s.grids, s.angular, reduced.asymptotic, fft, LapackEigensolver());
  RotationalBasis basis{s.angular, s.config.input.j_in, s.config.input.mj_in};
  AngularTransform transform{s.angular, basis, fft};
  PropagationTables tables = make_propagation_tables(s.config, s.grids, basis, reduced.asymptotic);
  ScatteringShape shape{basis.size(), s.config.input.nx, s.config.input.ny, s.config.input.nr, s.config.input.nz};
  AlignedVector<double> reduced_a;  // -V dt at the representatives
  AlignedVector<double> full_a;     // the same values gathered onto every (x, y) point

  OperatorCase() {
    reduced_a = reduced.values;
    premultiply_with_time_step(reduced_a, s.config.input.time_step);
    const std::size_t n_tp = reduced.index.n0;
    const std::size_t n_xy = std::size_t(shape.nx) * shape.ny;
    const std::size_t n_zr = reduced.index.n2;
    full_a.resize(n_tp * n_xy * n_zr);
    for (std::size_t izr = 0; izr < n_zr; ++izr)
      for (std::size_t ixy = 0; ixy < n_xy; ++ixy)
        for (std::size_t itp = 0; itp < n_tp; ++itp)
          full_a[itp + n_tp * (ixy + n_xy * izr)] = -reduced(itp, ixy, izr) * s.config.input.time_step;
  }

  [[nodiscard]] WavePacket wave() const {
    WavePacket psi{ScatteringWave(shape.n_rot, shape.nx * shape.ny, shape.nr * shape.nz),
                   SpecularWave(shape.nr, s.grids.z_specular.size())};
    fill_pseudo_random(psi.scattering.values);
    fill_pseudo_random(psi.specular.values);
    return psi;
  }
};

const OperatorCase& operator_case() {
  static const OperatorCase c;
  return c;
}

}  // namespace

// The gather through the symmetry maps gives bitwise the same result as the full table holding
// the same values, with precomputed exponentials and on the fly (table limit 0), for any thread
// count.
TEST(SymmetricPotentialOperator, BitwiseEqualToTheExpandedPotential) {
  const OperatorCase& c = operator_case();
  const int nzpot = c.s.grids.nz_potential;
  const int ichn = c.basis.initial_index();
  for (const std::size_t limit : {PotentialOperator::default_table_limit, std::size_t(0)}) {
    const PotentialOperator full(c.shape, nzpot, ichn, c.tables, c.transform, c.full_a, limit);
    const PotentialOperator reduced(c.shape, nzpot, ichn, c.tables, c.transform, c.reduced_a, c.reduced.symmetry,
                                    limit);
    EXPECT_EQ(full.uses_precomputed_exponentials(), limit > 0);
    EXPECT_EQ(reduced.uses_precomputed_exponentials(), limit > 0);
    WavePacket a = c.wave();
    WavePacket b = c.wave();
    omp_set_num_threads(1);
    full.apply(a);
    omp_set_num_threads(4);
    reduced.apply(b);
    EXPECT_EQ(a.scattering.values, b.scattering.values) << "table limit " << limit;
    EXPECT_EQ(a.specular.values, b.specular.values) << "table limit " << limit;
  }
}

TEST(SymmetricPotentialOperator, RejectsMismatchingSizes) {
  const OperatorCase& c = operator_case();
  const int nzpot = c.s.grids.nz_potential;
  const int ichn = c.basis.initial_index();
  EXPECT_THROW(PotentialOperator(c.shape, nzpot, ichn, c.tables, c.transform, c.full_a, c.reduced.symmetry),
               std::invalid_argument);
  EXPECT_THROW(PotentialOperator(c.shape, nzpot, ichn, c.tables, c.transform, c.reduced_a), std::invalid_argument);
  EXPECT_THROW(PotentialOperator(c.shape, nzpot, ichn, c.tables, c.transform, c.reduced_a, GridSymmetry(8, 8, 10)),
               std::invalid_argument);
}

// --- propagation ----------------------------------------------------------------------------------

// Ten analysis steps (50 time steps) of the initial wave packet: with symmetry bitwise equal to
// the expanded potential, and within 1e-13 of the potential computed at every grid point.
TEST(SymmetricPropagation, AgreesWithTheFullPotential) {
  const OperatorCase& c = operator_case();
  const RunSetup& s = c.s;
  const PotentialGrid direct =
      compute_potential_grid(DaiLightCu111Potential(), s.config, s.grids, s.angular,
                             GridSymmetry(s.config.input.nx, s.config.input.ny, s.angular.n_phi()));
  AlignedVector<double> direct_a = direct.values;
  premultiply_with_time_step(direct_a, s.config.input.time_step);

  const int nzpot = s.grids.nz_potential;
  const int ichn = c.basis.initial_index();
  const ScatteringKineticOperator kinetic(c.shape, c.tables, c.fft);
  const SpecularOperator specular(c.shape.nr, s.grids.z_specular.size(), nzpot, ichn, *s.grids.z_specular_optical,
                                  s.grids.r_optical, c.tables, c.fft);
  const InitialWavePacket initial = make_initial_wave_packet(s.config, s.grids, c.gas, c.fft);

  auto propagate = [&](const PotentialOperator& potential) {
    const SplitOperatorPropagator propagator(s.config.derived.steps_per_analysis, kinetic, potential, specular);
    WavePacket psi{ScatteringWave(c.shape.n_rot, c.shape.nx * c.shape.ny, c.shape.nr * c.shape.nz), initial.specular};
    for (int step = 0; step < 10; ++step) propagator.propagate_analysis_step(psi);
    return psi;
  };
  const WavePacket symmetric =
      propagate(PotentialOperator(c.shape, nzpot, ichn, c.tables, c.transform, c.reduced_a, c.reduced.symmetry));
  const WavePacket expanded = propagate(PotentialOperator(c.shape, nzpot, ichn, c.tables, c.transform, c.full_a));
  const WavePacket full = propagate(PotentialOperator(c.shape, nzpot, ichn, c.tables, c.transform, direct_a));

  EXPECT_EQ(symmetric.scattering.values, expanded.scattering.values);
  EXPECT_EQ(symmetric.specular.values, expanded.specular.values);
  double worst = 0.0;
  double largest = 0.0;
  for (std::size_t i = 0; i < full.scattering.values.size(); ++i) {
    worst = nan_aware_max(worst, std::abs(symmetric.scattering.values[i] - full.scattering.values[i]));
    largest = nan_aware_max(largest, std::abs(full.scattering.values[i]));
  }
  for (std::size_t i = 0; i < full.specular.values.size(); ++i) {
    worst = nan_aware_max(worst, std::abs(symmetric.specular.values[i] - full.specular.values[i]));
    largest = nan_aware_max(largest, std::abs(full.specular.values[i]));
  }
  EXPECT_GT(largest, 0.01);
  EXPECT_LT(worst, 1e-13 * largest);
  std::cout << "  symmetric vs full: max |difference| " << worst << " (max |psi| " << largest << ")\n";
}

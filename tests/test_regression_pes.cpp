// Regression gates 1 and 2 (docs/validation.md) for the corrugated PES cases against the Fortran reference
// files produced by tools/fortran_reference/run_setup_reference.sh <case>:
//   fortran_setup_stdout.txt  vasy table (gate 1)
//   pes_reference.txt         6D subsample of the potential and the angular grid (gate 1)
//   debug_asyp.dat            rovibrational energies and eigenvectors (gate 2)
//   fortran_setup_stdout.txt  the analysis limits nana, nstata, inchna (gate 2, with the gas phase)
// Tests are skipped when a reference file is missing, except those of the SCM cases and cu111_dailight_oblique,
// which fail.
// ctest runs this program as two tests, RegressionPes and RegressionPesScmNegative, each in one
// process so that a case's calculation is computed once (calculation()); a skipped case fails
// them (tests/CMakeLists.txt). The SCM cases also have the tests behind their looser gate 1b (the
// section at the end).
#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis_reference.hpp"
#include "app/fine_analysis_run.hpp"
#include "comparison.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "pes_gate.hpp"
#include "potential/pes_cut.hpp"
#include "potential/potential.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/scm/scm_correction.hpp"
#include "potential/scm/scm_crp_potential.hpp"
#include "reference_case_input.hpp"
#include "setup_reference.hpp"

using namespace spodvr;
using namespace spodvr::testing;
using spodvr::test::relative_difference;
namespace fs = std::filesystem;

namespace {

// The Fortran grids are evaluated in quad precision and differ from the C++ grids by up to 1 ulp
// (accepted, see docs/fortran_deviations.md). The potential is compared with an absolute floor of 1 eV
// (scaled_difference, pes_gate.hpp).
constexpr double grid_tolerance = 4e-16;        // relative, ~2 ulp
constexpr double vasy_tolerance = 1e-14;        // relative to max(|V|, 1 eV)
constexpr double potential_tolerance = 1e-13;   // relative to max(|V|, 1 eV), full 6D potential
// The real CRP data set (crp_cu111_srp48_small), relative to max(|V|, 1 eV): see potential_tolerance_of.
constexpr double crp_tolerance = 3e-13;
// The SCM cases, relative to max(|V|, 1 eV): V = V_CRP + ΔV, whose terms reach ±125 eV where they
// cancel to a few tenths of an eV, so the CRP's absolute error weighs more (docs/validation.md).
constexpr double scm_potential_tolerance = 1e-12;

struct Calculation {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  std::unique_ptr<Potential> pes;  // the PES the grid was computed from
  PotentialGrid potential;
  GasPhase gas;
};

// A case is a directory name, optionally with "+UseSymmetry": the same input with UseSymmetry
// (the C++ potential on the symmetry-reduced grid, compared with the full Fortran potential).
std::string directory_of(const std::string& name) { return name.substr(0, name.find('+')); }
bool with_symmetry(const std::string& name) { return name.ends_with("+UseSymmetry"); }
// The SCM cases (SystemID CRPDistrtSurface).
const std::set<std::string> scm_cases = {"scm_cu111_srp48_small", "scm_cu111_srp48_small_rescaling_off"};
// The cases whose reference files must exist: a missing file fails their tests where the other cases'
// tests skip. These cases were committed together with their reference files and tests.
const std::set<std::string> required_reference_cases = {"scm_cu111_srp48_small", "scm_cu111_srp48_small_rescaling_off",
                                                        "cu111_dailight_oblique"};
bool references_required(const std::string& name) { return required_reference_cases.contains(directory_of(name)); }

// The 6D tolerance of a case. The real CRP data set (crp_cu111_srp48_small) holds total energies
// up to 77 Hartree (CutsInterpolation = F): CRP D4 (scalar cos where gfortran uses libmvec in
// the reduction of the cuts) gives up to 1.5e-13 there; with libmvec's cos it is 1.5e-14
// (2026-11-11), the level of the other cases. docs/fortran_deviations.md "CRP library".
double potential_tolerance_of(const std::string& name) {
  if (scm_cases.contains(directory_of(name))) return scm_potential_tolerance;
  return directory_of(name) == "crp_cu111_srp48_small" ? crp_tolerance : potential_tolerance;
}

fs::path case_dir(const std::string& name) { return fs::path(SPODVR_REFERENCE_DIR) / directory_of(name); }

// The calculation of a case, computed once per process and shared by every test of the case.
const Calculation& calculation(const std::string& name) {
  static std::map<std::string, std::unique_ptr<Calculation>> cache;
  auto& entry = cache[name];
  if (!entry) {
    std::ifstream file(case_dir(name) / "input.json");
    InputParameters input = parse_input(nlohmann::json::parse(file));
    resolve_repository_paths(input);
    if (with_symmetry(name)) input.use_symmetry = true;
    std::unique_ptr<Potential> pes = make_potential(input);
    Config config = build_config(input, *pes);
    SpatialGrids grids = make_spatial_grids(config);
    AngularGrid angular = make_angular_grid(config);
    PotentialGrid potential = compute_potential_grid(*pes, config, grids, angular);
    GasPhase gas = compute_gas_phase(config, grids, angular, potential.asymptotic, FftwFft(), LapackEigensolver());
    entry = std::make_unique<Calculation>(Calculation{std::move(config), std::move(grids), std::move(angular),
                                                      std::move(pes), std::move(potential), std::move(gas)});
  }
  return *entry;
}

// The 6D reference of a case (pes_reference.txt), read once.
const test::PesReference& pes_reference(const std::string& name) {
  static std::map<std::string, std::unique_ptr<test::PesReference>> cache;
  auto& entry = cache[name];
  if (!entry)
    entry = std::make_unique<test::PesReference>(test::read_pes_reference_6d(case_dir(name) / "pes_reference.txt"));
  return *entry;
}

class RegressionPes : public ::testing::TestWithParam<std::string> {
 protected:
  [[nodiscard]] fs::path reference(const std::string& file) const { return case_dir(GetParam()) / file; }
};

// A missing reference file of a case skips the test, or fails it if the case's references are required
// (ctest fails on a skip too: see the top of the file).
#define REQUIRE_CASE_REFERENCE(name, file)                                                     \
  if (const fs::path missing_file = case_dir(name) / (file); !fs::exists(missing_file)) {      \
    if (references_required(name)) FAIL() << "reference file " << missing_file << " missing";  \
    GTEST_SKIP() << "reference file " << missing_file                                          \
                 << " missing; run tools/fortran_reference/run_setup_reference.sh " << (name); \
  }
#define REQUIRE_REFERENCE(file) REQUIRE_CASE_REFERENCE(GetParam(), file)

}  // namespace

// Gate 1a: the Fortran vasy table (r, V in eV).
TEST_P(RegressionPes, Gate1AsymptoticPotentialTable) {
  REQUIRE_REFERENCE("fortran_setup_stdout.txt");
  const auto table = test::read_vasy_table(reference("fortran_setup_stdout.txt"));
  ASSERT_FALSE(table.empty()) << "vasy table not found";
  const Calculation& c = calculation(GetParam());
  ASSERT_EQ(table.size(), c.grids.r.size());
  double max_difference = 0.0;
  for (int ir = 0; ir < c.grids.r.size(); ++ir) {
    const auto [r_ref, v_ref] = table[ir];
    EXPECT_LE(relative_difference(c.grids.r[ir], r_ref), grid_tolerance) << "ir = " << ir;
    const double v = units::legacy.hartree_to_ev(c.potential.asymptotic[ir]);
    const double difference = scaled_difference(v, v_ref, 1.0);
    max_difference = nan_aware_max(max_difference, difference);
    EXPECT_LE(difference, vasy_tolerance) << "ir = " << ir << ": " << v << " vs " << v_ref;
  }
  std::cout << "  " << GetParam() << ": vasy max scaled difference " << max_difference << '\n';
}

// Gate 1b: the 6D potential (hartree, before the -dt premultiplication) at the subsampled points
// of the Fortran .pes file, and the angular grid.
TEST_P(RegressionPes, Gate1FullPotentialAndAngularGrid) {
  REQUIRE_REFERENCE("pes_reference.txt");
  const test::PesReference ref = test::read_pes_reference_6d(reference("pes_reference.txt"));
  const Calculation& c = calculation(GetParam());
  ASSERT_EQ(ref.nzpot, c.grids.nz_potential);
  ASSERT_EQ(ref.nr, c.grids.r.size());
  ASSERT_EQ(ref.nx, c.grids.x.size());
  ASSERT_EQ(ref.ny, c.grids.y.size());
  ASSERT_EQ(ref.ntheta, c.angular.n_theta());
  ASSERT_EQ(ref.nphi, c.angular.n_phi());
  for (int t = 0; t < ref.ntheta; ++t) EXPECT_NEAR(c.angular.thetas()[t], ref.thetas[t], 1e-15) << t;
  for (int p = 0; p < ref.nphi; ++p) EXPECT_EQ(c.angular.phis()[p], ref.phis[p]) << p;

  // With surface symmetry the Fortran computed V on its triangle 1 only (C6v, nphi a multiple of
  // 6): the same number of (x, y) points as the C++ representatives (one per orbit).
  if (ref.use_symmetry) {
    EXPECT_EQ(ref.n_xy_computed, c.potential.symmetry.n_representatives());
    EXPECT_FALSE(c.potential.symmetry.is_trivial());
  } else {
    EXPECT_EQ(ref.n_xy_computed, ref.nx * ref.ny);
  }
  ASSERT_GT(ref.points.size(), 10000U);
  expect_pes_gate(GetParam(), ref, potential_at_points(c.potential, ref), potential_limits(c.config),
                  units::legacy.ev_to_hartree(1.0), potential_tolerance_of(GetParam()));
}

// Gate 2: rovibrational energies and eigenvectors vs debug_asyp.dat.
TEST_P(RegressionPes, Gate2RovibrationalStates) {
  REQUIRE_REFERENCE("debug_asyp.dat");
  const auto reference_states = test::read_rovibrational_states(reference("debug_asyp.dat"));
  const Calculation& c = calculation(GetParam());
  const auto& states = c.gas.states;
  ASSERT_EQ(reference_states.size(), (states.n_v() * (states.j_max() + 1)));
  for (const auto& [key, ref] : reference_states) {
    const auto [v, j] = key;
    const double e = states.energy(v, j);
    const auto psi = states.wave_function(v, j);
    ASSERT_EQ(ref.psi.size(), psi.size());
    double max_abs_psi = 0.0;
    for (std::size_t i = 0; i < psi.size(); ++i)
      max_abs_psi = nan_aware_max(max_abs_psi, std::abs(psi[i] - ref.psi[i]));
    if (v < c.gas.n_v_analysis) {
      EXPECT_LE(relative_difference(e, ref.energy), 1e-12) << "v = " << v << ", j = " << j;
      EXPECT_LE(max_abs_psi, 1e-10) << "v = " << v << ", j = " << j;
    }
    // Higher levels: LAPACK vs Jacobi; only a loose check.
    EXPECT_LE(relative_difference(e, ref.energy), 1e-10) << "v = " << v << ", j = " << j;
  }
}

// Gate 2: the analysis limits the Fortran set-up prints (RotDifBasis.f90: nana, nstata, inchna; a potential-only run
// stops before jana and nviba, which gate 4 compares) against the C++ analysis basis built on the gas phase.
TEST_P(RegressionPes, Gate2AnalysisLimits) {
  REQUIRE_REFERENCE("fortran_setup_stdout.txt");
  const Calculation& c = calculation(GetParam());
  const AnalysisSetup setup(c.config, c.angular, c.gas);
  expect_analysis_limits(
      compare_analysis_limits(analysis_limits(setup.basis),
                              read_fortran_analysis_limits(reference("fortran_setup_stdout.txt")), limits_of_a_setup),
      limits_of_a_setup);
}

// cu111_dailight_short_sym: C6v surface symmetry (UseSymmetry, nPhisOptimal); the Fortran .pes
// stores the reduced potential and symmetryIndex, expanded by extract_pes.py, and the C++
// potential is gathered through its own maps: the first Fortran gate of the symmetry maps.
INSTANTIATE_TEST_SUITE_P(Cu111, RegressionPes,
                         ::testing::Values("cu111_dailight_small", "cu111_persson_small", "cu111_dailight_short_sym",
                                           "cu111_dailight_oblique"));

// The CRP library (SystemID CRPlibrary) on the synthetic data sets of its tests (pairwise model,
// tests/crplib/reference/<sym>_pairwise/dataset); r and Z partly outside the data (clamped).
INSTANTIATE_TEST_SUITE_P(Crp, RegressionPes,
                         ::testing::Values("crp_c6v_synthetic_small", "crp_c3v_synthetic_small",
                                           "crp_c4v_synthetic_small"));

// The real data set (H2/Cu(111) SRP48, C3v, potentials/H2Cu111SRP_C3v.dir): the upstream input on
// coarse grids, VMax imposed, r mostly beyond Rmax (clamped).
INSTANTIATE_TEST_SUITE_P(CrpReal, RegressionPes, ::testing::Values("crp_cu111_srp48_small"));

// Speed-up S1 for CRP: the C++ potential on the symmetry-reduced grid (the subgroup usable with
// nphi = 10: C2v for C6v and C4v, the mirror for C3v) against the full Fortran potential (which the
// Fortran computes without symmetry for C3v and C4v).
INSTANTIATE_TEST_SUITE_P(CrpSymmetric, RegressionPes,
                         ::testing::Values("crp_c6v_synthetic_small+UseSymmetry", "crp_c3v_synthetic_small+UseSymmetry",
                                           "crp_c4v_synthetic_small+UseSymmetry", "crp_cu111_srp48_small+UseSymmetry"));

// The SCM distorted-surface PES (SystemID CRPDistrtSurface): the real data set plus the correction
// of the example slab (potentials/H2Cu111SRP_C3v.dir), with the slab rescaled to the CRP lattice
// and without; gate 1b at scm_potential_tolerance.
INSTANTIATE_TEST_SUITE_P(Scm, RegressionPes,
                         ::testing::Values("scm_cu111_srp48_small", "scm_cu111_srp48_small_rescaling_off"));

// --- Gate 1c: the write2DPES cut of the PES against the Fortran's 2dpes.dat ------------------------------
// tools/fortran_reference/run_pes_cut_reference.sh writes pes_cut_<k>.dat for each cut k of the case's
// pes_cuts.json, whose keys it adds to DVR.inp; here the same keys go into input.json and through the
// program's input reader and checks. V is the PES itself in eV (no switching, no VMin/VMax), printed by the
// Fortran as e15.8: a value passes within half a printed unit plus potential_tolerance_of(case) x
// max(|V|, 1 eV) (a value on a rounding boundary can print either way). r and Z pass within half a printed
// unit plus 1e-14 bohr: the Fortran's grid points are rounded from quad precision (its Z = 0 is 4.4e-16 in
// the C++'s plain double, docs/fortran_deviations.md).

namespace {

constexpr double cut_grid_floor = 1e-14;  // bohr

// One Fortran cut: its header's x, y (bohr), theta, phi (rad; '#u= x v= y theta= t phi= p', f20.10) and r, Z
// (bohr) and V (eV) per point, in the file's order (r outer, Z inner).
struct FortranPesCut {
  std::vector<double> header;  // x, y, theta, phi
  std::vector<double> r, z, v_ev;
};

FortranPesCut read_fortran_pes_cut(const fs::path& path) {
  std::ifstream file(path);
  FortranPesCut cut;
  std::string line;
  while (std::getline(file, line)) {
    if (line.starts_with("#u=")) {
      std::istringstream fields(line);
      std::string label;
      for (int k = 0; k < 4; ++k) {
        double value = 0.0;
        fields >> label >> value;
        cut.header.push_back(value);
      }
      continue;
    }
    if (line.empty() || line.starts_with("#")) continue;
    std::istringstream fields(line);
    double r = 0.0, z = 0.0, v = 0.0;
    if (!(fields >> r >> z >> v)) throw std::runtime_error(path.string() + ": not 'r Z V': " + line);
    cut.r.push_back(r);
    cut.z.push_back(z);
    cut.v_ev.push_back(v);
  }
  return cut;
}

// The cuts of a case: the objects of its pes_cuts.json.
nlohmann::json pes_cuts(const std::string& name) {
  std::ifstream file(case_dir(name) / "pes_cuts.json");
  return nlohmann::json::parse(file).at("cuts");
}

// The input of cut k (from 1) of a case: its input.json with write2DPES and the cut's keys, read and checked
// by the program's input code.
Config cut_config(const std::string& name, int k) {
  std::ifstream file(case_dir(name) / "input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["write2DPES"] = true;
  const nlohmann::json cut = pes_cuts(name).at(k - 1);  // a copy: the range below must not refer into a temporary
  for (const auto& [key, value] : cut.items()) json[key] = value;
  InputParameters input = parse_input(json);
  resolve_repository_paths(input);
  return build_config(input, *calculation(name).pes);  // the case's PES (the same data set), built once
}

// How a C++ cut compares with a Fortran one: the largest deviation of V and of r, Z in units of its
// tolerance (NaN if a value is NaN) and the number of points beyond the tolerance.
struct CutComparison {
  double worst_v = 0.0;
  double worst_grid = 0.0;
  int failing = 0;
};

CutComparison compare_cut(const PesCut& cpp, const FortranPesCut& fortran, double tolerance) {
  CutComparison result;
  const int nz = static_cast<int>(cpp.z.size());
  for (std::size_t i = 0; i < fortran.v_ev.size(); ++i) {
    const int ir = static_cast<int>(i) / nz;
    const int iz = static_cast<int>(i) % nz;
    const double v = cpp.v(ir, iz);
    const double v_tolerance =
        printed_half_unit(fortran.v_ev[i]) + tolerance * std::max(std::abs(fortran.v_ev[i]), 1.0);
    const double v_ratio = std::abs(v - fortran.v_ev[i]) / v_tolerance;
    const double r_ratio = std::abs(cpp.r[ir] - fortran.r[i]) / (printed_half_unit(fortran.r[i]) + cut_grid_floor);
    const double z_ratio = std::abs(cpp.z[iz] - fortran.z[i]) / (printed_half_unit(fortran.z[i]) + cut_grid_floor);
    const double grid_ratio = nan_aware_max(r_ratio, z_ratio);
    result.worst_v = nan_aware_max(result.worst_v, v_ratio);
    result.worst_grid = nan_aware_max(result.worst_grid, grid_ratio);
    result.failing += (v_ratio <= 1.0 && grid_ratio <= 1.0) ? 0 : 1;
  }
  return result;
}

// The number of the header's four fields (x, y, theta, phi) that differ from the cut's beyond the header's 10
// printed decimals (4 when the file has no header).
int header_mismatches(const PesCut& cpp, const FortranPesCut& fortran) {
  if (fortran.header.size() != 4) return 4;
  const std::vector<double> where = {cpp.x, cpp.y, cpp.theta, cpp.phi};
  int mismatches = 0;
  for (int i = 0; i < 4; ++i) mismatches += (std::abs(fortran.header[i] - where[i]) <= 0.6e-10) ? 0 : 1;
  return mismatches;
}

// Cut k of a case, computed by the C++ from the case's PES and grids, and read from the Fortran file.
struct CutPair {
  PesCut cpp;
  FortranPesCut fortran;
};

CutPair cut_pair(const std::string& name, int k) {
  const Calculation& c = calculation(name);
  const Config config = cut_config(name, k);
  return {compute_pes_cut(*c.pes, config.input, c.grids),
          read_fortran_pes_cut(case_dir(name) / ("pes_cut_" + std::to_string(k) + ".dat"))};
}

}  // namespace

class RegressionPesCut : public ::testing::TestWithParam<std::string> {};

// Gate 1c: every cut of the case (the Fortran defaults x = y = 0, theta = 90, phi = 0; and an off-symmetry
// point and orientation), every point within its tolerance. The references are required.
TEST_P(RegressionPesCut, Gate1cPotentialCut) {
  const std::string name = GetParam();
  const int n_cuts = static_cast<int>(pes_cuts(name).size());
  ASSERT_GE(n_cuts, 2);
  for (int k = 1; k <= n_cuts; ++k) {
    const fs::path file = case_dir(name) / ("pes_cut_" + std::to_string(k) + ".dat");
    ASSERT_TRUE(fs::exists(file)) << file << " missing; run tools/fortran_reference/run_pes_cut_reference.sh " << name;
    const CutPair pair = cut_pair(name, k);
    ASSERT_EQ(pair.fortran.v_ev.size(), pair.cpp.v_ev.size()) << "cut " << k;
    // The Fortran file is of this cut: its header (10 decimals) gives the cut's position and orientation.
    EXPECT_EQ(header_mismatches(pair.cpp, pair.fortran), 0) << "cut " << k;
    const CutComparison result = compare_cut(pair.cpp, pair.fortran, potential_tolerance_of(name));
    EXPECT_EQ(result.failing, 0) << "cut " << k;
    EXPECT_LE(result.worst_v, 1.0) << "cut " << k;
    EXPECT_LE(result.worst_grid, 1.0) << "cut " << k;
    std::cout << "  " << name << " cut " << k << ": " << pair.cpp.v_ev.size() << " points, worst V " << result.worst_v
              << ", worst r/Z " << result.worst_grid << " of the tolerance\n";
  }
}

INSTANTIATE_TEST_SUITE_P(Cases, RegressionPesCut,
                         ::testing::Values("hardwall_small", "cu111_dailight_small", "cu111_persson_small",
                                           "crp_cu111_srp48_small", "scm_cu111_srp48_small"));

namespace {

// Whether V (hartree) is the VMin or VMax clamp of a case.
bool clamped(const Config& config, double v_hartree) {
  return (config.input.impose_v_min && v_hartree <= config.derived.v_min) ||
         (config.input.impose_v_max && v_hartree >= config.derived.v_max);
}

// A cut (eV) against the potential grid at (itp, ixy) (hartree), at the points of Z <= ZminSwitch (not
// switched) where the cut is not clamped: the number of points compared and of those beyond
// tolerance x max(|V|, 1 eV).
struct GridComparison {
  int compared = 0;
  int failing = 0;
};

GridComparison compare_cut_with_grid(const Calculation& c, const PesCut& cut, std::size_t itp, int ixy,
                                     double tolerance) {
  const double ev2au = units::legacy.ev_to_hartree(1.0);  // as compute_potential_grid converts
  const int nr = c.grids.r.size();
  GridComparison result;
  for (int iz = 0; iz < c.grids.nz_potential; ++iz) {
    if (!fortran_compare::rle(c.grids.z_scattering[iz], c.config.input.z_min_switch)) continue;
    for (int ir = 0; ir < nr; ++ir) {
      const double v_hartree = cut.v(ir, iz) * ev2au;
      if (clamped(c.config, v_hartree)) continue;
      const double grid = c.potential(itp, static_cast<std::size_t>(ixy), static_cast<std::size_t>(ir + nr * iz));
      const double allowed = tolerance * std::max(std::abs(cut.v(ir, iz)), 1.0) * ev2au;
      ++result.compared;
      result.failing += (std::abs(v_hartree - grid) <= allowed) ? 0 : 1;
    }
  }
  return result;
}

}  // namespace

// Negative tests of gate 1c on the Dai-Light cut 2: one Fortran value moved by one printed unit away from the
// C++ value fails at exactly that point; a NaN in the C++ cut fails at exactly that point and makes the worst
// deviation NaN; the angles passed to the PES in degrees (a realistic porting slip) fail at most points; the file
// of cut 1 (x = y = 0, theta 90, phi 0) mismatches cut 2 in all four header fields.
TEST(RegressionPesCutChecks, OnePrintedUnitOrNaNFails) {
  const std::string name = "cu111_dailight_small";
  CutPair pair = cut_pair(name, 2);
  ASSERT_GT(pair.fortran.v_ev.size(), 1000U);
  ASSERT_EQ(pair.fortran.v_ev.size(), pair.cpp.v_ev.size());
  ASSERT_EQ(compare_cut(pair.cpp, pair.fortran, potential_tolerance_of(name)).failing, 0);
  EXPECT_EQ(header_mismatches(pair.cpp, read_fortran_pes_cut(case_dir(name) / "pes_cut_1.dat")), 4);
  const Calculation& c = calculation(name);
  Config degrees = cut_config(name, 2);
  degrees.input.theta_2d_pes = rad_to_deg(degrees.input.theta_2d_pes);  // so that the cut gets degrees as radians
  degrees.input.phi_2d_pes = rad_to_deg(degrees.input.phi_2d_pes);
  const PesCut in_degrees = compute_pes_cut(*c.pes, degrees.input, c.grids);
  const int slipped = compare_cut(in_degrees, pair.fortran, potential_tolerance_of(name)).failing;
  EXPECT_GT(slipped, static_cast<int>(pair.fortran.v_ev.size()) / 2);
  std::cout << "  angles in degrees: " << slipped << " of " << pair.fortran.v_ev.size() << " points fail\n";
  const std::size_t i = pair.fortran.v_ev.size() / 2;
  const double unit = 2.0 * printed_half_unit(pair.fortran.v_ev[i]) / 1.01;
  const int nz = static_cast<int>(pair.cpp.z.size());
  const double cpp_value = pair.cpp.v(static_cast<int>(i) / nz, static_cast<int>(i) % nz);
  FortranPesCut moved = pair.fortran;
  moved.v_ev[i] += (cpp_value >= moved.v_ev[i] ? -1.0 : 1.0) * unit;
  EXPECT_EQ(compare_cut(pair.cpp, moved, potential_tolerance_of(name)).failing, 1);
  pair.cpp.v_ev[i] = std::numeric_limits<double>::quiet_NaN();
  const CutComparison with_nan = compare_cut(pair.cpp, pair.fortran, potential_tolerance_of(name));
  EXPECT_EQ(with_nan.failing, 1);
  EXPECT_TRUE(std::isnan(with_nan.worst_v));
}

// The cut against the potential the propagation uses (C++ only, not part of the Fortran gate): on the
// Dai-Light case, a cut at a grid point (x, y) and a DVR node (theta, phi), given in degrees, equals the
// PotentialGrid there at every point that is neither switched (Z > ZminSwitch) nor clamped, to
// potential_tolerance x max(|V|, 1 eV). Negative test: x and y swapped (the grid is 8 x 8) fails at exactly
// the points where the grid's V at (iy, ix) differs from V at (ix, iy) by more than the tolerance.
TEST(RegressionPesCutChecks, CutAtANodeEqualsThePotentialGrid) {
  const std::string name = "cu111_dailight_small";
  const Calculation& c = calculation(name);
  const int ix = 1, iy = 2, itheta = 1, iphi = 3;
  const int ny = c.grids.y.size();
  ASSERT_EQ(c.grids.x.size(), ny);
  const std::size_t itp = iphi + c.angular.n_phi() * itheta;
  const int ixy = iy + ny * ix;
  const double tolerance = potential_tolerance_of(name);
  Config config = cut_config(name, 1);
  config.input.theta_2d_pes = rad_to_deg(c.angular.thetas()[itheta]);
  config.input.phi_2d_pes = rad_to_deg(c.angular.phis()[iphi]);
  config.input.x_2d_pes = c.grids.x[ix];
  config.input.y_2d_pes = c.grids.y[iy];
  const GridComparison at_node =
      compare_cut_with_grid(c, compute_pes_cut(*c.pes, config.input, c.grids), itp, ixy, tolerance);
  EXPECT_GT(at_node.compared, 1000);
  EXPECT_EQ(at_node.failing, 0);

  // The slip: the cut at (x, y) = (y[iy], x[ix]), i.e. the grid's (iy, ix), compared with the grid at (ix, iy).
  config.input.x_2d_pes = c.grids.y[iy];
  config.input.y_2d_pes = c.grids.x[ix];
  const PesCut swapped_cut = compute_pes_cut(*c.pes, config.input, c.grids);
  const GridComparison swapped = compare_cut_with_grid(c, swapped_cut, itp, ixy, tolerance);
  // The expected failures: the points where the grid itself differs between (iy, ix) and (ix, iy).
  const GridComparison grid_swap = compare_cut_with_grid(c, swapped_cut, itp, ix + ny * iy, tolerance);
  EXPECT_EQ(grid_swap.failing, 0);  // the swapped cut is the grid at (iy, ix)
  int differing = 0;
  const double ev2au = units::legacy.ev_to_hartree(1.0);
  const int nr = c.grids.r.size();
  for (int iz = 0; iz < c.grids.nz_potential; ++iz) {
    if (!fortran_compare::rle(c.grids.z_scattering[iz], c.config.input.z_min_switch)) continue;
    for (int ir = 0; ir < nr; ++ir) {
      const std::size_t izr = static_cast<std::size_t>(ir + nr * iz);
      if (clamped(c.config, swapped_cut.v(ir, iz) * ev2au)) continue;
      const double here = c.potential(itp, static_cast<std::size_t>(ixy), izr);
      const double there = c.potential(itp, static_cast<std::size_t>(ix + ny * iy), izr);
      const double allowed = tolerance * std::max(std::abs(swapped_cut.v(ir, iz)), 1.0) * ev2au;
      differing += (std::abs(there - here) <= allowed) ? 0 : 1;
    }
  }
  EXPECT_GT(differing, 0);
  EXPECT_EQ(swapped.failing, differing);
  std::cout << "  cut at a node vs the potential grid: " << at_node.compared << " points, " << at_node.failing
            << " failing; x and y swapped: " << swapped.failing << " failing\n";
}

// --- The SCM cases: the worst point of gate 1b at 1e-12, and that the gate still fails on small errors -
// The two terms are checked separately where both programs evaluate them at the same points (the
// SCM driver's point lists, tests/scm: ΔV bit for bit, V_CRP within 3e-13 of max(|V_CRP|, 1 eV)).
// On the grid they cannot be told apart: the Fortran's grid points are rounded from quad precision
// (1 ulp from the C++'s, docs/fortran_deviations.md), which moves the steep ΔV of the deep holes
// by up to 1.4e-12 eV (|ΔV| ≈ 330 eV at r = 8.6 bohr). The compared subsample shows no cancellation;
// on the whole .pes gate 1b is 2.98e-13 (docs/validation.md: measured 2026-11-15 on the .pes of the
// 8222b8d reference runs, extract_pes.py at stride 1, C++ Release build).

namespace {

// The case of the negative tests (rescaling on: s ≠ 1, the ideal frame differs from the slab's).
const std::string scm_case = "scm_cu111_srp48_small";

// The potential of an SCM case.
const ScmCrpPotential& scm_potential(const Calculation& c) { return dynamic_cast<const ScmCrpPotential&>(*c.pes); }

// The correction of an SCM potential with the given pair-function parameters, built from its slab,
// rescaling factor and skew angle as the potential builds its own.
ScmCorrection correction_with(const ScmCrpPotential& pes, const ScmParameters& parameters) {
  return ScmCorrection(parameters, pes.slab(), pes.rescaling_factor(), pes.crp().system_parameters().skewing_angle);
}

// ΔV (hartree) at the Cartesian (x, y), Z, r, θ, φ: the correction, or a porting mistake of it.
using DeltaV = std::function<double(double x, double y, double z, double r, double theta, double phi)>;

// An SCM potential with another ΔV, composed as ScmCrpPotential::value: (x, y) of (u, v), the CRP
// at (u/s, v/s), then ΔV/Constants_eV2Ha added (eV). Test-only: the negative tests run the program's
// grid code (compute_potential_grid: hartree, gas-phase switching, VMax) on it.
class ScmPotentialWithDeltaV final : public Potential {
 public:
  ScmPotentialWithDeltaV(const ScmCrpPotential& scm, DeltaV delta_v)
      : scm_(scm),
        correction_(correction_with(scm, scm.scm_parameters())),
        crp_scale_(scm.rescaling_factor().value_or(1.0)),
        delta_v_(std::move(delta_v)) {}

  [[nodiscard]] double value(double u, double v, double z, double r, double theta, double phi, int iz,
                             int ir) const override {
    const CartesianVector xy = correction_.to_cartesian(u, v);
    const double v_crp = scm_.crp().value(u / crp_scale_, v / crp_scale_, z, r, theta, phi, iz, ir);
    return v_crp + units::precise.hartree_to_ev(delta_v_(xy.x, xy.y, z, r, theta, phi));
  }
  [[nodiscard]] const SystemParameters& system_parameters() const override { return scm_.system_parameters(); }
  [[nodiscard]] std::string name() const override { return scm_.name(); }

 private:
  const ScmCrpPotential& scm_;
  ScmCorrection correction_;  // for (x, y) of (u, v) only
  double crp_scale_;          // s with rescaling, else 1
  DeltaV delta_v_;
};

// The correction's own ΔV, and two porting mistakes: P17/P18 (the r clamp of the pair-function
// parameters) taken as Å and converted to bohr, and θ and φ swapped once (the Fortran swaps them
// twice, scm_correction.cpp).
DeltaV correct_delta_v(const ScmCorrection& correction) {
  return [&correction](double x, double y, double z, double r, double theta, double phi) {
    return correction.value_hartree(x, y, z, r, theta, phi);
  };
}
DeltaV swapped_angles_delta_v(const ScmCorrection& correction) {
  return [&correction](double x, double y, double z, double r, double theta, double phi) {
    return correction.value_hartree(x, y, z, r, phi, theta);
  };
}
ScmParameters with_r_clamp_in_angstrom(ScmParameters parameters) {
  for (ScmLayerParameters& layer : parameters.layers) {
    layer.r_clamp_min = units::legacy.angstrom_to_bohr(layer.r_clamp_min);
    layer.r_clamp_max = units::legacy.angstrom_to_bohr(layer.r_clamp_max);
  }
  return parameters;
}

// The two terms of V (eV) at the reference points, and the point's coordinates.
struct ScmTermsAtPoint {
  double u, v, z, r, theta, phi;  // lattice (u, v) and Z, r (bohr); θ, φ (rad)
  double crp;                     // V_CRP(u/s, v/s, …)
  double delta_v;                 // ΔV/Constants_eV2Ha
};

// The two terms of V at each reference point of case c, each computed on its own.
std::vector<ScmTermsAtPoint> scm_terms_at_points(const Calculation& c, const test::PesReference& ref) {
  const ScmCrpPotential& pes = scm_potential(c);
  const ScmCorrection correction = correction_with(pes, pes.scm_parameters());
  const double s = pes.rescaling_factor().value_or(1.0);
  std::vector<ScmTermsAtPoint> points;
  points.reserve(ref.points.size());
  for (const test::PesReference::Point& p : ref.points) {
    ScmTermsAtPoint q{c.grids.x[p.ix],
                      c.grids.y[p.iy],
                      c.grids.z_scattering[p.iz],
                      c.grids.r[p.ir],
                      c.angular.thetas()[p.itheta],
                      c.angular.phis()[p.iphi],
                      0.0,
                      0.0};
    const CartesianVector xy = correction.to_cartesian(q.u, q.v);
    q.crp = pes.crp().value(q.u / s, q.v / s, q.z, q.r, q.theta, q.phi, p.iz, p.ir);
    q.delta_v = units::precise.hartree_to_ev(correction.value_hartree(xy.x, xy.y, q.z, q.r, q.theta, q.phi));
    points.push_back(q);
  }
  return points;
}

// The values (hartree) at the reference points of the grid the program computes from a potential.
std::vector<double> grid_values_at_points(const Potential& pes, const Calculation& c, const test::PesReference& ref) {
  return potential_at_points(compute_potential_grid(pes, c.config, c.grids, c.angular), ref);
}

// Gate 1b of a case on a list of values: the floor 1 eV, the limits and the tolerance of the case.
PesGateResult scm_gate(const std::string& name, const std::vector<double>& values) {
  const Calculation& c = calculation(name);
  return compare_pes_values(pes_reference(name), values, potential_limits(c.config), units::legacy.ev_to_hartree(1.0),
                            potential_tolerance_of(name));
}

// Gate 1b of the SCM cases, one instance per case name.
class RegressionPesScm : public RegressionPes {};

}  // namespace

// The total at its tolerance, with the worst point printed: its coordinates, the two terms, the
// margin (tolerance / largest difference) and, as a diagnostic, the largest difference relative to
// max(|V|, |V_CRP|, |ΔV|, 1 eV), the scale of the terms that cancel.
TEST_P(RegressionPesScm, Gate1WorstPointAndCancellationAwareDifference) {
  REQUIRE_REFERENCE("pes_reference.txt");
  const test::PesReference& ref = pes_reference(GetParam());
  const Calculation& c = calculation(GetParam());
  const std::vector<double> values = potential_at_points(c.potential, ref);
  const PesGateResult result = scm_gate(GetParam(), values);
  const std::vector<ScmTermsAtPoint> points = scm_terms_at_points(c, ref);
  const PotentialLimits limits = potential_limits(c.config);

  // In hartree like the gate, so that it never exceeds the gate's value.
  double cancellation_aware = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    const double expected = with_limits(ref.points[i].v, limits);
    const double scale =
        std::max({std::abs(values[i]), std::abs(expected), units::legacy.ev_to_hartree(std::abs(points[i].crp)),
                  units::legacy.ev_to_hartree(std::abs(points[i].delta_v)), units::legacy.ev_to_hartree(1.0)});
    cancellation_aware = nan_aware_max(cancellation_aware, std::abs(values[i] - expected) / scale);
  }
  const ScmTermsAtPoint& w = points.at(result.worst);
  std::cout << "  " << GetParam() << ": max scaled difference " << result.max_difference << " (margin x"
            << scm_potential_tolerance / result.max_difference << " to " << scm_potential_tolerance << ") at u " << w.u
            << ", v " << w.v << ", Z " << w.z << ", r " << w.r << " bohr, theta " << w.theta << ", phi " << w.phi
            << ": V " << units::legacy.hartree_to_ev(values[result.worst]) << " eV (Fortran "
            << units::legacy.hartree_to_ev(ref.points[result.worst].v) << "), V_CRP " << w.crp << ", dV " << w.delta_v
            << " eV; relative to max(|V|, |V_CRP|, |dV|, 1 eV): " << cancellation_aware << '\n';
  EXPECT_LE(result.max_difference, scm_potential_tolerance);
  EXPECT_TRUE(result.failing.empty());
}

INSTANTIATE_TEST_SUITE_P(Scm, RegressionPesScm,
                         ::testing::Values("scm_cu111_srp48_small", "scm_cu111_srp48_small_rescaling_off"));

// The potential the negative tests use is the program's when given the correction's own ΔV: equal
// bit for bit at every reference point, so their values differ from the program's by the mistake only.
TEST(RegressionPesScmNegative, TestPotentialWithTheOwnDeltaVIsTheProgramsBitwise) {
  REQUIRE_CASE_REFERENCE(scm_case, "pes_reference.txt");
  const Calculation& c = calculation(scm_case);
  const ScmCrpPotential& pes = scm_potential(c);
  const ScmCorrection correction = correction_with(pes, pes.scm_parameters());
  const ScmPotentialWithDeltaV same(pes, correct_delta_v(correction));
  int differ = 0;
  for (const test::PesReference::Point& p : pes_reference(scm_case).points) {
    const double args[] = {c.grids.x[p.ix],
                           c.grids.y[p.iy],
                           c.grids.z_scattering[p.iz],
                           c.grids.r[p.ir],
                           c.angular.thetas()[p.itheta],
                           c.angular.phis()[p.iphi]};
    const double expected = pes.value(args[0], args[1], args[2], args[3], args[4], args[5], p.iz, p.ir);
    const double actual = same.value(args[0], args[1], args[2], args[3], args[4], args[5], p.iz, p.ir);
    if (std::bit_cast<std::uint64_t>(actual) != std::bit_cast<std::uint64_t>(expected)) ++differ;
  }
  EXPECT_EQ(differ, 0);
}

// Porting mistakes in ΔV fail the gate at many points: P17/P18 read as Å (the r clamp [0.725, 1.847]
// becomes [1.370, 3.490] bohr), and θ and φ swapped once.
TEST(RegressionPesScmNegative, Gate1FailsForTheRClampInAngstrom) {
  REQUIRE_CASE_REFERENCE(scm_case, "pes_reference.txt");
  const Calculation& c = calculation(scm_case);
  const ScmCrpPotential& pes = scm_potential(c);
  const ScmCorrection mistaken = correction_with(pes, with_r_clamp_in_angstrom(pes.scm_parameters()));
  const ScmPotentialWithDeltaV potential(pes, correct_delta_v(mistaken));
  const PesGateResult result = scm_gate(scm_case, grid_values_at_points(potential, c, pes_reference(scm_case)));
  std::cout << "  r clamp in Angstrom: " << result.failing.size() << " failing points, max scaled difference "
            << result.max_difference << '\n';
  EXPECT_EQ(result.failing.size(), 35244u);  // of 42620
}

TEST(RegressionPesScmNegative, Gate1FailsForAOneSidedThetaPhiSwap) {
  REQUIRE_CASE_REFERENCE(scm_case, "pes_reference.txt");
  const Calculation& c = calculation(scm_case);
  const ScmCrpPotential& pes = scm_potential(c);
  const ScmCorrection correction = correction_with(pes, pes.scm_parameters());
  const ScmPotentialWithDeltaV potential(pes, swapped_angles_delta_v(correction));
  const PesGateResult result = scm_gate(scm_case, grid_values_at_points(potential, c, pes_reference(scm_case)));
  std::cout << "  theta and phi swapped once: " << result.failing.size() << " failing points, max scaled difference "
            << result.max_difference << '\n';
  EXPECT_EQ(result.failing.size(), 39499u);  // of 42620
}

// The boundary: the program's values with one error of 1.01 times the tolerance at the point of the
// largest |V| fail there and nowhere else; with the worst point moved to 0.99 times the tolerance
// (its error's sign kept) they pass. The errors use the documented 1e-12, which must be the case's.
constexpr double documented_scm_tolerance = 1e-12;  // docs/validation.md

TEST(RegressionPesScmNegative, Gate1FailsJustAboveTheToleranceAtTheLargestV) {
  REQUIRE_CASE_REFERENCE(scm_case, "pes_reference.txt");
  const Calculation& c = calculation(scm_case);
  const test::PesReference& ref = pes_reference(scm_case);
  const PotentialLimits limits = potential_limits(c.config);
  const double floor = units::legacy.ev_to_hartree(1.0);
  ASSERT_EQ(potential_tolerance_of(scm_case), documented_scm_tolerance);
  std::size_t largest = 0;
  for (std::size_t i = 0; i < ref.points.size(); ++i)
    if (std::abs(with_limits(ref.points[i].v, limits)) > std::abs(with_limits(ref.points[largest].v, limits)))
      largest = i;
  const double expected = with_limits(ref.points[largest].v, limits);
  std::vector<double> values = potential_at_points(c.potential, ref);
  values[largest] =
      expected + std::copysign(1.01 * documented_scm_tolerance * std::max(std::abs(expected), floor), expected);
  const PesGateResult result = scm_gate(scm_case, values);
  std::cout << "  error 1.01 x tolerance at V = " << units::legacy.hartree_to_ev(expected) << " eV: scaled difference "
            << result.max_difference << '\n';
  ASSERT_EQ(result.failing.size(), 1u);
  EXPECT_EQ(result.failing[0], largest);
}

TEST(RegressionPesScmNegative, Gate1PassesJustBelowTheToleranceAtTheWorstPoint) {
  REQUIRE_CASE_REFERENCE(scm_case, "pes_reference.txt");
  const Calculation& c = calculation(scm_case);
  const test::PesReference& ref = pes_reference(scm_case);
  const double floor = units::legacy.ev_to_hartree(1.0);
  ASSERT_EQ(potential_tolerance_of(scm_case), documented_scm_tolerance);
  std::vector<double> values = potential_at_points(c.potential, ref);
  const std::size_t worst = scm_gate(scm_case, values).worst;
  const double expected = with_limits(ref.points[worst].v, potential_limits(c.config));
  values[worst] = expected + std::copysign(0.99 * documented_scm_tolerance * std::max(std::abs(expected), floor),
                                           values[worst] - expected);
  const PesGateResult result = scm_gate(scm_case, values);
  std::cout << "  worst point at 0.99 x tolerance: max scaled difference " << result.max_difference << '\n';
  EXPECT_TRUE(result.failing.empty());
  EXPECT_EQ(result.worst, worst);
  EXPECT_GT(result.max_difference, 0.98 * documented_scm_tolerance);
}

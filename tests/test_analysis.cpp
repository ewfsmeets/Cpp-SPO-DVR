// Tests: analysis basis (diffraction channels, analysed rotational and rotation-diffraction
// states, channel energies).
#include <gtest/gtest.h>

#include <omp.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <numbers>
#include <set>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/analysis_basis.hpp"
#include "analysis/balint_kurti.hpp"
#include "analysis/s_matrix.hpp"
#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "propagation/analysis_surface.hpp"
#include "propagation/initial_wave_packet.hpp"
#include "propagation/propagation_tables.hpp"

using namespace spodvr;

namespace {

nlohmann::json reference_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return nlohmann::json::parse(file);
}

// Setup up to the analysis basis.
struct AnalysisSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  PotentialGrid potential;
  FftwFft fft;
  GasPhase gas;
  RotationalBasis basis;
  PropagationTables tables;
  AnalysisBasis analysis;

  explicit AnalysisSetup(const nlohmann::json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(config.input.j_max, config.derived.mj_max_basis, config.input.homonuclear, config.input.j_in,
                config.input.mj_in),
        potential(compute_potential_grid(*make_potential(config.input.system_id), config, grids, angular)),
        gas(compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver())),
        basis(angular, config.input.j_in, config.input.mj_in),
        tables(make_propagation_tables(config, grids, basis, potential.asymptotic)),
        analysis(config, angular, basis, gas, tables.xy_kinetic_energies) {}
};

const AnalysisSetup& reference_setup() {
  static const AnalysisSetup setup(reference_json());
  return setup;
}

}  // namespace

TEST(DiffractionOrder, PointGroups) {
  EXPECT_EQ(diffraction_order("C2v", 2, -1), 3);
  EXPECT_EQ(diffraction_order(" c4v ", -2, -2), 4);
  // Hexagonal lattices: n m > 0 counts max(|n|, |m|), otherwise |n| + |m|.
  EXPECT_EQ(diffraction_order("C6v", 2, 1), 2);
  EXPECT_EQ(diffraction_order("C3v", -1, -3), 3);
  EXPECT_EQ(diffraction_order("C3", 2, -1), 3);
  EXPECT_EQ(diffraction_order("C6v", 0, -2), 2);
  EXPECT_THROW((void)diffraction_order("D4h", 0, 0), NotSupportedError);
}

// Fortran stdout of the reference case: "INFO FROM SETTING UP THE DIFFRACTION / ROTATION /
// ROTATION-DIFFRACTION BASIS" and "INFO FROM SETTING UP THE ANALYSIS".
TEST(AnalysisBasis, ReferenceCaseMatchesFortran) {
  const AnalysisBasis& a = reference_setup().analysis;
  EXPECT_EQ(a.diffraction_channels().size(), 16u);  // ndftow
  EXPECT_EQ(a.n_diffraction(), 16);                  // ndifa
  EXPECT_EQ(a.highest_order(), 4);                   // maxord
  EXPECT_EQ(a.max_channels_per_order(), 6);          // maxnum
  EXPECT_EQ(a.max_order(), 4);                       // mdifa
  EXPECT_TRUE(a.diffraction_basis_too_small());      // "Diffractional basis too small" warning
  EXPECT_EQ(a.n_rotation(), 15);                     // nana
  EXPECT_EQ(a.n_states(), 240);                      // nstata
  EXPECT_EQ(a.initial_state() + 1, 97);              // inchna
  EXPECT_EQ(a.j_max(), 4);                           // jana
  EXPECT_EQ(a.n_v(), 2);                             // nviba
}

TEST(AnalysisBasis, DiffractionChannelsAreSortedStablyByOrder) {
  const auto& channels = reference_setup().analysis.diffraction_channels();
  EXPECT_EQ(channels.front().n, 0);
  EXPECT_EQ(channels.front().m, 0);
  for (std::size_t i = 1; i < channels.size(); ++i) {
    EXPECT_LE(channels[i - 1].order, channels[i].order);
    // Within one order the Fortran generation order (n outer, m inner) is kept.
    if (channels[i - 1].order == channels[i].order) {
      EXPECT_LT(std::make_pair(channels[i - 1].n, channels[i - 1].m), std::make_pair(channels[i].n, channels[i].m));
    }
  }
  std::set<std::pair<int, int>> slots;
  for (const DiffractionChannel& c : channels) {
    EXPECT_EQ(c.order, std::abs(c.n) + std::abs(c.m));
    EXPECT_EQ(c.ix, c.n >= 0 ? c.n : 4 + c.n);
    EXPECT_EQ(c.iy, c.m >= 0 ? c.m : 4 + c.m);
    slots.insert({c.ix, c.iy});
  }
  EXPECT_EQ(slots.size(), 16u);  // every FFT slot exactly once
}

TEST(AnalysisBasis, StatesCombineRotationOuterAndDiffractionInner) {
  const AnalysisSetup& s = reference_setup();
  const AnalysisBasis& a = s.analysis;
  for (int is = 0; is < a.n_states(); ++is) {
    const int irot = a.rotation_index(a.state_rotation(is));
    EXPECT_EQ(a.j(is), s.basis.j(irot));
    EXPECT_EQ(a.mj(is), s.basis.mj(irot));
    EXPECT_EQ(&a.channel(is), &a.diffraction_channels()[is % a.n_diffraction()]);
  }
  const int is0 = a.initial_state();
  EXPECT_EQ(a.j(is0), s.config.input.j_in);
  EXPECT_EQ(a.mj(is0), s.config.input.mj_in);
  EXPECT_EQ(a.channel(is0).n, 0);
  EXPECT_EQ(a.channel(is0).m, 0);
}

TEST(AnalysisBasis, ChannelEnergiesAreDiffractionPlusRovibrationalEnergy) {
  const AnalysisSetup& s = reference_setup();
  const AnalysisBasis& a = s.analysis;
  for (int v = 0; v < a.n_v(); ++v)
    for (int is = 0; is < a.n_states(); ++is) {
      const DiffractionChannel& c = a.channel(is);
      const double e_dif = s.tables.xy_kinetic_energies[c.iy + s.config.input.ny * c.ix];
      EXPECT_EQ(a.channel_energy(is, v), e_dif + s.gas.states.energy(v, a.j(is)));
    }
  EXPECT_EQ(a.diffraction_energy(0), s.tables.parallel_energy);
}

TEST(AnalysisBasis, RegularOrderIsJThenMjAscending) {
  const AnalysisBasis& a = reference_setup().analysis;
  const std::vector<int>& rotations = a.regular_rotation_order();
  ASSERT_EQ(static_cast<int>(rotations.size()), a.n_rotation());
  std::vector<std::pair<int, int>> labels;
  for (const int irotana : rotations)
    labels.push_back({a.j(irotana * a.n_diffraction()), a.mj(irotana * a.n_diffraction())});
  EXPECT_TRUE(std::is_sorted(labels.begin(), labels.end()));
  EXPECT_EQ(labels.front(), std::make_pair(0, 0));
  EXPECT_EQ(labels.back(), std::make_pair(4, 4));

  const std::vector<int> states = a.regular_state_order();
  ASSERT_EQ(static_cast<int>(states.size()), a.n_states());
  EXPECT_EQ(std::set<int>(states.begin(), states.end()).size(), states.size());  // a permutation
  for (std::size_t k = 0; k < states.size(); ++k)
    EXPECT_EQ(a.state_diffraction(states[k]), static_cast<int>(k) % a.n_diffraction());
}

TEST(AnalysisBasis, DifOrdAnaLimitsTheAnalysedChannels) {
  nlohmann::json json = reference_json();
  json["DifOrdAna"] = 1;
  const AnalysisSetup s(json);
  EXPECT_EQ(s.analysis.max_order(), 1);
  EXPECT_EQ(s.analysis.n_diffraction(), 5);  // (0,0), (0,±1), (±1,0)
  EXPECT_EQ(s.analysis.n_states(), 5 * 15);
  EXPECT_EQ(s.analysis.initial_state(), 6 * 5);
}

namespace {

using cplx = std::complex<double>;

// Empty wave packet of the reference case.
WavePacket empty_wave_packet(const AnalysisSetup& s) {
  const InputParameters& p = s.config.input;
  return {ScatteringWave(s.basis.size(), p.nx * p.ny, p.nr * p.nz), SpecularWave(p.nr, s.grids.z_specular.size())};
}

double sum_of_squares(std::span<const double> phi, double dr) {
  double sum = 0.0;
  for (const double x : phi) sum += x * x;
  return sum * dr;
}

}  // namespace

// psi(Zana) = phi_vj(r) x (diffraction plane wave (n, m) with the parallel momentum k0 shifted
// in) in one rotational state: exactly one BK coefficient is non-zero, equal to
// nx ny <phi|phi> sqrt(dx dy sin(skew) / (nx ny)). Checks the FFT sign, the slots and the scale.
TEST(BalintKurti, PlaneWaveGivesOneCoefficient) {
  const AnalysisSetup& s = reference_setup();
  const InputParameters& p = s.config.input;
  const AnalysisBasis& a = s.analysis;
  const BalintKurtiAnalyzer analyzer(s.config, s.grids, a, s.gas, s.tables, s.fft);
  const InitialWavePacket wp = make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);
  const AnalysisSurface surface(s.config, s.grids, wp.specular);
  ASSERT_FALSE(surface.fixed());

  const int v = 1;
  const int is_target = 5 * a.n_diffraction() + 7;  // some rotation, some diffraction channel
  const DiffractionChannel& c = a.channel(is_target);
  const int irot = a.rotation_index(a.state_rotation(is_target));
  const auto phi = s.gas.states.wave_function(v, a.j(is_target));

  WavePacket psi = empty_wave_packet(s);
  for (int ir = 0; ir < p.nr; ++ir) {
    auto block = psi.scattering.block(ir + p.nr * s.grids.iz_analysis);
    for (int ix = 0; ix < p.nx; ++ix)
      for (int iy = 0; iy < p.ny; ++iy) {
        const int ixy = iy + p.ny * ix;
        const double phase = 2.0 * std::numbers::pi * (double(c.n * ix) / p.nx + double(c.m * iy) / p.ny);
        block[irot + s.basis.size() * ixy] = phi[ir] * s.tables.shift_backward[ixy] * std::polar(1.0, phase);
      }
  }

  const BkCoefficients bk = analyzer.coefficients(psi, surface);
  const DerivedParameters& d = s.config.derived;
  const double expected = p.nx * p.ny * sum_of_squares(phi, p.dr) *
                          std::sqrt(d.dx * d.dy * std::sin(d.system.skewing_angle) / (p.nx * p.ny));
  for (int vv = 0; vv < a.n_v(); ++vv)
    for (int is = 0; is < a.n_states(); ++is) {
      if (vv == v && is == is_target) {
        EXPECT_NEAR(std::abs(bk(is, vv) - expected), 0.0, 1e-12 * expected);
      } else {
        EXPECT_NEAR(std::abs(bk(is, vv)), 0.0, 1e-10 * expected) << is << ' ' << vv;
      }
    }
}

// With a fixed specular analysis surface the initial channel's coefficient is replaced by the
// projection of the specular wave there (not added to the scattering grid part).
TEST(BalintKurti, FixedSpecularSurfaceReplacesTheInitialChannel) {
  nlohmann::json json = reference_json();
  json["DynamZSpAna"] = false;
  const AnalysisSetup s(json);
  const InputParameters& p = s.config.input;
  const AnalysisBasis& a = s.analysis;
  const BalintKurtiAnalyzer analyzer(s.config, s.grids, a, s.gas, s.tables, s.fft);
  const InitialWavePacket wp = make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);
  const AnalysisSurface surface(s.config, s.grids, wp.specular);
  ASSERT_TRUE(surface.fixed());

  WavePacket psi = empty_wave_packet(s);
  for (auto& x : psi.scattering.values) x = cplx(0.3, -0.1);  // must not show up in the initial channel
  const auto phi0 = s.gas.states.wave_function(0, p.j_in);
  for (int ir = 0; ir < p.nr; ++ir) psi.specular(ir, surface.index()) = cplx(0.0, 2.0) * phi0[ir];

  const BkCoefficients bk = analyzer.coefficients(psi, surface);
  const cplx expected(0.0, 2.0 * sum_of_squares(phi0, p.dr));
  EXPECT_NEAR(std::abs(bk(a.initial_state(), 0) - expected), 0.0, 1e-13);
  EXPECT_NEAR(std::abs(bk(a.initial_state(), 1)), 0.0, 1e-12);  // phi_1 is orthogonal to phi_0
}

TEST(BalintKurti, ResultDoesNotDependOnTheThreadCount) {
  const AnalysisSetup& s = reference_setup();
  const BalintKurtiAnalyzer analyzer(s.config, s.grids, s.analysis, s.gas, s.tables, s.fft);
  const InitialWavePacket wp = make_initial_wave_packet(s.config, s.grids, s.gas, s.fft);
  const AnalysisSurface surface(s.config, s.grids, wp.specular);
  WavePacket psi = empty_wave_packet(s);
  const int n = static_cast<int>(psi.scattering.values.size());
  for (int i = 0; i < n; ++i) psi.scattering.values[i] = cplx(std::sin(0.37 * i), std::cos(0.11 * i));

  const int threads = omp_get_max_threads();
  omp_set_num_threads(1);
  const BkCoefficients one = analyzer.coefficients(psi, surface);
  omp_set_num_threads(4);
  const BkCoefficients four = analyzer.coefficients(psi, surface);
  omp_set_num_threads(threads);
  EXPECT_EQ(one.values, four.values);
}

// Fortran stdout: "total energy = 0.45837825D+00" for the first of 41 energies, 0.01 eV apart.
TEST(SMatrix, AnalysisEnergiesOfTheReferenceCase) {
  const AnalysisSetup& s = reference_setup();
  const AnalysisEnergies e = make_analysis_energies(s.config, s.gas);
  ASSERT_EQ(e.size(), 41);
  EXPECT_NEAR(e.collision.front(), 0.2, 1e-15);
  EXPECT_NEAR(e.collision.back(), 0.6, 1e-14);
  EXPECT_NEAR(e.total.front(), 0.45837825, 0.5e-8);
  EXPECT_NEAR(e.total[1] - e.total[0], 0.01, 1e-14);
  EXPECT_EQ(e.total_hartree[3], e.total[3] / units::legacy.energy_factor());
}

TEST(SMatrix, Nen95UsesTheSpectralResolution) {
  nlohmann::json json = reference_json();
  json["NrOfEs"] = -95;
  const AnalysisSetup s(json);
  const AnalysisEnergies e = make_analysis_energies(s.config, s.gas);
  const double step = Constants::pi / 14000.0 * units::legacy.energy_factor();
  EXPECT_NEAR(e.collision[1] - e.collision[0], step, 1e-14);
  EXPECT_GT(e.collision.front(), 0.2);
  EXPECT_LE(e.collision.front(), 0.2 + step);
  EXPECT_LE(e.collision.back(), 0.6);
  EXPECT_GT(e.collision.back() + step, 0.6);
}

// A constant BK series c accumulates to the geometric sum c dT sum_k exp(i E k dT); open
// channels then get the analytic prefactor, closed channels are zero.
TEST(SMatrix, ConstantSeriesGivesTheGeometricSumAndThePrefactor) {
  const AnalysisSetup& s = reference_setup();
  const AnalysisBasis& a = s.analysis;
  const InputParameters& p = s.config.input;
  const DerivedParameters& d = s.config.derived;
  SMatrixAccumulator acc(s.config, a, s.gas);
  BkCoefficients bk(a.n_states(), a.n_v());
  for (auto& c : bk.values) c = cplx(0.25, -0.5);
  const double dT = d.steps_per_analysis * p.time_step;
  const int n_steps = 37;
  for (int k = 1; k <= n_steps; ++k) acc.add(bk, k * dT);

  const ScatteringResults r = acc.results();
  ASSERT_EQ(static_cast<int>(r.energies.size()), acc.energies().size());
  const double chaev = units::legacy.energy_factor();
  for (const int ien : {0, 20, 40}) {
    const EnergyResult& e = r.energies[ien];
    const double energy = acc.energies().total_hartree[ien];
    const cplx q = std::exp(cplx(0.0, energy * dT));
    const cplx sum = cplx(0.25, -0.5) * dT * q * (1.0 - std::pow(q, n_steps)) / (1.0 - q);

    const double k_total = std::sqrt(2.0 * d.total_mass * e.collision_energy / chaev);
    const double k0 = k_total;  // no parallel momentum in the reference case
    const double weight = std::exp(-(d.kz - k0) * (d.kz - k0) / (2.0 * d.wave_packet_alpha));
    int n_open = 0;
    for (int v = 0; v < a.n_v(); ++v)
      for (int is = 0; is < a.n_states(); ++is) {
        const int k = is + a.n_states() * v;
        const double kz_sq = 2.0 * d.total_mass * (e.total_energy / chaev - a.channel_energy(is, v));
        if (kz_sq <= 0.0) {
          EXPECT_EQ(e.s_matrix[k], cplx(0.0, 0.0));
          EXPECT_EQ(e.open[k], 0);
          continue;
        }
        n_open = n_open + 1;
        const double kz = std::sqrt(kz_sq);
        const cplx expected = sum * std::pow(d.wave_packet_alpha / (2.0 * Constants::pi), 0.25) / d.total_mass *
                              std::exp(cplx(0.0, -kz * p.z_analysis + (d.kz - k0) * p.z0 - Constants::pi)) *
                              std::sqrt(kz * k0 / weight);
        EXPECT_NEAR(std::abs(e.s_matrix[k] - expected), 0.0, 1e-12 * std::abs(expected)) << ien << ' ' << k;
        EXPECT_NEAR(e.probability[k], std::norm(expected), 1e-11 * std::norm(expected));
        EXPECT_NEAR(e.phase[k], std::arg(expected), 1e-10);
      }
    EXPECT_GT(n_open, 0);
    EXPECT_LT(n_open, a.n_states() * a.n_v());
  }
}

// The summed probabilities are consistent with each other and with the state probabilities.
TEST(SMatrix, SumsAndAveragesAreConsistent) {
  const AnalysisSetup& s = reference_setup();
  const AnalysisBasis& a = s.analysis;
  SMatrixAccumulator acc(s.config, a, s.gas);
  BkCoefficients bk(a.n_states(), a.n_v());
  for (int k = 0; k < static_cast<int>(bk.values.size()); ++k)
    bk.values[k] = cplx(std::sin(1.0 + k), std::cos(2.0 * k));
  acc.add(bk, 10.0);
  acc.add(bk, 20.0);
  const EnergyResult e = acc.results().energies[30];

  double total = 0.0;
  double e_vib = 0.0;
  for (int v = 0; v < a.n_v(); ++v) {
    double p_v = 0.0;
    double p_rot = 0.0;
    double p_ord = 0.0;
    for (int is = 0; is < a.n_states(); ++is) p_v += e.probability[is + a.n_states() * v];
    for (int irotana = 0; irotana < a.n_rotation(); ++irotana) p_rot += e.p_rotation[irotana + a.n_rotation() * v];
    for (int i = 0; i < (a.max_order() + 1) * a.n_rotation(); ++i)
      p_ord += e.p_order[i + (a.max_order() + 1) * a.n_rotation() * v];
    EXPECT_NEAR(e.p_vibration[v], p_v, 1e-12 * p_v);
    EXPECT_NEAR(p_rot, p_v, 1e-12 * p_v);
    EXPECT_NEAR(p_ord, p_v, 1e-12 * p_v);
    total += p_v;
    e_vib += s.gas.vibrational_energies[v] * p_v;
  }
  EXPECT_NEAR(e.p_total, total, 1e-12 * total);
  EXPECT_EQ(e.p_reaction, 1.0 - e.p_total);
  EXPECT_NEAR(e.e_vibration, e_vib * units::legacy.energy_factor() / total, 1e-12);
}

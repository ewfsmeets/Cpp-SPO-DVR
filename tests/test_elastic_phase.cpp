// The elastic S-matrix element of a hard-wall run against the exact 1D reflection amplitude: a check
// of the analysis independent of the Fortran, including the phase of the initial state, which is
// measured at the specular analysis line (docs/fortran_bugs/elastic_phase.md).
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

#include <nlohmann/json.hpp>

#include "core/constants.hpp"
#include "core/unit_conversion.hpp"
#include "potential/potential_factory.hpp"
#include "propagation_reference.hpp"

namespace {

using cplx = std::complex<double>;
using spodvr::Constants;

// hardwall_small reduced to JMax 2 (6 rotational states; JIn = JMax, i.e. JMax 0, is an input error since
// r1726) and 2 x 2 (x, y) points: the hard wall has no coupling, so the Z, r dynamics and the elastic element
// are those of the gate case (7000 steps).
nlohmann::json reduced_hard_wall() {
  nlohmann::json json =
      nlohmann::json::parse(std::ifstream(std::filesystem::path(SPODVR_REFERENCE_DIR) / "hardwall_small/input.json"));
  json["nX"] = 2;
  json["nY"] = 2;
  json["JMax"] = 2;
  json["MJMax"] = 2;
  json["writeProbs"] = false;
  json["writeRDProbs"] = false;
  return json;
}

// R(k) of the stationary state psi ~ e^{-ikZ} + R e^{ikZ} in the wall V_Z(Z) = V(Z, r0) - V(far, r0)
// of the hard-wall PES (separable: V = V(r) + V_Z(Z), V_Z = 0 from Z = 7 on), by Numerov
// integration from Z = 0.3 inside the wall (V_Z ≈ 110 eV), matched to the free waves beyond Z = 7.
cplx reflection_amplitude(const spodvr::Potential& pes, double mass, double k) {
  const double chaev = spodvr::units::legacy.energy_factor();
  const double r0 = 1.4;
  const double v_far = pes.value(0.0, 0.0, 20.0, r0, 0.0, 0.0, 0, 0);
  const double h = 1.0e-4;
  const double z_first = 0.3;
  const int n = static_cast<int>((7.2 - z_first) / h) + 1;
  std::vector<double> z(n), w(n), psi(n, 0.0);
  for (int i = 0; i < n; ++i) {
    z[i] = z_first + i * h;
    const double v_z = (pes.value(0.0, 0.0, z[i], r0, 0.0, 0.0, 0, 0) - v_far) / chaev;
    w[i] = 1.0 - h * h * (2.0 * mass * v_z - k * k) / 12.0;  // psi'' = (2 M V_Z - k^2) psi
  }
  psi[1] = 1e-30;
  for (int i = 1; i + 1 < n; ++i) {
    psi[i + 1] = ((12.0 - 10.0 * w[i]) * psi[i] - w[i - 1] * psi[i - 1]) / w[i + 1];
    if (std::abs(psi[i + 1]) > 1e200)
      for (int j = 0; j <= i + 1; ++j) psi[j] *= 1e-200;
  }
  // psi = A e^{-ikZ} + B e^{ikZ} at two points beyond Z = 7; R = B / A.
  const int i1 = n - 400;
  const int i2 = n - 1;
  const cplx i_unit(0.0, 1.0);
  const cplx in1 = std::exp(-i_unit * k * z[i1]), out1 = std::exp(i_unit * k * z[i1]);
  const cplx in2 = std::exp(-i_unit * k * z[i2]), out2 = std::exp(i_unit * k * z[i2]);
  const cplx det = in1 * out2 - out1 * in2;
  const cplx a = (psi[i1] * out2 - psi[i2] * out1) / det;
  const cplx b = (in1 * psi[i2] - in2 * psi[i1]) / det;
  return b / a;
}

double wrapped(double phase) { return std::remainder(phase, 2.0 * Constants::pi); }  // to [-pi, pi]

}  // namespace

// The program's convention makes S = -R for a pure reflection. Between 0.25 and 0.55 eV (inside the packet's
// 0.2-0.6 eV) the phase of the elastic element agrees with arg(-R) to 1e-2 rad and |S| with 1 to 4e-3; the
// residuals measured at dt = 2 were 6.2e-3 rad and 1.8e-3 (the time-step error, O(dt^2), and the norm not yet
// across the line at tEnd). Negative test: the phase referred to Zana although the initial state is measured at the
// specular line Z_sp = 8.64 (the Fortran's fine analysis, and its on-the-fly analysis up to r1719), i.e.
// S e^{ik(Z_sp - Zana)}, fails at every energy.
TEST(ElasticPhase, HardWallElasticElementIsTheExactReflection) {
  const spodvr::testing::CaseSetup setup = spodvr::testing::case_setup(reduced_hard_wall());
  const spodvr::Config& config = setup.config;
  const std::unique_ptr<spodvr::Potential>& pes = setup.pes;
  const spodvr::testing::PropagationRun run = spodvr::testing::run_cpp_propagation(config, *pes);
  ASSERT_TRUE(run.specular_line.has_value());
  EXPECT_NEAR(*run.specular_line, 8.64, 1e-12);  // the dynamic line of hardwall_small (fortran_stdout.txt)
  const double z_offset = *run.specular_line - config.input.z_analysis;

  const double chaev = spodvr::units::legacy.energy_factor();
  const double mass = config.derived.total_mass;
  const int elastic = run.initial_state + run.n_states * (config.derived.v_in_index - 1);
  int n_compared = 0;
  int n_fortran_failing = 0;
  double worst_phase = 0.0;
  double worst_modulus = 0.0;
  for (const spodvr::EnergyResult& e : run.scattering.energies) {
    if (e.z_energy < 0.25 - 1e-9 || e.z_energy > 0.55 + 1e-9) continue;
    const double k = std::sqrt(2.0 * mass * e.z_energy / chaev);
    const cplx exact = -reflection_amplitude(*pes, mass, k);
    const cplx s = e.s_matrix[elastic];
    const double phase_error = std::abs(wrapped(std::arg(s) - std::arg(exact)));
    worst_phase = spodvr::testing::nan_aware_max(worst_phase, phase_error);
    worst_modulus = spodvr::testing::nan_aware_max(worst_modulus, std::abs(std::abs(s) - std::abs(exact)));
    const cplx fortran = s * std::exp(cplx(0.0, k * z_offset));
    if (!(std::abs(wrapped(std::arg(fortran) - std::arg(exact))) <= 1e-2)) n_fortran_failing = n_fortran_failing + 1;
    n_compared = n_compared + 1;
  }
  std::cout << "  elastic element at " << n_compared << " energies: worst |arg S - arg(-R)| " << worst_phase
            << " rad, worst ||S| - |R|| " << worst_modulus << ", Z_sp " << *run.specular_line << '\n';
  EXPECT_EQ(n_compared, 31);
  EXPECT_LE(worst_phase, 1e-2);
  EXPECT_LE(worst_modulus, 4e-3);
  EXPECT_EQ(n_fortran_failing, 31);  // every energy: k (Z_sp - Zana) is 5.26 to 7.80 rad, at least 0.027 from 2 pi
}

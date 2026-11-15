// Tests: the SCM potential (scm_crp_potential.*) against the Fortran example runs (tests/scm/reference/example,
// 768 points): parameters and composition bitwise (SCM D4), V within the CRP's tolerance, V − V_CRP within its
// rounding bound, each with negative controls; the synthetic slabs refused (SCM D9).
#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "comparison.hpp"
#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "potential/scm/scm_correction.hpp"
#include "potential/scm/scm_crp_potential.hpp"
#include "potential/scm/scm_files.hpp"
#include "scm_driver_output.hpp"
#include "scm_reference_runs.hpp"

using namespace spodvr;
using scm_test::EvaluationPoint;
using scm_test::ReferenceRun;
using scm_test::ScmDriverOutput;
using spodvr::testing::nan_aware_max;

namespace {

// The bit pattern of x, for bitwise comparisons (−0 differs from +0, a NaN equals itself).
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

// The CRP data set every reference run used: potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v (SRP48).
const std::filesystem::path data_set = scm_test::example_directory / "H2Cu111SRP_C3v";

// The tolerance of the CRP term (its few-ulp difference: std::cos where gfortran uses libmvec in the
// reduction of the cuts), relative to max(|V_CRP|, 1 eV): that of the program's potential on this
// data set (gate 1 of crp_cu111_srp48_small). Measured here: 4.4e-14.
constexpr double crp_tolerance = 3e-13;

// The input of a run, as the program reads it from the run's files and settings.
ScmInput input_of(const ReferenceRun& run, const ScmDriverOutput& output) {
  ScmInput input;
  input.parameter_file = run.parameters.string();
  input.position_file = run.slab.string();
  input.configuration = std::stoi(output.text.at("slab_number"));
  input.cutoff = output.value("cutoff");
  input.rescale = output.text.at("rescaling") == "T";
  return input;
}

// One reference point: the C++ values next to the Fortran's (eV unless named otherwise).
struct PointValues {
  double total = 0.0;            // C++ V
  double crp = 0.0;              // C++ CRP(u/s, v/s)
  double total_fortran = 0.0;    // the driver's v_total
  double crp_fortran = 0.0;      // the driver's v_crp
  double delta_v_fortran = 0.0;  // the driver's delta_v (hartree)
  double crp_unscaled = 0.0;     // a porting mistake: the CRP at (u, v) instead of (u/s, v/s)
  double delta_v_at_uv = 0.0;    // a porting mistake: ΔV (hartree) with (u, v) taken as Cartesian (x, y)
};

// An example run: its potential's parameters and every point evaluated.
struct RunValues {
  std::string name;
  ScmDriverOutput output;
  SystemParameters parameters;
  std::optional<double> rescaling_factor;
  std::vector<PointValues> points;
};

// The program's values (and those of the controls) at the points of a reference run, next to the
// values the Fortran driver wrote.
RunValues evaluate(const ReferenceRun& run) {
  RunValues result;
  result.name = run.output.string();
  result.output = scm_test::read_scm_driver_output(scm_test::reference_directory / run.output);
  const ScmInput input = input_of(run, result.output);
  const ScmCrpPotential potential(data_set, input);
  result.parameters = potential.system_parameters();
  result.rescaling_factor = potential.rescaling_factor();
  const double s = potential.rescaling_factor().value_or(1.0);
  // ΔV alone, for the control that skips the u, v → x, y conversion.
  const ScmCorrection correction(read_scm_parameters(input.parameter_file), potential.slab(),
                                 potential.rescaling_factor(), potential.crp().system_parameters().skewing_angle);
  const std::vector<EvaluationPoint> points = scm_test::read_points(scm_test::points_file(run));
  for (std::size_t k = 0; k < points.size(); ++k) {
    const EvaluationPoint& p = points[k];
    const int ir = static_cast<int>(k);  // unused by the CRP, as in the driver
    PointValues v;
    v.total = potential.value(p.u, p.v, p.z, p.r, p.theta, p.phi, 0, ir);
    v.crp = potential.crp().value(p.u / s, p.v / s, p.z, p.r, p.theta, p.phi, 0, ir);
    v.total_fortran = result.output.values.at("v_total").at(k);
    v.crp_fortran = result.output.values.at("v_crp").at(k);
    v.delta_v_fortran = result.output.values.at("delta_v").at(k);
    v.crp_unscaled = potential.crp().value(p.u, p.v, p.z, p.r, p.theta, p.phi, 0, ir);
    v.delta_v_at_uv = correction.value_hartree(p.u, p.v, p.z, p.r, p.theta, p.phi);
    result.points.push_back(v);
  }
  return result;
}

// The two example runs (rescaling on and off), evaluated once.
const std::vector<RunValues>& example_runs() {
  static const std::vector<RunValues> runs = [] {
    std::vector<RunValues> result;
    for (const ReferenceRun& run : scm_test::reference_runs())
      if (run.output.parent_path() == "example") result.push_back(evaluate(run));
    return result;
  }();
  return runs;
}

// The CRP tolerance at a point (eV).
double crp_tolerance_at(const PointValues& v) { return crp_tolerance * std::max(std::abs(v.crp_fortran), 1.0); }

// The rounding bound of V − V_CRP (eV): fl(V_CRP + ΔV) and fl(V − V_CRP) each err by at most ε/2 of their
// result, so with the same ΔV the two codes differ by at most ε/2·(|V| + |V − V_CRP|) summed over both.
double rounding_bound_at(const PointValues& v) {
  const double fortran = std::abs(v.total_fortran) + std::abs(v.total_fortran - v.crp_fortran);
  const double cpp = std::abs(v.total) + std::abs(v.total - v.crp);
  return 0.5 * std::numeric_limits<double>::epsilon() * (fortran + cpp);
}

// A value of a point (eV): the C++ one, a porting mistake, or the reference.
using Evaluation = std::function<double(const PointValues&)>;

// The points of a run where |value − reference| exceeds the tolerance (NaN counts as exceeding),
// and the largest ratio of deviation to tolerance.
struct Failures {
  int count = 0;
  double worst_ratio = 0.0;
};
Failures failures(const RunValues& run, const Evaluation& value, const Evaluation& reference,
                  const Evaluation& tolerance) {
  Failures f;
  for (const PointValues& v : run.points) {
    const double ratio = std::abs(value(v) - reference(v)) / tolerance(v);
    f.worst_ratio = nan_aware_max(f.worst_ratio, ratio);
    if (!(ratio <= 1.0)) ++f.count;
  }
  return f;
}

// V as the potential computes it, and with three porting mistakes.
double total(const PointValues& v) { return v.total; }
double total_with_legacy_factor(const PointValues& v) { return v.crp + units::legacy.hartree_to_ev(v.delta_v_fortran); }
double total_without_u_over_s(const PointValues& v) {
  return v.crp_unscaled + units::precise.hartree_to_ev(v.delta_v_fortran);
}
double total_with_delta_v_at_uv(const PointValues& v) { return v.crp + units::precise.hartree_to_ev(v.delta_v_at_uv); }

}  // namespace

// The system parameters are the Fortran's GetCRPDistSurfParams bit for bit: the CRP's, with both
// lattice constants the slab's a_T (4.8896663475 bohr) when rescaling; s = a_T/a₀ is the Fortran's.
TEST(ScmPotentialReference, SystemParametersAreTheFortrans) {
  ASSERT_EQ(example_runs().size(), 2u);
  for (const RunValues& run : example_runs()) {
    SCOPED_TRACE(run.name);
    const ScmDriverOutput& out = run.output;
    EXPECT_EQ(bits(run.parameters.lattice_constant_x), bits(out.value("lattice_x")));
    EXPECT_EQ(bits(run.parameters.lattice_constant_y), bits(out.value("lattice_y")));
    EXPECT_EQ(bits(run.parameters.z_asymptotic), bits(out.value("z_asymptotic")));
    EXPECT_EQ(bits(run.parameters.skewing_angle), bits(out.value("skew_angle")));
    EXPECT_EQ(run.parameters.point_group, out.text.at("point_group"));
    if (out.text.at("rescaling") == "T") {
      ASSERT_TRUE(run.rescaling_factor.has_value());
      EXPECT_EQ(bits(*run.rescaling_factor), bits(out.value("scale")));
      EXPECT_EQ(bits(run.parameters.lattice_constant_x), bits(out.value("slab_lattice_constant")));
    } else {
      EXPECT_FALSE(run.rescaling_factor.has_value());
      EXPECT_EQ(bits(run.parameters.lattice_constant_x), bits(out.value("ideal_lattice_constant")));
    }
  }
}

// V = CRP(u/s, v/s) + ΔV/Constants_eV2Ha bit for bit, with the Fortran's ΔV (equal to the C++ ΔV,
// test ScmCorrectionReference): the composition is the Fortran's (the driver checked the same
// identity for the Fortran's values).
TEST(ScmPotentialReference, ScmD4ValueIsCrpPlusDeltaVInEvBitwise) {
  int compared = 0, differ = 0;
  for (const RunValues& run : example_runs()) {
    for (const PointValues& v : run.points) {
      ++compared;
      if (bits(v.total) != bits(v.crp + units::precise.hartree_to_ev(v.delta_v_fortran))) ++differ;
    }
  }
  EXPECT_EQ(compared, 768);
  EXPECT_EQ(differ, 0);
}

// V and its CRP term are the Fortran's within the CRP tolerance; the ΔV term adds nothing to the
// deviation. Negative controls: ΔV converted with the legacy factor 27.2114 (a relative error of
// 1.3e-7 of ΔV), and the CRP evaluated at (u, v) instead of (u/s, v/s) when rescaling.
TEST(ScmPotentialReference, TotalAndCrpTermWithinTheCrpTolerance) {
  const Evaluation reference_total = [](const PointValues& v) { return v.total_fortran; };
  const Evaluation reference_crp = [](const PointValues& v) { return v.crp_fortran; };
  const Evaluation crp_term = [](const PointValues& v) { return v.crp; };
  for (const RunValues& run : example_runs()) {
    SCOPED_TRACE(run.name);
    const Failures f_total = failures(run, total, reference_total, crp_tolerance_at);
    const Failures f_crp = failures(run, crp_term, reference_crp, crp_tolerance_at);
    std::cout << run.name << ": V within " << f_total.worst_ratio << " of the tolerance, its CRP term within "
              << f_crp.worst_ratio << " (" << run.points.size() << " points)\n";
    EXPECT_LE(f_total.worst_ratio, 1.0);
    EXPECT_LE(f_crp.worst_ratio, 1.0);
    EXPECT_EQ(f_total.count, 0);

    const Failures legacy = failures(run, total_with_legacy_factor, reference_total, crp_tolerance_at);
    const Failures unscaled = failures(run, total_without_u_over_s, reference_total, crp_tolerance_at);
    std::cout << "  controls: legacy factor fails at " << legacy.count << ", CRP at (u, v) at " << unscaled.count
              << "\n";
    // 357 of the 384 points have ΔV ≠ 0, 3 of them too small for the factor's 1.3e-7 to show; without
    // rescaling s = 1, and the CRP control is the CRP itself.
    const bool rescaling = run.rescaling_factor.has_value();
    EXPECT_EQ(legacy.count, 354);
    EXPECT_EQ(unscaled.count, rescaling ? 293 : 0);
  }
}

// V − V_CRP equals the Fortran's within the rounding bound: x, y, u/s and the factor are the Fortran's.
// Negative controls: the legacy factor, the CRP at (u, v) when rescaling, ΔV with (u, v) as Cartesian.
TEST(ScmPotentialReference, DeltaVTermWithinItsRoundingBound) {
  const Evaluation reference = [](const PointValues& v) { return v.total_fortran - v.crp_fortran; };
  const auto term = [](const Evaluation& value) {
    return Evaluation([value](const PointValues& v) { return value(v) - v.crp; });
  };
  for (const RunValues& run : example_runs()) {
    SCOPED_TRACE(run.name);
    const Failures f = failures(run, term(total), reference, rounding_bound_at);
    std::cout << run.name << ": V - V_CRP within " << f.worst_ratio << " of the rounding bound\n";
    EXPECT_LE(f.worst_ratio, 1.0);
    EXPECT_EQ(f.count, 0);

    const Failures legacy = failures(run, term(total_with_legacy_factor), reference, rounding_bound_at);
    const Failures unscaled = failures(run, term(total_without_u_over_s), reference, rounding_bound_at);
    const Failures cartesian = failures(run, term(total_with_delta_v_at_uv), reference, rounding_bound_at);
    std::cout << "  controls: legacy factor fails at " << legacy.count << ", CRP at (u, v) at " << unscaled.count
              << ", ΔV at (u, v) at " << cartesian.count << "\n";
    const bool rescaling = run.rescaling_factor.has_value();
    EXPECT_EQ(legacy.count, 354);
    EXPECT_EQ(unscaled.count, rescaling ? 297 : 0);
    EXPECT_EQ(cartesian.count, 275);
  }
}

// The synthetic runs' slabs have their ideal atoms on their own lattices (a_T = 4.95 and 4.85 bohr),
// not on the SRP48 lattice (a₀ = 4.916 bohr): the Fortran accepted them, the potential refuses them
// (SCM D9), naming PosFilepath. Their ΔV is checked in ScmCorrectionReference.
TEST(ScmPotentialReference, SyntheticSlabsAreNotOnTheCrpLattice) {
  int refused = 0;
  for (const ReferenceRun& run : scm_test::reference_runs()) {
    if (run.output.parent_path() != "synthetic") continue;
    SCOPED_TRACE(run.output.string());
    const ScmDriverOutput output = scm_test::read_scm_driver_output(scm_test::reference_directory / run.output);
    try {
      (void)ScmCrpPotential(data_set, input_of(run, output));
      ADD_FAILURE() << "the slab was accepted";
    } catch (const InputError& e) {
      EXPECT_EQ(e.parameter(), "PosFilepath") << e.what();
      EXPECT_NE(std::string(e.what()).find("CRP lattice"), std::string::npos) << e.what();
      ++refused;
    }
  }
  EXPECT_EQ(refused, 2);
}

// Tests: the SCM correction ΔV (src/potential/scm/scm_correction.*) is bitwise the Fortran's at
// every point of the reference runs (tests/scm/reference/: the driver's delta_v, from
// DistortedSurfaceWrapper.f90 GetDistPESEnergy → SCM_Wrapper.f90 GetSCMPotential).
#include <gtest/gtest.h>

#include <bit>
#include <cstdint>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "potential/scm/scm_correction.hpp"
#include "potential/scm/scm_files.hpp"
#include "scm_driver_output.hpp"
#include "scm_reference_runs.hpp"

using namespace spodvr;
using scm_test::EvaluationPoint;
using scm_test::ReferenceRun;
using scm_test::ScmDriverOutput;

namespace {

// The bit pattern of a double: a NaN differs from every finite reference value.
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

// The correction of a run, built as the program will build it: the files read with the C++
// readers (the run's configuration and cut-off), s = a_T/a₀ from the slab's lattice constant and
// the CRP's (the driver's ideal_lattice_constant) when the run rescales, and the CRP's skew angle.
ScmCorrection correction_of(const ReferenceRun& run, const ScmDriverOutput& output) {
  const ScmParameters parameters = read_scm_parameters(run.parameters);
  const ScmSlab slab = read_scm_slab(run.slab, std::stoi(output.text.at("slab_number")), output.value("cutoff"));
  std::optional<double> rescaling_factor;
  if (output.text.at("rescaling") == "T") {
    rescaling_factor = slab.header.lattice_constant / output.value("ideal_lattice_constant");
    EXPECT_EQ(bits(*rescaling_factor), bits(output.value("scale")));
  }
  return ScmCorrection(parameters, slab, rescaling_factor, output.value("skew_angle"));
}

}  // namespace

// ΔV equals the Fortran's bit for bit at every point of every run. x, y are computed from the
// driver's input (u, v) with ScmCorrection::to_cartesian, the conversion the program will use, and
// are themselves checked against the driver's x, y; ΔV is evaluated at the computed x, y.
TEST(ScmCorrectionReference, DeltaVIsBitwiseTheFortranAtEveryPoint) {
  int runs = 0, total_compared = 0;
  for (const ReferenceRun& run : scm_test::reference_runs()) {
    SCOPED_TRACE(run.output.string());
    const ScmDriverOutput output = scm_test::read_scm_driver_output(scm_test::reference_directory / run.output);
    ASSERT_TRUE(output.done);
    const ScmCorrection correction = correction_of(run, output);
    const std::vector<EvaluationPoint> points = scm_test::read_points(scm_test::points_file(run));
    const std::vector<double>& x = output.values.at("x");
    const std::vector<double>& y = output.values.at("y");
    const std::vector<double>& delta_v = output.values.at("delta_v");
    ASSERT_EQ(x.size(), points.size());
    ASSERT_EQ(y.size(), points.size());
    ASSERT_EQ(delta_v.size(), points.size());

    int compared = 0, xy_differ = 0, delta_v_differ = 0, nonzero = 0;
    for (std::size_t k = 0; k < points.size(); ++k) {
      const EvaluationPoint& p = points[k];
      const CartesianVector xy = correction.to_cartesian(p.u, p.v);
      if (bits(xy.x) != bits(x[k]) || bits(xy.y) != bits(y[k])) ++xy_differ;
      const double value = correction.value_hartree(xy.x, xy.y, p.z, p.r, p.theta, p.phi);
      ++compared;
      if (delta_v[k] != 0.0) ++nonzero;
      if (bits(value) == bits(delta_v[k])) continue;
      if (++delta_v_differ <= 5)
        ADD_FAILURE() << "point " << k + 1 << ": C++ " << value << " (0x" << std::hex << bits(value) << "), Fortran "
                      << delta_v[k] << " (0x" << bits(delta_v[k]) << std::dec << ")";
    }
    std::cout << run.output.string() << ": " << compared << " points compared bitwise, " << delta_v_differ
              << " differ (" << nonzero << " with ΔV ≠ 0)\n";
    RecordProperty(run.output.stem().string() + "_points_compared", compared);
    EXPECT_EQ(xy_differ, 0);
    EXPECT_EQ(delta_v_differ, 0);
    // Most points have a nonzero ΔV, so the comparison is not one of zeros.
    EXPECT_GT(nonzero, compared / 2);
    ++runs;
    total_compared += compared;
  }
  // Every point of every run: the reference data hold 384 + 384 + 284 + 284 points.
  EXPECT_EQ(runs, 4);
  EXPECT_EQ(total_compared, 1336);
}

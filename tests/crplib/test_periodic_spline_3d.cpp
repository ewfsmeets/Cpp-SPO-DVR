// Tests: the 3D spline periodic in x, y and clamped in z (src/crplib/periodic_spline_3d.*), used by
// the C4v atom surface: properties and bitwise equality with the Fortran SPLINE3D/SPLINT3D.
#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "crplib/errors.hpp"
#include "crplib/periodic_spline_3d.hpp"
#include "fortran_reference.hpp"

using crplib::EvaluationError;
using crplib::detail::PeriodicSpline3D;
using crplib_test::DriverInput;
using crplib_test::DriverOutput;
using crplib_test::reference_file;

namespace {

constexpr double period = 4.8;

// Periodic in x, y; cubic in z (reproduced exactly by the clamped spline with the exact slopes).
double model(double x, double y, double z) {
  const double k = 2.0 * M_PI / period;
  return (1.0 + 0.3 * std::cos(k * x) * std::cos(k * y)) * (0.2 + 0.5 * z - 0.1 * z * z + 0.02 * z * z * z);
}
double model_dz(double x, double y, double z) {
  const double k = 2.0 * M_PI / period;
  return (1.0 + 0.3 * std::cos(k * x) * std::cos(k * y)) * (0.5 - 0.2 * z + 0.06 * z * z);
}

PeriodicSpline3D make_model_spline() {
  std::vector<double> x(5), z = {-1.0, -0.3, 0.4, 1.2, 2.0, 3.1, 4.0};
  for (int i = 0; i < 5; ++i) x[i] = i * period * 0.25;
  std::vector<double> f, first, last;
  for (double zz : z)
    for (double y : x)
      for (double xx : x) f.push_back(model(xx, y, zz));
  for (double y : x) {
    for (double xx : x) {
      first.push_back(model_dz(xx, y, z.front()));
      last.push_back(model_dz(xx, y, z.back()));
    }
  }
  return PeriodicSpline3D(x, x, z, f, first, last);
}

}  // namespace

TEST(PeriodicSpline3D, KnotsAndPeriodicity) {
  const PeriodicSpline3D s = make_model_spline();
  for (double z : {-1.0, 0.4, 2.0, 4.0}) {
    EXPECT_NEAR(s.evaluate(0.0, 1.2, z).value, model(0.0, 1.2, z), 1e-13);
    EXPECT_NEAR(s.evaluate(2.4, 2.4, z).value, model(2.4, 2.4, z), 1e-13);
    // x = 0 and x = period are the same point of the periodic function.
    const auto a = s.evaluate(0.0, 0.7, z), b = s.evaluate(period, 0.7, z);
    EXPECT_NEAR(a.value, b.value, 1e-13);
    EXPECT_NEAR(a.d_dx, b.d_dx, 1e-12);
  }
}

TEST(PeriodicSpline3D, ExactInZAndGradientMatchesDifferences) {
  const PeriodicSpline3D s = make_model_spline();
  for (double z = -1.0; z <= 4.0; z += 0.37) EXPECT_NEAR(s.evaluate(1.2, 3.6, z).value, model(1.2, 3.6, z), 1e-12);
  const double h = 1e-5;
  for (const auto& [x, y, z] : {std::tuple{0.7, 1.9, 0.3}, std::tuple{3.3, 0.2, 2.6}, std::tuple{4.1, 4.4, -0.6}}) {
    const auto v = s.evaluate(x, y, z);
    EXPECT_NEAR(v.d_dx, (s.evaluate(x + h, y, z).value - s.evaluate(x - h, y, z).value) / (2 * h), 1e-7);
    EXPECT_NEAR(v.d_dy, (s.evaluate(x, y + h, z).value - s.evaluate(x, y - h, z).value) / (2 * h), 1e-7);
    EXPECT_NEAR(v.d_dz, (s.evaluate(x, y, z + h).value - s.evaluate(x, y, z - h).value) / (2 * h), 1e-7);
  }
}

TEST(PeriodicSpline3D, OutsideIsAnError) {
  const PeriodicSpline3D s = make_model_spline();
  EXPECT_THROW((void)s.evaluate(1.0, 1.0, -1.01), EvaluationError);  // Fortran SPLINT3D: STOP
  EXPECT_THROW((void)s.evaluate(-0.1, 1.0, 0.0), EvaluationError);
}

TEST(FortranSpline3D, BitwiseEqual) {
  const auto output_file = reference_file("kernels/spline3d.out");
  REQUIRE_CRP_REFERENCE(output_file);
  DriverInput in(reference_file("kernels/spline3d.in"));
  const int nx = in.integer(), ny = in.integer(), nz = in.integer(), m = in.integer();
  const std::vector<double> x = in.reals(nx), y = in.reals(ny), z = in.reals(nz);
  const std::vector<double> f = in.reals(nx * ny * nz), d1 = in.reals(nx * ny), d2 = in.reals(nx * ny);
  const PeriodicSpline3D s(x, y, z, f, d1, d2);
  DriverOutput out(output_file);
  const int n = nx * ny * nz;
  const auto& c = s.coefficients();
  auto compare = [&](const std::string& label, const std::vector<double>& values) {
    const auto fortran = out.next_bits(label, n);
    for (int i = 0; i < n; ++i) EXPECT_BITS(values[i], fortran[i]) << " " << label << " i = " << i;
  };
  compare("cx", c.xx);
  compare("cy", c.yy);
  compare("cxy", c.xxyy);
  compare("cz", c.zz);
  compare("czx", c.zzxx);
  compare("czy", c.zzyy);
  compare("czxy", c.zzxxyy);
  for (int q = 0; q < m; ++q) {
    const double qx = in.real(), qy = in.real(), qz = in.real();
    const auto v = s.evaluate(qx, qy, qz);
    const auto fortran = out.next_bits("Q", 4);
    EXPECT_BITS(v.value, fortran[0]);
    EXPECT_BITS(v.d_dx, fortran[1]);
    EXPECT_BITS(v.d_dy, fortran[2]);
    EXPECT_BITS(v.d_dz, fortran[3]);
  }
}

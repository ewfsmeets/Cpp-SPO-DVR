// Tests: the angular interpolation at one site (src/crplib/site_interpolation.*): each site function
// reproduces its cut values (and their Z, r derivatives) at the cut orientations of Appendix A of
// docs/crp_library.md, is symmetric under the exchange of the two atoms, and has angular
// derivatives that match finite differences. No data set needed.
#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <numbers>
#include <string>
#include <vector>

#include "crplib/site_interpolation.hpp"

using namespace crplib::detail;

namespace {

using SiteFunction = std::function<SiteValue(std::span<const BicubicValue>, const ThetaTrig&, double)>;

struct Orientation {
  double theta, phi;  // degrees
};

struct Case {
  std::string name;
  SiteFunction f;
  std::vector<Orientation> cuts;
};

const std::vector<Case>& cases() {
  static const std::vector<Case> all = {
      {"hexagonal_top", hexagonal_top, {{0, 0}, {90, 0}, {90, 30}}},
      {"hexagonal_bridge", hexagonal_bridge, {{0, 0}, {90, 0}, {90, 60}, {90, 90}}},
      {"hexagonal_hollow", hexagonal_hollow, {{0, 0}, {90, 0}, {90, 30}, {45, 30}, {45, 210}}},
      {"c6v_t2h", c6v_t2h, {{0, 0}, {90, 30}, {90, 120}, {45, 30}, {45, 210}, {45, 120}}},
      {"c3v_t2h", c3v_t2h, {{0, 0}, {90, 30}, {90, 120}, {45, 30}, {45, 120}, {45, 210}}},
      {"c3v_t2f", c3v_t2f, {{0, 0}, {90, 240}, {90, 330}, {45, 150}, {45, 240}, {45, 330}}},
      {"square_top", square_top, {{0, 0}, {90, 0}, {90, 45}}},
      {"square_bridge", square_bridge, {{0, 0}, {90, 0}, {90, 45}, {90, 90}}},
      {"square_hollow", square_hollow, {{0, 0}, {90, 0}, {90, 45}}},
      {"square_t2h", square_t2h, {{0, 0}, {45, 225}, {45, 45}, {90, 45}, {90, 135}, {45, 135}}},
  };
  return all;
}

// Arbitrary cut values (value, d/dZ, d/dr) for n cuts.
std::vector<BicubicValue> cut_values(int n) {
  std::vector<BicubicValue> v;
  for (int k = 0; k < n; ++k) v.push_back({0.3 + 0.7 * std::sin(1.3 * k + 0.2), -0.4 + 0.11 * k, 1.1 - 0.23 * k * k});
  return v;
}

constexpr double degree = std::numbers::pi / 180.0;

SiteValue at(const Case& c, std::span<const BicubicValue> v, double theta, double phi) {
  return c.f(v, theta_trig(theta), phi);
}

}  // namespace

class SiteInterpolation : public testing::TestWithParam<int> {
 protected:
  const Case& c() const { return cases()[GetParam()]; }
};

TEST_P(SiteInterpolation, ReproducesTheCutsAtTheirOrientations) {
  const std::vector<BicubicValue> v = cut_values(static_cast<int>(c().cuts.size()));
  for (std::size_t k = 0; k < c().cuts.size(); ++k) {
    const SiteValue s = at(c(), v, c().cuts[k].theta * degree, c().cuts[k].phi * degree);
    EXPECT_NEAR(s.v, v[k].value, 1e-14) << c().name << " cut " << k;
    EXPECT_NEAR(s.d_z, v[k].d_dy, 1e-14) << c().name << " cut " << k;
    EXPECT_NEAR(s.d_r, v[k].d_dx, 1e-14) << c().name << " cut " << k;
  }
}

TEST_P(SiteInterpolation, HomonuclearSymmetry) {
  const std::vector<BicubicValue> v = cut_values(static_cast<int>(c().cuts.size()));
  for (double theta : {0.1, 0.7, 1.3, 2.0}) {
    for (double phi : {-1.0, 0.4, 2.5}) {
      const SiteValue a = at(c(), v, theta, phi);
      const SiteValue b = at(c(), v, std::numbers::pi - theta, phi + std::numbers::pi);
      EXPECT_NEAR(a.v, b.v, 1e-14) << c().name;
      EXPECT_NEAR(a.d_z, b.d_z, 1e-14) << c().name;
      EXPECT_NEAR(a.d_r, b.d_r, 1e-14) << c().name;
    }
  }
}

TEST_P(SiteInterpolation, AngularDerivativesMatchFiniteDifferences) {
  const std::vector<BicubicValue> v = cut_values(static_cast<int>(c().cuts.size()));
  const double h = 1e-6;
  for (double theta : {0.2, 0.9, 1.4, 2.6}) {
    for (double phi : {-0.8, 0.3, 1.9}) {
      const SiteValue s = at(c(), v, theta, phi);
      EXPECT_NEAR(s.d_phi, (at(c(), v, theta, phi + h).v - at(c(), v, theta, phi - h).v) / (2 * h), 1e-8) << c().name;
      EXPECT_NEAR(s.d_theta, (at(c(), v, theta + h, phi).v - at(c(), v, theta - h, phi).v) / (2 * h), 1e-8) << c().name;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Sites, SiteInterpolation, testing::Range(0, 10),
                         [](const testing::TestParamInfo<int>& param) { return cases()[param.param].name; });

// theta_trig gives exactly the Fortran's per-routine COS(2.00 * Theta) etc. (speed-up S5).
TEST(SiteInterpolationTrig, SameValuesAsPerCall) {
  for (double theta : {0.0, 0.3, 1.5707963267948966, 2.9}) {
    const ThetaTrig t = theta_trig(theta);
    EXPECT_EQ(t.cos2, std::cos(2.00 * theta));
    EXPECT_EQ(t.sin2, std::sin(2.00 * theta));
    EXPECT_EQ(t.cos4, std::cos(4.00 * theta));
    EXPECT_EQ(t.sin4, std::sin(4.00 * theta));
    EXPECT_EQ(t.cos2, std::cos(2.0 * theta));  // the C4v routines write 2.0, 4.0: the same doubles
  }
}

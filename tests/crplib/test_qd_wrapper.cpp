// Tests: the quantum-dynamics wrapper (src/crplib/qd_wrapper.*): clamps, the skewed-coordinate
// transform, the energy shift, the silent V = 0 branch in eV, parameters and the range report.
#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
#include <vector>

#include "crp_synthetic.hpp"
#include "crplib/qd_wrapper.hpp"
#include "crplib/wrapper.hpp"

using namespace crplib;

namespace {

constexpr double ev_per_hartree = 27.2113961;  // the CRP library's factor (UnitConversion.f90)

}  // namespace

class QdWrapperTest : public testing::TestWithParam<Symmetry> {
 protected:
  void SetUp() override {
    data_ = crplib_test::pairwise_data_set(GetParam());
    wrapper_ = std::make_unique<QdWrapper>(data_);
  }
  DataSet data_;
  std::unique_ptr<QdWrapper> wrapper_;
};

TEST_P(QdWrapperTest, Parameters) {
  const QdParameters p = wrapper_->parameters();
  EXPECT_EQ(p.lattice_constant, data_.lattice_constant);  // bohr data: no conversion
  EXPECT_EQ(p.z_asymptotic, data_.z_dependence.x.back());
  EXPECT_EQ(p.point_group, symmetry_name(GetParam()));
  EXPECT_EQ(p.skewing_angle, GetParam() == Symmetry::c4v ? std::numbers::pi / 2.0 : std::numbers::pi / 3.0);
}

// r below Rmin, above Rmax and Z below Zmin are evaluated at Rmin + 1e-12, Rmax - 1e-12, Zmin + 1e-12.
TEST_P(QdWrapperTest, ClampsAreBitwise) {
  const MoleculeSurface& s = wrapper_->surface();
  const double u = 1.3, v = 0.7, theta = 0.9, phi = 0.4;
  EXPECT_EQ(wrapper_->potential_ev(u, v, 2.0, 0.1, theta, phi),
            wrapper_->potential_ev(u, v, 2.0, s.r_min() + 1e-12, theta, phi));
  EXPECT_EQ(wrapper_->potential_ev(u, v, 2.0, 9.0, theta, phi),
            wrapper_->potential_ev(u, v, 2.0, s.r_max() - 1e-12, theta, phi));
  EXPECT_EQ(wrapper_->potential_ev(u, v, -3.0, 1.4, theta, phi),
            wrapper_->potential_ev(u, v, s.z_min() + 1e-12, 1.4, theta, phi));
  EXPECT_EQ(wrapper_->clamped_r(1.4), 1.4);
  EXPECT_EQ(wrapper_->clamped_z(12.0), 12.0);  // no upper clamp
}

// V(eV) = 27.2113961 (V_molecule - VH2(r_last)) with the centre at u a1 + v a2 / a.
TEST_P(QdWrapperTest, SkewedCoordinatesAndEnergyShift) {
  const bool square = GetParam() == Symmetry::c4v;
  const double u = 1.9, v = 0.8, z = 3.2, r = 1.2, theta = 1.1, phi = -0.6;
  const double x = square ? u : u + 0.5 * v, y = square ? v : std::sqrt(3.0) / 2.0 * v;
  const auto [a, b] = crplib_test::atom_pair(x, y, z, r, theta, phi);
  const double molecule = wrapper_->surface().evaluate({a.x, a.y, a.z}, {b.x, b.y, b.z})->v;
  const double expected = ev_per_hartree * molecule - ev_per_hartree * wrapper_->surface().h2_asymptote();
  EXPECT_NEAR(wrapper_->potential_ev(u, v, z, r, theta, phi), expected, 1e-13);
}

TEST_P(QdWrapperTest, RangeReport) {
  const std::vector<double> z = {-2.0, 0.3, 1.0, 3.0, 6.0, 9.0}, r = {0.2, 0.6, 1.4, 3.0, 5.0};
  const RangeReport report = wrapper_->range_report(z, r);
  EXPECT_EQ(report.z_points, 6);
  EXPECT_EQ(report.r_points, 5);
  EXPECT_EQ(report.z_clamped, 2);  // -2.0 and 0.3 are below Zmin = 0.5
  EXPECT_EQ(report.r_clamped, 3);  // 0.2, 0.6 (< Rmin + 1e-12) and 5.0
  EXPECT_EQ(report.zr_points, 30);
  EXPECT_EQ(report.zr_uncovered, 0);  // the clamps keep every point inside the cuts
}

INSTANTIATE_TEST_SUITE_P(Symmetries, QdWrapperTest, testing::Values(Symmetry::c3v, Symmetry::c4v, Symmetry::c6v));

// A cut whose Z grid ends below ZsoMax: there the Fortran returns V = 0 silently, which the wrapper
// shifts to -VH2(r_last) eV; the range report counts these points.
TEST(QdWrapperSilentBranch, ShiftedZeroAndCounted) {
  DataSet data = crplib_test::pairwise_data_set(Symmetry::c6v);
  data.h2_potential.f.back() = 0.002;  // a nonzero asymptote, to see the shift
  MoleculeCut& cut = data.molecule_cuts[4];
  const std::size_t nr = cut.r.size();
  cut.z.resize(cut.z.size() - 3);  // the cut now ends below ZsoMax
  cut.v.resize(nr * cut.z.size());
  ASSERT_LT(cut.z.back(), data.zso_max);
  const QdWrapper wrapper(data);
  const double z = 0.5 * (cut.z.back() + data.zso_max);
  EXPECT_EQ(wrapper.potential_ev(0.3, 0.4, z, 1.4, 0.7, 0.2), -0.002 * ev_per_hartree);
  const std::vector<double> zs = {1.0, z}, rs = {1.4};
  EXPECT_EQ(wrapper.range_report(zs, rs).zr_uncovered, 1);
}

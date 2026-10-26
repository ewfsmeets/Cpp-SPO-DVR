// Tests: associated Legendre functions, rotational basis and the FBR <-> DVR angular transform.
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <vector>

#include "core/constants.hpp"
#include "core/fortran_math.hpp"
#include "fft/fftw_fft.hpp"
#include "grids/angular_grid.hpp"
#include "grids/associated_legendre.hpp"
#include "grids/rotational_basis.hpp"
#include "propagation/angular_transform.hpp"

using namespace spodvr;
using cplx = std::complex<double>;

namespace {

// Angular setups: the reference case and variants (heteronuclear, odd J, MjMax < JMax).
struct AngularCase {
  int j_max;
  int mj_max;
  bool homonuclear;
  int j_in;
};
const AngularCase angular_cases[] = {{4, 4, true, 0}, {3, 2, false, 1}, {5, 3, true, 1}, {6, 1, true, 2}};

// <J m| cos(theta) |J+1 m> for normalised functions.
double cos_coupling(int j, int m) {
  return std::sqrt(((j + 1.0) * (j + 1.0) - m * m) / ((2.0 * j + 1.0) * (2.0 * j + 3.0)));
}

// <J m| cos^2(theta) |J' m> from the cos(theta) couplings (exact, no basis truncation).
double cos2_matrix_element(int j, int jp, int m) {
  if (j == jp) {
    const double up = cos_coupling(j, m);
    const double down = j - 1 >= std::abs(m) ? cos_coupling(j - 1, m) : 0.0;
    return up * up + down * down;
  }
  if (std::abs(j - jp) == 2) {
    const int lower = std::min(j, jp);
    return cos_coupling(lower, m) * cos_coupling(lower + 1, m);
  }
  return 0.0;
}

}  // namespace

TEST(FortranMath, IntegerPowerMatchesRepeatedSquaring) {
  EXPECT_EQ(fortran_integer_power(0.3, 0), 1.0);
  EXPECT_EQ(fortran_integer_power(0.3, 1), 0.3);
  EXPECT_EQ(fortran_integer_power(0.3, 2), 0.3 * 0.3);
  const double x2 = 0.3 * 0.3;
  EXPECT_EQ(fortran_integer_power(0.3, 3), 0.3 * x2);
  EXPECT_EQ(fortran_integer_power(0.3, 5), 0.3 * (x2 * x2));
  EXPECT_EQ(fortran_integer_power(2.0, -2), 0.25);
}

TEST(AssociatedLegendre, ClosedFormsAndPhaseConvention) {
  const double theta = 0.7;
  const double x = std::cos(theta);
  const double s = std::sin(theta);
  const auto p0 = associated_legendre(2, 0, theta);
  EXPECT_NEAR(p0[0], std::sqrt(0.5), 1e-15);
  EXPECT_NEAR(p0[1], std::sqrt(1.5) * x, 1e-15);
  EXPECT_NEAR(p0[2], std::sqrt(2.5) * (3.0 * x * x - 1.0) / 2.0, 1e-15);
  const auto p1 = associated_legendre(2, 1, theta);
  EXPECT_NEAR(p1[0], -std::sqrt(0.75) * s, 1e-15);
  EXPECT_NEAR(p1[1], -std::sqrt(3.75) * x * s, 1e-15);
  const auto p2 = associated_legendre(2, 2, theta);
  EXPECT_NEAR(p2[0], std::sqrt(15.0 / 16.0) * s * s, 1e-15);
  // Negative m: same functions without the Condon-Shortley phase.
  const auto pm1 = associated_legendre(2, -1, theta);
  EXPECT_EQ(pm1[0], -p1[0]);
  EXPECT_EQ(pm1[1], -p1[1]);
  EXPECT_EQ(associated_legendre(2, -2, theta)[0], p2[0]);
}

TEST(AssociatedLegendre, OrthonormalOnTheGaussGrid) {
  for (const AngularCase& c : angular_cases) {
    const AngularGrid grid(c.j_max, c.mj_max, c.homonuclear, c.j_in, 0);
    for (int m = -grid.mj_max(); m <= grid.mj_max(); ++m) {
      for (int j = grid.j_min(); j <= grid.j_max(); j += grid.delta_j()) {
        for (int jp = grid.j_min(); jp <= grid.j_max(); jp += grid.delta_j()) {
          if (std::abs(m) > std::min(j, jp)) continue;
          double overlap = 0.0;
          for (int t = 0; t < grid.n_theta(); ++t) {
            const auto p = associated_legendre(grid.j_max(), m, grid.thetas()[t]);
            overlap += grid.theta_weights()[t] * p[j - std::abs(m)] * p[jp - std::abs(m)];
          }
          EXPECT_NEAR(overlap, j == jp ? 1.0 : 0.0, 1e-14) << "J=" << j << " J'=" << jp << " m=" << m;
        }
      }
    }
  }
}

TEST(RotationalBasis, ReferenceCaseMatchesFortran) {
  // Fortran stdout: 15 rotation channels (nrotw), initial rotational state ichn = 7.
  const AngularGrid grid(4, 4, true, 0, 0);
  const RotationalBasis basis(grid, 0, 0);
  EXPECT_EQ(basis.size(), 15);
  EXPECT_EQ(basis.initial_index(), 6);
  EXPECT_EQ(basis.j(0), 4);
  EXPECT_EQ(basis.mj(0), -4);
  EXPECT_EQ(basis.mj_block(-2).first, 2);
  EXPECT_EQ(basis.mj_block(-2).last, 3);
  EXPECT_EQ(basis.mj_block(0).size(), 3);
  EXPECT_EQ(basis.j(basis.size() - 1), 4);
  EXPECT_EQ(basis.mj(basis.size() - 1), 4);
}

TEST(RotationalBasis, MjMaxBelowJMax) {
  const AngularGrid grid(3, 1, false, 1, 0);
  const RotationalBasis basis(grid, 1, -1);
  // mJ = -1: J = 1, 2, 3; mJ = 0: J = 0..3; mJ = 1: J = 1..3.
  EXPECT_EQ(basis.size(), 10);
  EXPECT_EQ(basis.initial_index(), 0);
}

TEST(AngularTransform, DvrValuesOfSingleBasisFunctions) {
  const FftwFft fft;
  for (const AngularCase& c : angular_cases) {
    const AngularGrid grid(c.j_max, c.mj_max, c.homonuclear, c.j_in, 0);
    const RotationalBasis basis(grid, c.j_in, 0);
    const AngularTransform transform(grid, basis, fft);
    for (int irot = 0; irot < basis.size(); ++irot) {
      std::vector<cplx> fbr(basis.size());
      fbr[irot] = 1.0;
      std::vector<cplx> dvr(transform.dvr_block_size());
      transform.fbr_to_dvr(fbr, dvr);
      const int m = basis.mj(irot);
      for (int t = 0; t < grid.n_theta(); ++t) {
        const auto p = associated_legendre(grid.j_max(), m, grid.thetas()[t]);
        const double radial = std::sqrt(grid.theta_weights()[t] / grid.n_phi()) * p[basis.j(irot) - std::abs(m)];
        for (int k = 0; k < grid.n_phi(); ++k) {
          const cplx expected = radial * std::exp(cplx(0.0, m * grid.phis()[k]));
          EXPECT_NEAR(std::abs(dvr[k + transform.leading_dimension() * t] - expected), 0.0, 1e-14)
              << "irot=" << irot << " t=" << t << " k=" << k;
        }
      }
    }
  }
}

TEST(AngularTransform, RoundTripIsIdentity) {
  const FftwFft fft;
  for (const AngularCase& c : angular_cases) {
    const AngularGrid grid(c.j_max, c.mj_max, c.homonuclear, c.j_in, 0);
    const RotationalBasis basis(grid, c.j_in, 0);
    const AngularTransform transform(grid, basis, fft);
    std::vector<cplx> fbr(basis.size());
    for (int i = 0; i < basis.size(); ++i) fbr[i] = cplx(std::cos(0.7 * i), std::sin(1.3 * i + 0.2));
    std::vector<cplx> dvr(transform.dvr_block_size());
    std::vector<cplx> back(basis.size());
    transform.fbr_to_dvr(fbr, dvr);
    transform.dvr_to_fbr(dvr, back);
    for (int i = 0; i < basis.size(); ++i) EXPECT_NEAR(std::abs(back[i] - fbr[i]), 0.0, 1e-14) << i;
  }
}

TEST(AngularTransform, CosSquaredThetaInDvrMatchesAnalyticMatrix) {
  // Multiplying by cos^2(theta) on the DVR grid represents <J m|cos^2|J' m> exactly as long as
  // the Gauss rule integrates the product (degree J + J' + 2 <= 2 JMax + 1).
  const FftwFft fft;
  for (const AngularCase& c : angular_cases) {
    const AngularGrid grid(c.j_max, c.mj_max, c.homonuclear, c.j_in, 0);
    const RotationalBasis basis(grid, c.j_in, 0);
    const AngularTransform transform(grid, basis, fft);
    const int ld = transform.leading_dimension();
    for (int col = 0; col < basis.size(); ++col) {
      std::vector<cplx> fbr(basis.size());
      fbr[col] = 1.0;
      std::vector<cplx> dvr(transform.dvr_block_size());
      transform.fbr_to_dvr(fbr, dvr);
      for (int t = 0; t < grid.n_theta(); ++t) {
        const double x = std::cos(grid.thetas()[t]);
        for (int k = 0; k < grid.n_phi(); ++k) dvr[k + ld * t] *= x * x;
      }
      std::vector<cplx> result(basis.size());
      transform.dvr_to_fbr(dvr, result);
      for (int row = 0; row < basis.size(); ++row) {
        if (basis.j(row) + basis.j(col) + 2 > 2 * grid.j_max() + 1) continue;
        const double expected =
            basis.mj(row) == basis.mj(col) ? cos2_matrix_element(basis.j(row), basis.j(col), basis.mj(col)) : 0.0;
        EXPECT_NEAR(std::abs(result[row] - expected), 0.0, 1e-14) << "row=" << row << " col=" << col;
      }
    }
  }
}

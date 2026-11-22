// Tests: associated Legendre functions, rotational basis and the FBR <-> DVR angular transform.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <span>
#include <vector>

#include "core/aligned_vector.hpp"
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

// The potential step's DVR -> FBR keeps only the basis's part of exp(-i V dt) psi (the norm budget's leak at the
// surface; docs/output_files.md "norm_budget/"): with a phase that is the same on every DVR point (an isotropic V) the
// norm is kept to rounding; with a phase that varies over (theta, phi) the norm falls by exactly the part of D B c
// outside the basis, |D B c - B B^+ D B c|^2 (B, FBR -> DVR, is an isometry and D unitary). The grid has more DVR
// points than basis functions (JMax = MJMax = 4, homonuclear: 3 x 10 points for 15 functions).
TEST(AngularTransform, AnisotropicPhaseLosesTheNormOutsideTheBasis) {
  const FftwFft fft;
  const AngularGrid grid(4, 4, true, 0, 0);
  const RotationalBasis basis(grid, 0, 0);
  const AngularTransform transform(grid, basis, fft);
  ASSERT_EQ(grid.n_theta() * grid.n_phi(), 30);
  ASSERT_EQ(basis.size(), 15);
  const int ld = transform.leading_dimension();
  std::vector<cplx> fbr(basis.size());
  double norm = 0.0;
  for (int i = 0; i < basis.size(); ++i) {
    fbr[i] = cplx(std::cos(0.7 * i), std::sin(1.3 * i + 0.2));
    norm += std::norm(fbr[i]);
  }
  for (const double anisotropy : {0.0, 1.5}) {  // the phase's variation over the DVR points (rad)
    std::vector<cplx> dvr(transform.dvr_block_size());
    transform.fbr_to_dvr(fbr, dvr);
    for (int t = 0; t < grid.n_theta(); ++t)
      for (int k = 0; k < grid.n_phi(); ++k) {
        const double phi = Constants::two_pi * k / grid.n_phi();
        dvr[k + ld * t] *= std::polar(1.0, -(0.3 + anisotropy * std::cos(grid.thetas()[t]) * std::cos(phi)));
      }
    const std::vector<cplx> phased = dvr;  // D B c (dvr_to_fbr overwrites its input)
    std::vector<cplx> back(basis.size());
    transform.dvr_to_fbr(dvr, back);
    std::vector<cplx> inside(transform.dvr_block_size());  // B B^+ D B c
    transform.fbr_to_dvr(back, inside);
    double phased_norm = 0.0;
    double outside = 0.0;
    for (int t = 0; t < grid.n_theta(); ++t)
      for (int k = 0; k < grid.n_phi(); ++k) {
        phased_norm += std::norm(phased[k + ld * t]);
        outside += std::norm(phased[k + ld * t] - inside[k + ld * t]);
      }
    double back_norm = 0.0;
    for (const cplx& v : back) back_norm += std::norm(v);
    EXPECT_NEAR(phased_norm, norm, 1e-13 * norm) << anisotropy;  // B is an isometry, D unitary
    EXPECT_NEAR(norm - back_norm, outside, 1e-13 * norm) << anisotropy;
    if (anisotropy == 0.0)
      EXPECT_LT(outside, 1e-26 * norm);
    else
      EXPECT_GT(outside, 1e-3 * norm);
  }
}

namespace {

// The same bits: in passes of fewer columns the phi FFT runs with a smaller batch count, which oneMKL and
// FFTW both compute bitwise the same (checked in the oneMKL, FFTW and no-MKL builds, 2026-11-22).
void expect_same_transform(std::span<const cplx> a, std::span<const cplx> b) {
  ASSERT_EQ(a.size(), b.size());
  EXPECT_TRUE(std::equal(a.begin(), a.end(), b.begin()));
}

}  // namespace

// The CPU transforms the (x, y) columns of a point in passes (AngularTransform::columns_per_pass).
// Passes of 1 or 3 columns give the bits of one pass over all 6: every theta sum keeps its order of
// terms (5 rows at a time with a remainder for 3 and 7 theta points, three full blocks for 15); only
// the batch count of the phi FFT changes.
TEST(AngularTransform, ColumnPassesGiveTheBitsOfOnePass) {
  const FftwFft fft;
  constexpr int n_columns = 6;
  for (const int j_max : {4, 12, 28}) {  // 3, 7 and 15 theta points
    SCOPED_TRACE(j_max);
    const AngularGrid grid(j_max, j_max, true, 0, 0);
    const RotationalBasis basis(grid, 0, 0);
    const AngularTransform one_pass(grid, basis, fft, n_columns);
    const std::size_t n_rot = basis.size();
    const std::size_t block = one_pass.dvr_block_size();
    // 64-byte aligned like the program's work arrays: the FFT library picks its algorithm by address.
    ComplexVector fbr(n_rot * n_columns, cplx(0.0, 0.0));
    for (std::size_t i = 0; i < fbr.size(); ++i)
      fbr[i] = cplx(std::cos(0.7 * double(i)), std::sin(1.3 * double(i) + 0.2));
    ComplexVector dvr_one(block * n_columns, cplx(0.0, 0.0));
    one_pass.fbr_to_dvr(fbr, dvr_one);
    ComplexVector work = dvr_one;
    ComplexVector back_one(fbr.size(), cplx(0.0, 0.0));
    one_pass.dvr_to_fbr(work, back_one);

    for (const int pass : {1, 3}) {
      SCOPED_TRACE(pass);
      const AngularTransform passes(grid, basis, fft, pass);
      ComplexVector dvr(dvr_one.size(), cplx(0.0, 0.0));
      ComplexVector back(fbr.size(), cplx(0.0, 0.0));
      for (int c0 = 0; c0 < n_columns; c0 += pass) {
        const std::span<cplx> d = std::span<cplx>(dvr).subspan(block * c0, block * pass);
        passes.fbr_to_dvr(std::span<const cplx>(fbr).subspan(n_rot * c0, n_rot * pass), d);
        ComplexVector w(d.begin(), d.end());
        passes.dvr_to_fbr(w, std::span<cplx>(back).subspan(n_rot * c0, n_rot * pass));
      }
      expect_same_transform(dvr, dvr_one);
      expect_same_transform(back, back_one);
    }
  }
}

// FBR -> DVR against the earlier code: the row-by-row loop over all columns (one theta row of one column
// at a time), the empty phi slots set to zero, then the same phi FFT. The five-row sums give the same
// bits: 3 theta points (the remainder loop only), 7 (one block of five plus a remainder), 15 (blocks only).
TEST(AngularTransform, FiveRowSumsEqualTheRowByRowLoop) {
  const FftwFft fft;
  constexpr int n_columns = 4;
  for (const int j_max : {4, 12, 28}) {
    SCOPED_TRACE(j_max);
    const AngularGrid grid(j_max, j_max, true, 0, 0);
    const RotationalBasis basis(grid, 0, 0);
    const AngularTransform transform(grid, basis, fft, n_columns);
    const int n_theta = transform.n_theta();
    const int n_phi = transform.n_phi();
    const int ld = transform.leading_dimension();
    const std::size_t n_rot = basis.size();
    const std::size_t block = transform.dvr_block_size();
    ComplexVector fbr(n_rot * n_columns, cplx(0.0, 0.0));
    for (std::size_t i = 0; i < fbr.size(); ++i)
      fbr[i] = cplx(std::cos(0.7 * double(i)), std::sin(1.3 * double(i) + 0.2));

    ComplexVector expected(block * n_columns, cplx(0.0, 0.0));
    int index = 0;
    for (const AngularTransform::MjBlock& mj_block : transform.mj_blocks()) {
      const int n_j = mj_block.fbr.size();
      for (int t = 0; t < n_theta; ++t) {
        for (int c = 0; c < n_columns; ++c) {
          cplx sum(0.0, 0.0);
          for (int j = 0; j < n_j; ++j)
            sum = sum + transform.matrix()[index + j] * fbr[mj_block.fbr.first + j + n_rot * c];
          expected[mj_block.phi_slot + ld * t + block * c] = sum;
        }
        index += n_j;
      }
    }
    for (int c = 0; c < n_columns; ++c)
      for (int t = 0; t < n_theta; ++t)
        for (int k = transform.mj_max() + 1; k <= n_phi - transform.mj_max() - 1; ++k)
          expected[k + ld * t + block * c] = 0.0;
    fft.plan_complex(FftLayout{n_phi, n_theta * n_columns, 1, n_phi + 1}, FftSign::Plus, 1.0)->execute(expected);

    ComplexVector dvr(block * n_columns, cplx(0.0, 0.0));
    transform.fbr_to_dvr(fbr, dvr);
    EXPECT_TRUE(std::equal(dvr.begin(), dvr.end(), expected.begin()));
  }
}

// The largest divisor of n_xy whose DVR blocks fit AngularTransform::cpu_pass_bytes (58 KiB), at least 1.
TEST(AngularTransform, ColumnsPerPass) {
  EXPECT_EQ(AngularTransform::columns_per_pass(144, 15, 60), 4);   // DZ: 915 values = 14 640 B per block
  EXPECT_EQ(AngularTransform::columns_per_pass(256, 15, 60), 4);   // ER, 16 x 16
  EXPECT_EQ(AngularTransform::columns_per_pass(4096, 15, 60), 4);  // 64 x 64
  EXPECT_EQ(AngularTransform::columns_per_pass(64, 3, 12), 64);    // Dai-Light 8 x 8: 39 values, all columns
  EXPECT_EQ(AngularTransform::columns_per_pass(16, 3, 8), 16);     // hard wall 4 x 4
  EXPECT_EQ(AngularTransform::columns_per_pass(1, 15, 60), 1);
  EXPECT_EQ(AngularTransform::columns_per_pass(16, 100, 100), 1);  // one block alone exceeds the limit
  EXPECT_EQ(AngularTransform::columns_per_pass(7, 3, 8), 7);       // prime n_xy, small blocks
  EXPECT_EQ(AngularTransform::columns_per_pass(7, 15, 60), 1);     // prime n_xy, large blocks
}

// Tests: the molecule-surface PES (src/crplib/molecule_surface.*, hexagonal_ and square_interaction.cpp,
// wrapper.*, qd_wrapper.*): bitwise equality with the Fortran library on synthetic data sets,
// impulse responses, exactness at the data points, the gas-phase switch, gradients, and a
// representable model reproduced off the grid.
//
// Data sets with CutsInterpolation = F (total energies) are compared to 4e-15 max(u, |value|), u the
// atomic unit (Hartree, Hartree/bohr) in the output units:
// the Fortran reduces their cuts with gfortran's vectorised libmvec cos, the C++ with std::cos
// (decision D4, docs/fortran_deviations.md); measured on the reference cases: |dV| <= 4.4e-16,
// |dF| <= 2.7e-15 at under 12 % of the points. A scratch build calling _ZGVbN2v_cos instead was
// bitwise equal everywhere. CutsInterpolation = T data sets are compared bit for bit.
#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "crp_synthetic.hpp"
#include "crplib/errors.hpp"
#include "crplib/molecule_surface.hpp"
#include "crplib/molecule_surface_detail.hpp"
#include "crplib/qd_wrapper.hpp"
#include "crplib/wrapper.hpp"
#include "fortran_reference.hpp"

using namespace crplib;
using crplib_test::DriverInput;
using crplib_test::DriverOutput;
using crplib_test::reference_data_set;
using crplib_test::reference_file;

namespace {

// The largest |V| of the cut data in Hartree, at least 1. With CutsInterpolation = F the data are
// reduced by the atom-surface PES with cos theta, cos phi from libmvec in the Fortran (decision D4):
// a few ulp of the data values, which the interpolation carries into every result.
double data_scale(const DataSet& data) {
  double largest = 1.0;
  for (const MoleculeCut& cut : data.molecule_cuts)
    for (double v : cut.v)
      largest = std::max(largest, std::abs(convert_energy(v, data.units.energy, EnergyUnit::hartree)));
  return largest;
}

// Bitwise, or within the D4 tolerance for CutsInterpolation = F data sets: 4e-15 of the larger of
// the value and `unit` (the atomic unit of the quantity in its units, times the data_scale).
void expect_fortran(double cpp, std::uint64_t fortran_bits, bool bitwise, const std::string& where, double unit = 1.0) {
  if (bitwise) {
    EXPECT_BITS(cpp, fortran_bits) << where;
  } else {
    const double fortran = std::bit_cast<double>(fortran_bits);
    EXPECT_NEAR(cpp, fortran, 4e-15 * std::max(unit, std::abs(fortran))) << where;
  }
}

std::string at_text(const Vec3& a, const Vec3& b) {
  return " atoms (" + std::to_string(a.x) + ", " + std::to_string(a.y) + ", " + std::to_string(a.z) + ") (" +
         std::to_string(b.x) + ", " + std::to_string(b.y) + ", " + std::to_string(b.z) + ")";
}

// Compares calculate_molecule_surface with the driver's "molecule" output of one point list.
void compare_molecule_points(const std::string& case_name, const std::string& list, Units units) {
  const auto output_file = reference_file(case_name + "/" + list + ".out");
  if (!std::filesystem::exists(output_file)) GTEST_SKIP() << output_file << " missing";
  const DataSet data = read_data_set(reference_data_set(case_name));
  const bool bitwise = data.cuts_interpolation;
  const double scale = data_scale(data);
  const MoleculeSurface pes(data);
  DriverInput in(reference_file(case_name + "/" + list + ".in"));
  DriverOutput out(output_file);
  const int n = in.integer();
  for (int i = 0; i < n; ++i) {
    const double za = in.real(), xa = in.real(), ya = in.real(), zb = in.real(), xb = in.real(), yb = in.real();
    const Vec3 a{xa, ya, za}, b{xb, yb, zb};
    const MoleculeSurfaceValue v = calculate_molecule_surface(pes, a, b, units);
    const auto value = out.next("V");
    const auto f = out.next_bits("F", 6);  // Fortran F = (dXa, dYa, dZa, dXb, dYb, dZb)
    const std::string where = at_text(a, b);
    const double energy_unit = convert_energy(1.0, EnergyUnit::hartree, units.energy) * scale;
    const double force_unit = convert_force(1.0, internal_units, units) * scale;
    expect_fortran(v.v, value.bits, bitwise, where, energy_unit);
    expect_fortran(v.gradient_a.x, f[0], bitwise, where, force_unit);
    expect_fortran(v.gradient_a.y, f[1], bitwise, where, force_unit);
    expect_fortran(v.gradient_a.z, f[2], bitwise, where, force_unit);
    expect_fortran(v.gradient_b.x, f[3], bitwise, where, force_unit);
    expect_fortran(v.gradient_b.y, f[4], bitwise, where, force_unit);
    expect_fortran(v.gradient_b.z, f[5], bitwise, where, force_unit);
    if (units == internal_units) {  // S3: V alone is bitwise the full evaluation's V (0: silent branch)
      EXPECT_EQ(std::bit_cast<std::uint64_t>(pes.value(a, b).value_or(0.00)), std::bit_cast<std::uint64_t>(v.v))
          << " V only " << where;
    }
  }
}

// The QD wrapper (V in eV at skewed coordinates, clamps, shift) and its parameters at the points
// of <case>/wrapper.in.
void compare_qd_wrapper_points(const std::string& case_name) {
  const auto output_file = reference_file(case_name + "/wrapper.out");
  REQUIRE_CRP_REFERENCE(output_file);
  const DataSet data = read_data_set(reference_data_set(case_name));
  const QdWrapper wrapper(data);
  DriverOutput out(output_file);
  const auto records = crplib_test::read_driver_output(output_file);
  EXPECT_EQ(records.at(0).text, "point_group " + wrapper.parameters().point_group);
  (void)out.next("point_group");
  const auto params = out.next_bits("params", 4);
  EXPECT_BITS(wrapper.parameters().lattice_constant, params[0]);
  EXPECT_BITS(wrapper.parameters().lattice_constant, params[1]);
  EXPECT_BITS(wrapper.parameters().z_asymptotic, params[2]);
  EXPECT_BITS(wrapper.parameters().skewing_angle, params[3]);
  DriverInput in(reference_file(case_name + "/wrapper.in"));
  const int n = in.integer();
  for (int i = 0; i < n; ++i) {
    const double u = in.real(), v = in.real(), z = in.real(), r = in.real(), theta = in.real(), phi = in.real();
    expect_fortran(wrapper.potential_ev(u, v, z, r, theta, phi), out.next("V").bits, data.cuts_interpolation,
                   " u, v, Z, r, theta, phi = " + std::to_string(u) + ", " + std::to_string(v) + ", " +
                       std::to_string(z) + ", " + std::to_string(r) + ", " + std::to_string(theta) + ", " +
                       std::to_string(phi),
                   convert_energy(1.0, EnergyUnit::hartree, EnergyUnit::ev) * data_scale(data));
  }
}

}  // namespace

// ---------------------------------------------------------------------------------------------
// Bitwise equality with the Fortran library (tests/crplib/reference/<case>/)
// ---------------------------------------------------------------------------------------------

class FortranMoleculeSurface : public testing::TestWithParam<std::string> {};

TEST_P(FortranMoleculeSurface, MoleculeInternalUnits) {
  compare_molecule_points(GetParam(), "molecule_0_2", internal_units);
}

TEST_P(FortranMoleculeSurface, QdWrapper) { compare_qd_wrapper_points(GetParam()); }

TEST_P(FortranMoleculeSurface, GasPhase) {
  const auto output_file = reference_file(GetParam() + "/gasphase.out");
  REQUIRE_CRP_REFERENCE(output_file);
  const MoleculeSurface pes(read_data_set(reference_data_set(GetParam())));
  DriverOutput out(output_file);
  const auto zso = out.next_bits("zso", 2);
  EXPECT_BITS(pes.zso_min(), zso[0]);
  EXPECT_BITS(pes.zso_max(), zso[1]);
  DriverInput in(reference_file(GetParam() + "/gasphase.in"));
  const int n = in.integer();
  for (int i = 0; i < n; ++i) {
    const double z = in.real(), r = in.real();
    const auto f = out.next_bits("G", 5);  // Fso, dFso, V2D, dV2D/dZ, dV2D/dr
    const MoleculeSurface::Switch s = pes.switch_function(z);
    const MoleculeSurface::GasPhaseValue g = pes.gas_phase(z, r);
    EXPECT_BITS(s.f, f[0]) << " Z = " << z;
    EXPECT_BITS(s.df_dz, f[1]) << " Z = " << z;
    EXPECT_BITS(g.v, f[2]) << " Z, r = " << z << ", " << r;
    EXPECT_BITS(g.d_z, f[3]) << " Z, r = " << z << ", " << r;
    EXPECT_BITS(g.d_r, f[4]) << " Z, r = " << z << ", " << r;
  }
}

TEST_P(FortranMoleculeSurface, Cambcord) {
  const auto output_file = reference_file(GetParam() + "/cambcord.out");
  REQUIRE_CRP_REFERENCE(output_file);
  DriverInput in(reference_file(GetParam() + "/cambcord.in"));
  DriverOutput out(output_file);
  const int n = in.integer();
  for (int i = 0; i < n; ++i) {
    const Vec3 a{in.real(), in.real(), in.real()};
    const Vec3 b{in.real(), in.real(), in.real()};
    const auto f = out.next_bits("C", 6);
    const detail::MoleculeGeometry g = detail::molecule_geometry(a, b);
    EXPECT_BITS(g.x, f[0]);
    EXPECT_BITS(g.y, f[1]);
    EXPECT_BITS(g.z, f[2]);
    EXPECT_BITS(g.r, f[3]);
    EXPECT_BITS(g.theta, f[4]);
    EXPECT_BITS(g.phi, f[5]);
  }
}

INSTANTIATE_TEST_SUITE_P(Cases, FortranMoleculeSurface,
                         testing::Values("c6v_pairwise", "c3v_pairwise", "c4v_pairwise", "c6v_pairwise_angstrom",
                                         "c4v_negative_spacing", "c3v_missing_extra_cut", "cu111_srp48"));

// The data set in Angstrom/eV, evaluated in Angstrom/eV (the unit conversions of the wrapper).
TEST(FortranMoleculeSurfaceUnits, AngstromAndEv) {
  compare_molecule_points("c6v_pairwise_angstrom", "molecule_1_3", Units{LengthUnit::angstrom, EnergyUnit::ev});
}

namespace {

// The atom-surface PES of a case in Angstrom and eV (the unit conversions of the wrapper).
void compare_atom_points_in_angstrom_and_ev(const std::string& case_name) {
  const auto output_file = reference_file(case_name + "/atom_1_3.out");
  REQUIRE_CRP_REFERENCE(output_file);
  const auto pes = make_atom_surface(read_data_set(reference_data_set(case_name)));
  DriverInput in(reference_file(case_name + "/atom_1_3.in"));
  DriverOutput out(output_file);
  const int n = in.integer();
  const Units units{LengthUnit::angstrom, EnergyUnit::ev};
  for (int i = 0; i < n; ++i) {
    const Vec3 p{in.real(), in.real(), in.real()};
    const AtomSurfaceValue v = calculate_atom_surface(*pes, p, units);
    EXPECT_BITS(v.v, out.next("V").bits);
    const auto f = out.next_bits("F", 3);
    EXPECT_BITS(v.gradient.x, f[0]);
    EXPECT_BITS(v.gradient.y, f[1]);
    EXPECT_BITS(v.gradient.z, f[2]);
  }
}

}  // namespace

TEST(FortranMoleculeSurfaceUnits, AtomSurfaceInAngstromAndEv) {
  compare_atom_points_in_angstrom_and_ev("c6v_pairwise_angstrom");
}

// The real data set (cu111_srp48) is itself in Angstrom and eV.
TEST(FortranMoleculeSurfaceUnits, RealDataSetInAngstromAndEv) {
  compare_molecule_points("cu111_srp48", "molecule_1_3", Units{LengthUnit::angstrom, EnergyUnit::ev});
}

TEST(FortranMoleculeSurfaceUnits, RealDataSetAtomSurfaceInAngstromAndEv) {
  compare_atom_points_in_angstrom_and_ev("cu111_srp48");
}

// Impulse data sets (one cut 1, all else 0): bitwise V and gradient at the Fourier sites and their
// rotated images, any orientation.
class FortranImpulse : public testing::TestWithParam<std::string> {};

TEST_P(FortranImpulse, BitwiseEqual) { compare_molecule_points(GetParam(), "molecule_0_2", internal_units); }

// The QD wrapper bitwise (with CutsInterpolation = F the comparison allows a few ulp, decision D4,
// which would hide e.g. a wrong skew constant).
TEST_P(FortranImpulse, QdWrapper) { compare_qd_wrapper_points(GetParam()); }

INSTANTIATE_TEST_SUITE_P(Cuts, FortranImpulse,
                         testing::Values("c6v_impulse_1", "c6v_impulse_5", "c6v_impulse_10", "c6v_impulse_15",
                                         "c3v_impulse_2", "c3v_impulse_6", "c3v_impulse_9", "c3v_impulse_15",
                                         "c3v_impulse_21", "c3v_impulse_27", "c4v_impulse_2", "c4v_impulse_6",
                                         "c4v_impulse_9", "c4v_impulse_13"));

// ---------------------------------------------------------------------------------------------
// Properties (no Fortran needed)
// ---------------------------------------------------------------------------------------------

namespace {

// The Fourier sites in Psite order: position (lattice units, skewed for C3v/C6v, Cartesian for
// C4v), the cut-site type (index of its first cut) and the azimuth offset at which the site function
// is evaluated (degrees). Appendix A of docs/crp_library.md.
struct FourierSite {
  double u, v;
  int first_cut;
  double offset;
};

std::vector<FourierSite> fourier_sites(Symmetry s) {
  constexpr double t = 1.0 / 3, h = 1.0 / 6;
  if (s == Symmetry::c4v) {
    return {{0, 0, 0, 0},    {0.5, 0, 3, 0},        {0.5, 0.5, 7, 0},       {0.25, 0.25, 10, 0},
            {0, 0.5, 3, 90}, {0.75, 0.25, 10, -90}, {0.75, 0.75, 10, -180}, {0.25, 0.75, 10, 90}};
  }
  if (s == Symmetry::c6v) {
    return {{0, 0, 0, 0},           {0.5, 0, 3, 0},          {t, t, 7, 0},         {h, h, 12, 0},
            {0.5, 0.5, 3, 60},      {0, 0.5, 3, -60},        {2 * t, h, 12, -120}, {h, 2 * t, 12, 120},
            {2 * t, 2 * t, 7, 180}, {5 * h, 5 * h, 12, 180}, {t, 5 * h, 12, 60},   {5 * h, t, 12, -60}};
  }
  return {{0, 0, 0, 0},          {0.5, 0, 3, 0},          {t, t, 7, 0},         {h, h, 12, 0},
          {0.5, 0.5, 3, 60},     {0, 0.5, 3, -60},        {2 * t, h, 12, -120}, {h, 2 * t, 12, 120},
          {2 * t, 2 * t, 24, 0}, {5 * h, 5 * h, 18, 120}, {t, 5 * h, 18, 0},    {5 * h, t, 18, -120}};
}

// First cut of the site that cut k belongs to.
int site_of_cut(Symmetry s, int k) {
  const std::vector<int> starts = s == Symmetry::c4v   ? std::vector<int>{0, 3, 7, 10}
                                  : s == Symmetry::c6v ? std::vector<int>{0, 3, 7, 12}
                                                       : std::vector<int>{0, 3, 7, 12, 18, 24};
  int first = 0;
  for (int start : starts)
    if (k >= start) first = start;
  return first;
}

Vec3 centre(Symmetry s, double u, double v, double a) {
  if (s == Symmetry::c4v) return {u * a, v * a, 0.0};
  return {(u + v / 2.0) * a, v * std::sqrt(3.0) / 2.0 * a, 0.0};
}

std::pair<Vec3, Vec3> atoms(const Vec3& c, double z, double r, double theta, double phi) {
  const auto [a, b] = crplib_test::atom_pair(c.x, c.y, z, r, theta, phi);
  return {Vec3{a.x, a.y, a.z}, Vec3{b.x, b.y, b.z}};
}

constexpr double deg = std::numbers::pi / 180.0;

}  // namespace

// Impulse on cut k: the interpolation is 1 at every Fourier site of k's site type when the molecule
// has k's orientation in that site's frame (the rotated images included), and 0 at the sites of the
// other types. All cuts of all symmetries.
TEST(MoleculeSurfaceImpulse, ResponseAtTheFourierSites) {
  for (Symmetry s : {Symmetry::c3v, Symmetry::c4v, Symmetry::c6v}) {
    const auto geometry = crplib_test::molecule_cut_geometry(s);
    for (int k = 0; k < molecule_cut_count(s); ++k) {
      const MoleculeSurface pes(crplib_test::impulse_data_set(s, k, crplib_test::impulse_grids()));
      for (const FourierSite& site : fourier_sites(s)) {
        const double theta = geometry[k].theta * deg;
        const double phi = (geometry[k].phi - site.offset) * deg;
        const auto [a, b] = atoms(centre(s, site.u, site.v, pes.lattice_constant()), 2.0, 1.4, theta, phi);
        const auto value = pes.evaluate(a, b);
        ASSERT_TRUE(value.has_value());
        const double expected = site.first_cut == site_of_cut(s, k) ? 1.0 : 0.0;
        EXPECT_NEAR(value->v, expected, 1e-14)
            << symmetry_name(s) << " cut " << k + 1 << " site (" << site.u << ", " << site.v << ")";
      }
    }
  }
}

class MoleculeSurfaceProperties : public testing::TestWithParam<Symmetry> {
 protected:
  void SetUp() override {
    data_ = crplib_test::pairwise_data_set(GetParam());
    pes_ = std::make_unique<MoleculeSurface>(data_);
  }
  DataSet data_;
  std::unique_ptr<MoleculeSurface> pes_;
};

// CutsInterpolation = F: V reproduces the total energies at the cut knots (below ZsoMin). Interior
// knots only: on the grid edges the r and Z recomputed from the atoms (CAMBCORD) may fall an ulp
// outside the grid (the QD wrapper clamps with a margin for this reason).
TEST_P(MoleculeSurfaceProperties, ExactAtTheCutKnots) {
  const auto geometry = crplib_test::molecule_cut_geometry(GetParam());
  for (int k = 0; k < molecule_cut_count(GetParam()); ++k) {
    const MoleculeCut& cut = data_.molecule_cuts[k];
    const auto [x, y] = crplib_test::cut_site_cartesian(GetParam(), geometry[k], data_.lattice_constant);
    for (std::size_t j = 1; j + 1 < cut.z.size(); j += 3) {
      if (cut.z[j] > data_.zso_min) continue;
      for (std::size_t i = 1; i + 1 < cut.r.size(); i += 2) {
        const auto [a, b] = atoms({x, y, 0.0}, cut.z[j], cut.r[i], geometry[k].theta * deg, geometry[k].phi * deg);
        const auto v = pes_->evaluate(a, b);
        ASSERT_TRUE(v.has_value());
        const double expected = cut.v[i + cut.r.size() * j];
        EXPECT_NEAR(v->v, expected, 1e-12 * std::max(1.0, std::abs(expected)))
            << "cut " << k + 1 << " r, Z = " << cut.r[i] << ", " << cut.z[j];
      }
    }
  }
}

// The same PES from total energies (F) and from the reduced interpolation function (T).
TEST_P(MoleculeSurfaceProperties, CutsInterpolationTrueAndFalseAgree) {
  DataSet reduced = data_;
  reduced.cuts_interpolation = true;
  const auto atom = make_atom_surface(data_);
  const auto geometry = crplib_test::molecule_cut_geometry(GetParam());
  for (std::size_t k = 0; k < reduced.molecule_cuts.size(); ++k) {
    MoleculeCut& cut = reduced.molecule_cuts[k];
    const auto [x, y] = crplib_test::cut_site_cartesian(GetParam(), geometry[k], data_.lattice_constant);
    for (std::size_t j = 0; j < cut.z.size(); ++j) {
      for (std::size_t i = 0; i < cut.r.size(); ++i) {
        const auto [a, b] = atoms({x, y, 0.0}, cut.z[j], cut.r[i], geometry[k].theta * deg, geometry[k].phi * deg);
        cut.v[i + cut.r.size() * j] -= atom->evaluate(a).v + atom->evaluate(b).v;
      }
    }
  }
  const MoleculeSurface from_reduced(reduced);
  for (const auto& [u, v, z, r, theta, phi] : {std::tuple{0.13, 0.29, 1.7, 1.3, 0.4, 0.9},
                                               {0.71, 0.66, 3.1, 0.9, 1.4, -2.0},
                                               {0.5, 0.1, 5.8, 2.0, 0.2, 2.2}}) {
    const auto [a, b] = atoms(centre(GetParam(), u, v, data_.lattice_constant), z, r, theta, phi);
    EXPECT_NEAR(pes_->evaluate(a, b)->v, from_reduced.evaluate(a, b)->v, 1e-13);
  }
}

// Above ZsoMax only the gas phase: VH2(r) + Zdep(Z), independent of the lateral position and the
// orientation (to rounding: the centre and r are recomputed from the atoms).
TEST_P(MoleculeSurfaceProperties, GasPhaseAboveZsoMax) {
  for (double z : {6.6, 8.0, 12.0}) {
    for (double r : {0.8, 1.4, 3.0}) {
      std::vector<double> values;
      for (const auto& [u, v, theta, phi] :
           {std::tuple{0.0, 0.0, 0.0, 0.0}, {0.31, 0.57, 1.1, 2.0}, {0.9, 0.2, 2.5, -0.7}}) {
        const auto [a, b] = atoms(centre(GetParam(), u, v, data_.lattice_constant), z, r, theta, phi);
        values.push_back(pes_->evaluate(a, b)->v);
      }
      EXPECT_NEAR(values[1], values[0], 1e-14);
      EXPECT_NEAR(values[2], values[0], 1e-14);
      EXPECT_NEAR(values[0], pes_->gas_phase(z, r).v, 1e-14);
    }
  }
}

// V is continuous across ZsoMin and ZsoMax (the switch function goes from 1 to 0 smoothly).
TEST_P(MoleculeSurfaceProperties, ContinuousAcrossTheSwitchingLimits) {
  for (double limit : {data_.zso_min, data_.zso_max}) {
    const double e = 1e-8;
    const auto [a1, b1] = atoms(centre(GetParam(), 0.23, 0.41, data_.lattice_constant), limit - e, 1.3, 0.8, 0.5);
    const auto [a2, b2] = atoms(centre(GetParam(), 0.23, 0.41, data_.lattice_constant), limit + e, 1.3, 0.8, 0.5);
    EXPECT_NEAR(pes_->evaluate(a1, b1)->v, pes_->evaluate(a2, b2)->v, 1e-7);
  }
}

TEST_P(MoleculeSurfaceProperties, GradientMatchesFiniteDifferences) {
  const double h = 1e-5;
  // below ZsoMin, in the switching region, above ZsoMax; a nearly vertical molecule
  for (const auto& [u, v, z, r, theta, phi] : {std::tuple{0.13, 0.29, 1.7, 1.3, 0.4, 0.9},
                                               {0.71, 0.66, 3.1, 0.9, 1.4, -2.0},
                                               {0.4, 0.1, 5.7, 2.0, 0.9, 2.2},
                                               {0.6, 0.35, 8.0, 1.5, 1.1, 0.3},
                                               {0.2, 0.2, 2.5, 1.4, 1e-3, 0.4}}) {
    const auto [a, b] = atoms(centre(GetParam(), u, v, data_.lattice_constant), z, r, theta, phi);
    const MoleculeSurfaceValue value = *pes_->evaluate(a, b);
    auto shifted = [&](int atom, int axis, double step) {
      Vec3 p = a, q = b;
      Vec3& m = atom == 0 ? p : q;
      (axis == 0 ? m.x : axis == 1 ? m.y : m.z) += step;
      return pes_->evaluate(p, q)->v;
    };
    const double analytic[2][3] = {{value.gradient_a.x, value.gradient_a.y, value.gradient_a.z},
                                   {value.gradient_b.x, value.gradient_b.y, value.gradient_b.z}};
    for (int atom = 0; atom < 2; ++atom) {
      for (int axis = 0; axis < 3; ++axis) {
        const double numeric = (shifted(atom, axis, h) - shifted(atom, axis, -h)) / (2 * h);
        EXPECT_NEAR(analytic[atom][axis], numeric, 1e-6 * std::max(1.0, std::abs(numeric)))
            << symmetry_name(GetParam()) << " atom " << atom << " axis " << axis << " at Z = " << z
            << ", theta = " << theta;
      }
    }
  }
}

// Fortran quirk, kept (docs/fortran_deviations.md): for a vertical molecule (theta = 0, rho = 0) the
// Fortran sets d theta / dX = d theta / dY = 0, but the sin(2 theta) terms of the hollow and t2 sites
// make V(dx) linear in |dx| times a sign that flips with phi -> phi + pi, so V has a nonzero lateral
// slope there. The z components are right; x and y miss that term.
TEST_P(MoleculeSurfaceProperties, LateralGradientOfAVerticalMoleculeIsTheFortranOne) {
  const auto [a, b] = atoms(centre(GetParam(), 0.2, 0.2, data_.lattice_constant), 2.5, 1.4, 0.0, 0.0);
  const MoleculeSurfaceValue value = *pes_->evaluate(a, b);
  const double h = 1e-5;
  auto v_at = [&](double dx, double dz) {
    Vec3 p = a;
    p.x += dx;
    p.z += dz;
    return pes_->evaluate(p, b)->v;
  };
  EXPECT_NEAR(value.gradient_a.z, (v_at(0, h) - v_at(0, -h)) / (2 * h), 1e-6);
  const double numeric_x = (v_at(h, 0) - v_at(-h, 0)) / (2 * h);
  EXPECT_GT(std::abs(value.gradient_a.x - numeric_x), 1e-6) << "the theta = 0 quirk is gone: update the deviation note";
}

TEST_P(MoleculeSurfaceProperties, SilentBranchOutsideTheCuts) {
  // r beyond the cuts below ZsoMax, and below the free-molecule potential above it
  const auto [a1, b1] = atoms(centre(GetParam(), 0.2, 0.3, data_.lattice_constant), 2.0, 3.4, 0.5, 0.5);
  EXPECT_FALSE(pes_->evaluate(a1, b1).has_value());
  EXPECT_FALSE(pes_->covers(3.4, 2.0));
  const auto [a2, b2] = atoms(centre(GetParam(), 0.2, 0.3, data_.lattice_constant), 8.0, 0.3, 0.5, 0.5);
  EXPECT_FALSE(pes_->evaluate(a2, b2).has_value());
  EXPECT_TRUE(pes_->covers(1.4, 2.0));
  const MoleculeSurfaceValue zero = calculate_molecule_surface(*pes_, a1, b1, internal_units);
  EXPECT_EQ(zero.v, 0.0);
  EXPECT_EQ(zero.gradient_b.z, 0.0);
}

INSTANTIATE_TEST_SUITE_P(Symmetries, MoleculeSurfaceProperties,
                         testing::Values(Symmetry::c3v, Symmetry::c4v, Symmetry::c6v));

// A data set in the span of the interpolation: I = c(Z, r) (A + B cos 2 theta) S(x, y), c a bicubic
// polynomial, S invariant under the plane group (C6v: cosine stars; C3v: plus the sine star; C4v:
// a p4mm combination), everything else 0. The PES must reproduce it at arbitrary points.
TEST(MoleculeSurfaceRepresentable, ReproducedOffTheGrid) {
  auto c = [](double z, double r) {
    return 0.01 * (1 + 0.2 * z - 0.03 * z * z + 0.002 * z * z * z) * (0.5 - 0.4 * r + 0.15 * r * r - 0.02 * r * r * r);
  };
  auto angular = [](double theta) { return 1.0 + 0.4 * std::cos(2.0 * theta); };
  for (Symmetry s : {Symmetry::c3v, Symmetry::c4v, Symmetry::c6v}) {
    DataSet data = crplib_test::impulse_data_set(s, 0);  // zero atom PES and gas phase, default grids
    const double a = data.lattice_constant;
    auto lateral = [&](double x, double y) {
      if (s == Symmetry::c4v) {
        const double f = 2 * std::numbers::pi / a;
        return 1 + 0.3 * (std::cos(f * x) + std::cos(f * y)) + 0.1 * std::cos(f * (x + y)) +
               0.1 * std::cos(f * (x - y)) + 0.05 * (std::cos(2 * f * x) + std::cos(2 * f * y));
      }
      const double v = y / (a * std::sqrt(3.0) / 2.0), u = x / a - v / 2.0;
      const double X = 2 * std::numbers::pi * u, Y = 2 * std::numbers::pi * v;
      double value = 1 + 0.3 * (std::cos(X) + std::cos(Y) + std::cos(X + Y)) +
                     0.1 * (std::cos(2 * X) + std::cos(2 * Y) + std::cos(2 * (X + Y)));
      if (s == Symmetry::c3v) value += 0.2 * (-std::sin(X) - std::sin(Y) + std::sin(X + Y));
      return value;
    };
    const auto geometry = crplib_test::molecule_cut_geometry(s);
    for (std::size_t k = 0; k < data.molecule_cuts.size(); ++k) {
      MoleculeCut& cut = data.molecule_cuts[k];
      const auto [x, y] = crplib_test::cut_site_cartesian(s, geometry[k], a);
      for (std::size_t j = 0; j < cut.z.size(); ++j)
        for (std::size_t i = 0; i < cut.r.size(); ++i)
          cut.v[i + cut.r.size() * j] = c(cut.z[j], cut.r[i]) * angular(geometry[k].theta * deg) * lateral(x, y);
    }
    const MoleculeSurface pes(data);
    for (const auto& [u, v, z, r, theta, phi] : {std::tuple{0.137, 0.291, 1.73, 1.31, 0.42, 0.93},
                                                 {0.717, 0.659, 3.11, 0.93, 1.37, -2.04},
                                                 {-0.41, 1.27, 4.47, 2.21, 2.71, 0.31},
                                                 {0.5, 0.5, 0.97, 1.79, 0.0, 0.0}}) {
      const Vec3 m = centre(s, u, v, a);
      const auto [atom_a, atom_b] = atoms(m, z, r, theta, phi);
      const double expected = c(z, r) * angular(theta) * lateral(m.x, m.y);
      EXPECT_NEAR(pes.evaluate(atom_a, atom_b)->v, expected, 1e-13) << symmetry_name(s) << " at " << u << ", " << v;
    }
  }
}

// ---------------------------------------------------------------------------------------------
// Speed-up S4: one grid lookup per group of molecular cuts with the same (r, Z) grid
// ---------------------------------------------------------------------------------------------

namespace {

bool same_bits(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }

bool same_bits(const std::optional<MoleculeSurfaceValue>& a, const std::optional<MoleculeSurfaceValue>& b) {
  if (a.has_value() != b.has_value()) return false;
  if (!a) return true;
  for (auto [p, q] : {std::pair{a->gradient_a, b->gradient_a}, std::pair{a->gradient_b, b->gradient_b}})
    if (!same_bits(p.x, q.x) || !same_bits(p.y, q.y) || !same_bits(p.z, q.z)) return false;
  return same_bits(a->v, b->v);
}

// The molecule-surface PES with and without the shared lookup at points over the cell, Z from
// below Zmin to above ZsoMax and r beyond both ends of the cut grids (V and gradients, bitwise;
// outside a grid both give nullopt). Returns the number of points and of nullopt results.
std::pair<int, int> compare_shared_lookup(const DataSet& data) {
  const MoleculeSurface shared(data);
  const MoleculeSurface per_cut(data, EvaluationOptions{.share_grid_lookup = false});
  const double a = data.lattice_constant;
  int n = 0, outside = 0;
  for (double u : {0.0, 0.31, 0.5, 0.77})
    for (double v : {0.0, 0.18, 0.5, 0.93})
      for (double z = 0.2; z < 7.0; z += 0.23)
        for (double r : {0.5, 0.6, 0.83, 1.4, 2.2, 3.2, 3.4})
          for (const auto& [theta, phi] : {std::pair{0.0, 0.0}, {0.6, 1.9}, {1.5, -0.7}}) {
            const auto [atom_a, atom_b] = crplib_test::atom_pair((u + 0.5 * v) * a, 0.866 * v * a, z, r, theta, phi);
            const Vec3 pa{atom_a.x, atom_a.y, atom_a.z}, pb{atom_b.x, atom_b.y, atom_b.z};
            const std::optional<MoleculeSurfaceValue> s = shared.evaluate(pa, pb), p = per_cut.evaluate(pa, pb);
            EXPECT_TRUE(same_bits(s, p)) << symmetry_name(data.symmetry) << " at Z " << z << ", r " << r;
            ++n;
            if (!s) ++outside;
          }
  return {n, outside};
}

}  // namespace

// The synthetic data sets: all cuts share one grid.
TEST(MoleculeSurfaceSharedGrid, BitwiseAsPerCut) {
  for (Symmetry s : {Symmetry::c3v, Symmetry::c4v, Symmetry::c6v}) {
    const DataSet data = crplib_test::pairwise_data_set(s);
    EXPECT_EQ(MoleculeSurface(data).cut_grids(), 1);
    const auto [n, outside] = compare_shared_lookup(data);
    EXPECT_GT(outside, 0);  // r and Z outside the grids are part of the test
    EXPECT_GT(n - outside, 1000);
  }
}

// One cut with its own grid (Z ending below ZsoMax: points between give nullopt through that cut
// only), and one whose grid differs from the others only in the sign of a zero knot.
TEST(MoleculeSurfaceSharedGrid, CutsWithTheirOwnGrid) {
  DataSet data = crplib_test::pairwise_data_set(Symmetry::c6v);
  MoleculeCut& short_cut = data.molecule_cuts[5];
  short_cut.z.resize(short_cut.z.size() - 2);
  short_cut.v.resize(short_cut.r.size() * short_cut.z.size());
  ASSERT_LT(short_cut.z.back(), data.zso_max);
  EXPECT_EQ(MoleculeSurface(data).cut_grids(), 2);
  const auto [n, outside] = compare_shared_lookup(data);
  EXPECT_GT(outside, 0);

  // One interior Z knot of cut 7 moved: its cells and weights differ from the others' there.
  DataSet moved = crplib_test::pairwise_data_set(Symmetry::c6v);
  std::vector<double>& z7 = moved.molecule_cuts[7].z;
  z7[5] = z7[5] + 0.3 * (z7[6] - z7[5]);
  EXPECT_EQ(MoleculeSurface(moved).cut_grids(), 2);
  (void)compare_shared_lookup(moved);

  // A Z grid through zero: +0 in every cut, -0 in cut 9 (equal by ==, but not the same bits). An
  // impulse data set (cut 9 = 1; CutsInterpolation = T, atom cuts over z = -10..20).
  DataSet zero = crplib_test::impulse_data_set(Symmetry::c6v, 9, crplib_test::impulse_grids());
  for (MoleculeCut& cut : zero.molecule_cuts) {
    const double shift = cut.z[1];
    for (double& z : cut.z) z -= shift;  // the second knot at exactly +0
  }
  zero.molecule_cuts[9].z[1] = -0.0;
  EXPECT_EQ(MoleculeSurface(zero).cut_grids(), 2);
  (void)compare_shared_lookup(zero);
  (void)n;
}

// An extra cut limits the ranges (Rmin, Rmax, Zmin) unless its file is missing: then it is
// skipped, as in the Fortran (Read_MoleculeSurface_Data, INQUIRE EXIST).
TEST(MoleculeSurfaceRanges, ExtraCutsAndMissingFiles) {
  const MoleculeSurface plain(crplib_test::pairwise_data_set(Symmetry::c3v));
  DataSet d = crplib_test::case_data_set("c3v_missing_extra_cut");
  const MoleculeSurface narrowed(d);
  EXPECT_EQ(narrowed.r_max(), d.molecule_cuts[29].r.back());  // the present extra cut (bohr data)
  EXPECT_LT(narrowed.r_max(), plain.r_max());
  EXPECT_EQ(narrowed.r_min(), plain.r_min());
  EXPECT_EQ(narrowed.z_min(), plain.z_min());
  d.molecule_cuts[29].present = false;  // both extra cuts without a file
  EXPECT_EQ(MoleculeSurface(d).r_max(), plain.r_max());
}

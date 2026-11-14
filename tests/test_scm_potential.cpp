// Tests: the SCM potential (scm_crp_potential.*: factory, system parameters, V = CRP(u/s, v/s) + ΔV) on the
// SRP48 data set with its example files, and SCM D9 (scm_diagnostics.*) on slabs built in memory. Against the
// Fortran: tests/scm/test_scm_potential_reference.cpp.
#include <gtest/gtest.h>

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "potential/potential_factory.hpp"
#include "potential/scm/scm_correction.hpp"
#include "potential/scm/scm_crp_potential.hpp"
#include "potential/scm/scm_diagnostics.hpp"
#include "potential/scm/scm_files.hpp"

using namespace spodvr;
namespace fs = std::filesystem;

namespace {

std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

const fs::path example_directory =
    fs::path(SPODVR_REFERENCE_DIR).parent_path().parent_path() / "potentials/H2Cu111SRP_C3v.dir";
const fs::path srp48 = example_directory / "H2Cu111SRP_C3v";

// The example SCM files of the SRP48 data set, configuration 1, cut-off 16 bohr.
ScmInput example_input(bool rescale) {
  ScmInput input;
  input.parameter_file = (example_directory / "params.dat").string();
  input.position_file = (example_directory / "examplePosfile.dat").string();
  input.rescale = rescale;
  return input;
}

// ---------------------------------------------------------------------------------------------
// The lattice check (SCM D9) on slabs built in memory
// ---------------------------------------------------------------------------------------------

constexpr double pi = std::numbers::pi;
constexpr double a = 4.9;        // bohr, the lattice constant of the CRP
constexpr double spacing = 3.9;  // bohr, its first interlayer spacing
constexpr double third = 1.0 / 3.0;

SystemParameters crp_parameters() { return {a, a, 10.0, pi / 3.0, "C3v"}; }

// A C3v summary: layer 2 under (1/3, 1/3), one spacing below; with c6v, a C6v one (no second layer).
CrpDataSetSummary data_set_summary(bool c6v = false) {
  CrpDataSetSummary summary;
  summary.symmetry = c6v ? "C6v" : "C3v";
  summary.interlayer_spacings = {spacing, spacing, spacing};
  if (!c6v) summary.second_layer_site = std::array<double, 2>{third, third};
  return summary;
}

// The ideal (and undistorted) position of site (i + offset, j + offset)·a at height z.
Vec3 site(double i, double j, double offset, double z) {
  const double u = (i + offset) * a, v = (j + offset) * a;
  return {u + std::cos(pi / 3.0) * v, std::sin(pi / 3.0) * v, z};
}

// A 3 × 3 slab of two layers on the C3v lattice: layer 1 (atoms 0–8) at the sites (i, j), layer 2
// (atoms 9–17) under `second_offset` (1/3: the CRP's stacking) at z = −spacing.
ScmSlab slab_on_lattice(double second_offset = third) {
  ScmSlab slab;
  slab.configuration = 1;
  for (int layer = 1; layer <= 2; ++layer)
    for (int i = -1; i <= 1; ++i)
      for (int j = -1; j <= 1; ++j) {
        const Vec3 p = layer == 1 ? site(i, j, 0.0, 0.0) : site(i, j, second_offset, -spacing);
        slab.atoms.push_back({p, p, layer});
      }
  return slab;
}

// The InputError of require_slab_lattice, or "" if it accepts the slab.
std::string lattice_error(const ScmSlab& slab, const CrpDataSetSummary& summary = data_set_summary()) {
  const ScmLatticeCheck check = check_slab_lattice(slab, crp_parameters(), summary);
  try {
    require_slab_lattice(check, slab, crp_parameters());
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), "PosFilepath");
    return e.what();
  }
  return "";
}

}  // namespace

// A slab whose ideal atoms are the CRP's sites passes, without a warning; both layers are checked.
TEST(ScmLatticeCheck, ScmD9IdealSlabPasses) {
  const ScmLatticeCheck check = check_slab_lattice(slab_on_lattice(), crp_parameters(), data_set_summary());
  EXPECT_EQ(check.top_layer_atoms, 9);
  EXPECT_EQ(check.second_layer_atoms, 9);
  EXPECT_LE(check.max_deviation(), 1e-14);  // rounding of the site positions only
  EXPECT_FALSE(check.warning());
  EXPECT_EQ(lattice_error(slab_on_lattice()), "");
}

// One top-layer atom moved off its site by δ: no warning at 5e-6 bohr, a warning (accepted) at 2e-5
// and 5e-4 bohr, an error at 2e-3 bohr naming the moved atom; the deviation is δ.
TEST(ScmLatticeCheck, ScmD9ThresholdsBothSides) {
  for (const double delta : {5e-6, 2e-5, 5e-4, 2e-3}) {
    SCOPED_TRACE(delta);
    ScmSlab slab = slab_on_lattice();
    slab.atoms[4].ideal.x += delta;  // the atom at the origin
    const ScmLatticeCheck check = check_slab_lattice(slab, crp_parameters(), data_set_summary());
    EXPECT_NEAR(check.max_deviation(), delta, 1e-14);
    EXPECT_EQ(check.worst_atom, 4);
    EXPECT_EQ(check.warning(), delta > scm_lattice_warning_distance);
    const std::string error = lattice_error(slab);
    EXPECT_EQ(error.empty(), delta < scm_lattice_error_distance) << error;
  }
}

// The second layer of a C3v data set: under (1/3, 1/3) one spacing down. Mirrored stacking (2/3,
// 2/3: fcc and hcp swapped) lies a/√3 from every site, a wrong height its offset: both errors.
TEST(ScmLatticeCheck, ScmD9SecondLayerOfC3v) {
  const ScmSlab mirrored = slab_on_lattice(2.0 * third);
  const ScmLatticeCheck check = check_slab_lattice(mirrored, crp_parameters(), data_set_summary());
  EXPECT_NEAR(check.second_layer_deviation, a / std::sqrt(3.0), 1e-12);
  EXPECT_LE(check.top_layer_deviation, 1e-14);
  EXPECT_GE(check.worst_atom, 9);
  EXPECT_NE(lattice_error(mirrored).find("layer 2"), std::string::npos);

  ScmSlab low = slab_on_lattice();
  for (ScmAtom& atom : low.atoms)
    if (atom.layer == 2) atom.ideal.z -= 0.01;
  EXPECT_NEAR(check_slab_lattice(low, crp_parameters(), data_set_summary()).second_layer_deviation, 0.01, 1e-12);
  EXPECT_NE(lattice_error(low), "");
}

// A C6v data set describes one layer: the second layer is not checked (mirrored stacking passes),
// the top layer is.
TEST(ScmLatticeCheck, ScmD9C6vChecksTheTopLayerOnly) {
  const ScmSlab mirrored = slab_on_lattice(2.0 * third);
  const ScmLatticeCheck check = check_slab_lattice(mirrored, crp_parameters(), data_set_summary(true));
  EXPECT_EQ(check.top_layer_atoms, 9);
  EXPECT_EQ(check.second_layer_atoms, 0);
  EXPECT_EQ(lattice_error(mirrored, data_set_summary(true)), "");

  ScmSlab shifted = slab_on_lattice();
  shifted.atoms[0].ideal.y += 0.01;
  EXPECT_NE(lattice_error(shifted, data_set_summary(true)), "");
}

// A NaN position fails the check, also when atoms with finite deviations follow it.
TEST(ScmLatticeCheck, ScmD9NanFails) {
  ScmSlab slab = slab_on_lattice();
  slab.atoms[0].ideal.x = std::numeric_limits<double>::quiet_NaN();
  const ScmLatticeCheck check = check_slab_lattice(slab, crp_parameters(), data_set_summary());
  EXPECT_TRUE(std::isnan(check.max_deviation()));
  EXPECT_EQ(check.worst_atom, 0);
  EXPECT_TRUE(check.warning());
  EXPECT_NE(lattice_error(slab), "");
}

// ---------------------------------------------------------------------------------------------
// The potential on the SRP48 data set with the example slab
// ---------------------------------------------------------------------------------------------

// The example slab lies on the SRP48 lattice within 1.2e-6 bohr (its bohr constant differs from the
// CRP library's by 7.2e-8 relative; its farthest kept atom is 15.8 bohr from the origin): accepted
// without a warning; both layers are checked.
TEST(ScmCrpPotential, ScmD9ExampleSlabIsOnTheSrp48Lattice) {
  const ScmCrpPotential potential(srp48, example_input(true));
  const ScmLatticeCheck& check = potential.lattice_check();
  std::cout << "example slab: " << check.top_layer_atoms << " + " << check.second_layer_atoms
            << " atoms, largest deviation " << check.max_deviation() << " bohr\n";
  EXPECT_GT(check.top_layer_atoms, 0);
  EXPECT_GT(check.second_layer_atoms, 0);
  EXPECT_EQ(check.top_layer_atoms + check.second_layer_atoms, std::ssize(potential.slab().atoms));
  EXPECT_LE(check.max_deviation(), 1.2e-6);
  EXPECT_FALSE(check.warning());
}

// The factory selects the SCM potential for CRPDistrtSurface in any spelling (SCM D1).
TEST(ScmCrpPotential, ScmD1FactoryIsCaseInsensitive) {
  InputParameters p;
  p.system_id = "crpdistrtsurface";
  p.crp_dataset = srp48.string();
  p.scm = example_input(true);
  EXPECT_TRUE(potential_is_implemented(p.system_id));
  const std::unique_ptr<Potential> potential = make_potential(p);
  EXPECT_EQ(potential->name(), "CRPDistrtSurface");
  EXPECT_NO_THROW((void)dynamic_cast<const ScmCrpPotential&>(*potential));
}

// With rescaling both lattice constants are the slab's a_T and s = a_T/a₀; without, the system
// parameters are the CRP's. Point group, z_asymptotic and skew angle are the CRP's in both.
TEST(ScmCrpPotential, SystemParametersWithAndWithoutRescaling) {
  const ScmCrpPotential rescaled(srp48, example_input(true));
  const SystemParameters& crp = rescaled.crp().system_parameters();
  const SystemParameters& p = rescaled.system_parameters();
  const double a_t = rescaled.slab().header.lattice_constant;
  EXPECT_EQ(a_t, 4.8896663475);  // line 3 of examplePosfile.dat
  EXPECT_EQ(bits(p.lattice_constant_x), bits(a_t));
  EXPECT_EQ(bits(p.lattice_constant_y), bits(a_t));
  ASSERT_TRUE(rescaled.rescaling_factor().has_value());
  EXPECT_EQ(bits(*rescaled.rescaling_factor()), bits(a_t / crp.lattice_constant_x));
  EXPECT_EQ(bits(p.z_asymptotic), bits(crp.z_asymptotic));
  EXPECT_EQ(bits(p.skewing_angle), bits(crp.skewing_angle));
  EXPECT_EQ(p.point_group, "C3v");
  EXPECT_EQ(rescaled.name(), "CRPDistrtSurface");

  const ScmCrpPotential plain(srp48, example_input(false));
  const SystemParameters& q = plain.system_parameters();
  const SystemParameters& crp_plain = plain.crp().system_parameters();
  EXPECT_FALSE(plain.rescaling_factor().has_value());
  EXPECT_EQ(bits(q.lattice_constant_x), bits(crp_plain.lattice_constant_x));
  EXPECT_EQ(bits(q.lattice_constant_y), bits(crp_plain.lattice_constant_y));
  EXPECT_EQ(bits(q.z_asymptotic), bits(crp_plain.z_asymptotic));
  EXPECT_EQ(bits(q.skewing_angle), bits(crp_plain.skewing_angle));
  EXPECT_EQ(q.point_group, crp_plain.point_group);
  EXPECT_NE(q.lattice_constant_x, a_t);
}

// V = CRP(u/s, v/s) + ΔV(x, y)/Constants_eV2Ha bit for bit, with (x, y) the Cartesian position of
// (u, v) and ΔV from a correction built separately from the same files, at points spread over the
// cell, with and without rescaling.
TEST(ScmCrpPotential, ValueIsCrpAtScaledCoordinatesPlusDeltaV) {
  for (const bool rescale : {true, false}) {
    SCOPED_TRACE(rescale);
    const ScmInput input = example_input(rescale);
    const ScmCrpPotential potential(srp48, input);
    const double s = potential.rescaling_factor().value_or(1.0);
    const double cell = potential.system_parameters().lattice_constant_x;
    const ScmCorrection correction(read_scm_parameters(input.parameter_file), potential.slab(),
                                   potential.rescaling_factor(), potential.crp().system_parameters().skewing_angle);
    int compared = 0, nonzero = 0;
    for (int iu = 0; iu < 4; ++iu)
      for (int iv = 0; iv < 4; ++iv)
        for (const double z : {1.0, 2.25, 4.0})
          for (const double r : {1.0, 1.4, 2.2}) {
            const double u = iu * cell / 4.0, v = iv * cell / 4.0, theta = 1.1, phi = 0.7;
            const CartesianVector xy = correction.to_cartesian(u, v);
            const double delta_v = correction.value_hartree(xy.x, xy.y, z, r, theta, phi);
            const double expected =
                potential.crp().value(u / s, v / s, z, r, theta, phi, 0, 0) + units::precise.hartree_to_ev(delta_v);
            EXPECT_EQ(bits(potential.value(u, v, z, r, theta, phi, 0, 0)), bits(expected))
                << "u, v, Z, r = " << u << ", " << v << ", " << z << ", " << r;
            ++compared;
            if (delta_v != 0.0) ++nonzero;
          }
    EXPECT_EQ(compared, 144);
    EXPECT_GT(nonzero, compared / 2);  // the comparison is not one of zeros
  }
}

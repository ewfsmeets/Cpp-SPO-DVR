#include <gtest/gtest.h>

#include <cmath>

#include "core/constants.hpp"
#include "core/machine_epsilon.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"

using namespace spodvr;

TEST(Core, MachineEpsilonIsTwoToMinus54) { EXPECT_EQ(fortran_machine_epsilon(), std::ldexp(1.0, -54)); }

TEST(Core, ConstantsMatchFortran) {
  EXPECT_EQ(Constants::pi, 3.141592653589793);
  EXPECT_EQ(Constants::two_pi, 2.0 * Constants::pi);
  EXPECT_EQ(Constants::nen95_sentinel, -95);
}

TEST(Core, LegacyUnitsDivideByChaev) {
  EXPECT_EQ(units::legacy.ev_to_hartree(1.0), 1.0 / 27.2114);
  EXPECT_EQ(units::legacy.hartree_to_ev(2.0), 2.0 * 27.2114);
  EXPECT_EQ(units::legacy.amu_to_me(1.0078), 1.0078 * 1822.89);
  EXPECT_EQ(units::legacy.angstrom_to_bohr(1.0), 1.0 / 0.529177);
}

TEST(Core, PreciseUnitsMultiplyByEv2Ha) {
  EXPECT_EQ(units::precise.ev_to_hartree(3.0), 3.0 * 0.036749308286);
  EXPECT_EQ(units::precise.hartree_to_ev(3.0), 3.0 / 0.036749308286);
  EXPECT_EQ(units::precise.amu_to_me(1.0), 1822.88853006);
}

TEST(Core, FortranComparisonsUseAbsoluteTolerance) {
  using namespace fortran_compare;
  EXPECT_TRUE(req(1.0, 1.0 + 5e-13));
  EXPECT_FALSE(req(1.0, 1.0 + 2e-12));
  EXPECT_TRUE(rge(1.0, 1.0 + 5e-13));
  EXPECT_FALSE(rgt(1.0 + 5e-13, 1.0));
  EXPECT_TRUE(rgt(1.0 + 2e-12, 1.0));
  EXPECT_TRUE(rlt(1.0, 1.0 + 2e-12));
  EXPECT_FALSE(rlt(1.0, 1.0 + 5e-13));
  EXPECT_TRUE(rle(1.0 + 5e-13, 1.0));
  EXPECT_TRUE(seq("  HardWallPES ", "hardwallpes"));
  EXPECT_FALSE(seq("HardWallPES", "HardWall"));
}

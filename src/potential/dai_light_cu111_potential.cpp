// H2/Cu(111) LEPS potential of Dai and Light. Fortran: PESForCu111DaiLight.f90
// All literals are copied verbatim and the Fortran evaluation order is kept (bitwise results).
#include "potential/dai_light_cu111_potential.hpp"

#include <cmath>

#include "core/constants.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"

namespace spodvr {

namespace {

// Hydrogen masses (amu): the atoms sit at -/+ r/2 around the centre of mass.
constexpr double mass_1 = 1.00794;  // Fortran: U1
constexpr double mass_2 = 1.00794;  // Fortran: U2

// H2 molecular LEPS term (hartree, bohr). Fortran: DAA, ALPAA, R0, DELAA
constexpr double d_molecule = .1744;
constexpr double alpha_molecule = 1.028;
constexpr double r_molecule = 1.40;
constexpr double delta_molecule = -0.15;

// Atom-surface parameters at the top, bridge and centre sites (identical for both atoms).
// Fortran: DASTO1/DASBR1/DASCE1, ALPTO1/..., Z0TO1/..., DELTO1/...
struct SiteValues {
  double top;
  double bridge;
  double centre;
};
constexpr SiteValues d_atom{.638e-1, .818e-1, .877e-1};  // Fortran: .638D-1 etc. (quad literals, rounded once)
constexpr SiteValues alpha_atom{.790, .558, .499};
constexpr SiteValues z_atom{3.02, 2.13, 1.88};
constexpr SiteValues delta_atom{0.10, 0.10, 0.10};

constexpr double angle = 1.04719755119659774614;  // 60 deg, the Fortran literal (not pi/3). Fortran: ANGLE
constexpr double lattice_constant = 4.84;          // bohr. Fortran: CLATIC = xLatticeConstant
constexpr double z_asymptotic = 15.0;              // Fortran: ZAsymptotic

// Coulomb (U) and exchange (A) integrals of a LEPS pair term with Sato parameter delta.
// Fortran: UPOT
double coulomb(double x, double d, double delta, double alpha, double x0) {
  return d / 4.0 / (1.0 + delta) *
         ((3.0 + delta) * std::exp(-2.0 * alpha * (x - x0)) - (2.0 + 6.0 * delta) * std::exp(-alpha * (x - x0)));
}

// Fortran: APOT
double exchange(double x, double d, double delta, double alpha, double x0) {
  return d / 4.0 / (1.0 + delta) *
         ((1.0 + 3.0 * delta) * std::exp(-2.0 * alpha * (x - x0)) - (6.0 + 2.0 * delta) * std::exp(-alpha * (x - x0)));
}

// A site parameter at the skewed position (x, y): the lowest hexagonal Fourier terms through the
// top, bridge and centre values. Fortran: AFUN
double site_interpolation(double x, double y, const SiteValues& v) {
  constexpr double pi2 = 6.28318530717958647688;   // Fortran literal
  constexpr double pi4 = 12.56637061435917295376;  // Fortran literal
  const double c = lattice_constant;
  const double c0 = v.top / 3.0 + v.centre * 2.0 / 3.0;
  const double c1 = v.top / 4.0 - v.bridge / 4.0;
  const double c2 = v.bridge / 4.0 - v.centre * 2.0 / 9.0 - v.top / 36.0;
  const double first_star = std::cos(pi2 / c * x) + std::cos(pi2 / c * y) + std::cos(pi2 / c * (x + y));
  const double second_star = std::cos(pi4 / c * x) + std::cos(pi4 / c * y) + std::cos(pi4 / c * (x + y));
  return c0 + c1 * first_star + c2 * second_star;
}

}  // namespace

// Fortran: module parameters and getParamsForH2Cu111DaiLight
DaiLightCu111Potential::DaiLightCu111Potential()
    : parameters_{lattice_constant, lattice_constant, z_asymptotic, Constants::pi / 3.0, "C6v"} {}

// Fortran: getVH2Cu111DaiLight / VH2Cu111DaiLight. (x, y) are skewed coordinates.
double DaiLightCu111Potential::value(double x, double y, double z, double r, double theta, double phi, int /*iz*/,
                                     int /*ir*/) const {
  const double weight_1 = mass_2 / (mass_1 + mass_2);   // Fortran: UU1
  const double weight_2 = -mass_1 / (mass_1 + mass_2);  // Fortran: UU2

  // Atom positions in skewed coordinates: the bond's projection r sin(theta) decomposed along
  // the two lattice directions.
  const double r_parallel = r * std::sin(theta) / std::sin(angle);  // Fortran: XYT
  double t = std::sin(angle - phi) * r_parallel;
  const double x1 = x + weight_1 * t;
  const double x2 = x + weight_2 * t;
  t = std::sin(phi) * r_parallel;
  const double y1 = y + weight_1 * t;
  const double y2 = y + weight_2 * t;
  t = r * std::cos(theta);
  const double z1 = z + weight_1 * t;
  const double z2 = z + weight_2 * t;

  const double d1 = site_interpolation(x1, y1, d_atom);
  const double d2 = site_interpolation(x2, y2, d_atom);
  const double alpha1 = site_interpolation(x1, y1, alpha_atom);
  const double alpha2 = site_interpolation(x2, y2, alpha_atom);
  const double z01 = site_interpolation(x1, y1, z_atom);
  const double z02 = site_interpolation(x2, y2, z_atom);
  const double delta1 = site_interpolation(x1, y1, delta_atom);
  const double delta2 = site_interpolation(x2, y2, delta_atom);

  const double u_molecule = coulomb(r, d_molecule, delta_molecule, alpha_molecule, r_molecule);
  const double a_molecule = exchange(r, d_molecule, delta_molecule, alpha_molecule, r_molecule);
  const double u_atom1 = coulomb(z1, d1, delta1, alpha1, z01);
  const double a_atom1 = exchange(z1, d1, delta1, alpha1, z01);
  const double u_atom2 = coulomb(z2, d2, delta2, alpha2, z02);
  const double a_atom2 = exchange(z2, d2, delta2, alpha2, z02);

  // LEPS: U_M + U_1 + U_2 - sqrt(A_M^2 + (A_1 + A_2)^2 - A_M (A_1 + A_2)); gas phase beyond zasymp.
  const double a_surface = a_atom1 + a_atom2;
  double v = u_molecule + u_atom1 + u_atom2 -
             std::sqrt(a_molecule * a_molecule + a_surface * a_surface - a_molecule * a_surface);
  if (fortran_compare::rge(z, z_asymptotic)) v = u_molecule + a_molecule;

  return units::legacy.hartree_to_ev(v);
}

}  // namespace spodvr

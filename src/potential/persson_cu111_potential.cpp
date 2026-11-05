// H2/Cu(111) corrugated LEPS potential of Persson et al. Fortran: PESForCu111Persson.f90
// All literals are copied verbatim and the Fortran evaluation order is kept. (The Fortran build
// vectorises some cos/exp calls with libmvec, so values agree with it to a few ulp only.)
#include "potential/persson_cu111_potential.hpp"

#include <array>
#include <cmath>

#include "core/constants.hpp"
#include "core/unit_conversion.hpp"

namespace spodvr {

namespace {

// Fourier coefficients of the atom-surface LEPS parameters (entries 0-3) and the molecular
// parameter (entry 4). Angstrom, eV. Fortran: d, alpha, ze, f0, xa
constexpr std::array<double, 5> d = {2.00654, -7.10400E-02, 2.35238E-02, 8.94111E-03, 4.505};
constexpr std::array<double, 5> alpha = {1.35951, 6.91558E-02, 3.15166E-03, -5.89917E-03, 2.11};
constexpr std::array<double, 5> ze = {1.19994, 7.69362E-02, -1.39699E-02, -1.00913E-02, 0.754};
constexpr std::array<double, 5> f0 = {0.866197, -5.14712E-02, -4.79762E-03, 1.26976E-02, 0.545};
constexpr std::array<double, 5> xa = {1.93833, 1.08333E-02, -8.33338E-04, 2.77771E-04, 1.5};

constexpr double bba = 4.0;    // switching steepness, atom-surface terms
constexpr double bbm = 4.0;    // switching steepness, molecular term
constexpr double dela = 0.01;  // Sato parameter, atom-surface terms
constexpr double delm = -0.11;  // Sato parameter, molecular term
constexpr double cc = 2.5173;  // lattice constant (Angstrom)

constexpr double pi = 3.141592653589793238462643383279502884197;  // Fortran literal (= double pi)
constexpr double bohr_to_angstrom_persson = 0.529167;  // Persson's own factor (defines the PES). Fortran: bohr2angstrom
constexpr double z_asymptotic = 8.0;  // Fortran: ZAsymptotic

struct Position {
  double x;
  double y;
  double z;
};

// Atom-surface parameters at the lateral position (x, y) (Cartesian, Angstrom; x along a lattice
// vector, y towards the next top-site row). Fortran: corrextCu111
struct SurfaceParameters {
  double alpha;  // Fortran: aa
  double depth;  // Fortran: dd
  double z_eq;  // Fortran: zz
  double z_a;  // Fortran: za
  double f;  // Fortran: ff
};

SurfaceParameters corrugation(double x, double y) {
  const double tpc = 2.0 * pi / cc;
  const double tpc3 = tpc / std::sqrt(3.0);
  const double c2x = 2.0 * std::cos(tpc * x);
  const double c4x = 2.0 * std::cos(2.0 * tpc * x);
  const double c2y = 2.0 * std::cos(tpc3 * y);
  const double c4y = 2.0 * std::cos(2.0 * tpc3 * y);
  const double c6y = 2.0 * std::cos(3.0 * tpc3 * y);
  const double c8y = 2.0 * std::cos(4.0 * tpc3 * y);

  // The three Fourier stars of the hexagonal cell.
  const double star1 = c4y + c2x * c2y;
  const double star2 = c4x + c2x * c6y;
  const double star3 = c8y + c4x * c4y;
  auto expand = [&](const std::array<double, 5>& p) { return p[0] + p[1] * star1 + p[2] * star2 + p[3] * star3; };
  return {expand(alpha), expand(d), expand(ze), expand(xa), expand(f0)};
}

// LEPS energy (eV) of the two H atoms at h1, h2 (Angstrom). Fortran: VH2Cu111Persson
double leps(const Position& h1, const Position& h2) {
  const double pa = 1.0 / (4.0 * (1.0 + dela));
  const double u1a = pa * (3.0 + dela);
  const double u2a = pa * (2.0 + 6.0 * dela);
  const double q1a = pa * (1.0 + 3.0 * dela);
  const double q2a = pa * (6.0 + 2.0 * dela);

  const double pm = 1.0 / (4.0 * (1.0 + delm));
  const double u1m = pm * (3.0 + delm);
  const double u2m = pm * (2.0 + 6.0 * delm);
  const double q1m = pm * (1.0 + 3.0 * delm);
  const double q2m = pm * (6.0 + 2.0 * delm);

  // Coulomb and exchange integrals of the atom-surface terms, with the switching xmod.
  std::array<double, 3> coulomb{};
  std::array<double, 3> exchange{};
  const std::array<Position, 2> atoms = {h1, h2};
  for (int i = 0; i < 2; ++i) {
    const SurfaceParameters s = corrugation(atoms[i].x, atoms[i].y);
    const double az = s.alpha * (atoms[i].z - s.z_a);
    const double xmod = std::exp(-az * s.f / (1.0 + std::exp(-bba * az)));
    const double xx = std::exp(-s.alpha * (atoms[i].z - s.z_eq));
    coulomb[i] = s.depth * (u1a * xx - u2a * xmod) * xx;
    exchange[i] = s.depth * (q1a * xx - q2a * xmod) * xx;
  }

  // Molecular term at the H-H distance. Fortran: sum((x(:,2) - x(:,1))**2)
  const double dx = h2.x - h1.x;
  const double dy = h2.y - h1.y;
  const double dz = h2.z - h1.z;
  const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
  const double az = alpha[4] * (r - xa[4]);
  const double xmod = std::exp(-az * f0[4] / (1.0 + std::exp(-bbm * az)));
  const double xx = std::exp(-alpha[4] * (r - ze[4]));
  coulomb[2] = d[4] * (u1m * xx - u2m * xmod) * xx;
  exchange[2] = d[4] * (q1m * xx - q2m * xmod) * xx;

  const double exc12 = exchange[0] + exchange[1];
  return coulomb[0] + coulomb[1] + coulomb[2] -
         std::sqrt(exchange[2] * exchange[2] + exc12 * exc12 - exchange[2] * exc12);
}

}  // namespace

// Fortran: module parameters and getParamsForH2Cu111Persson (lattice constant cc / cboan)
PerssonCu111Potential::PerssonCu111Potential()
    : parameters_{units::legacy.angstrom_to_bohr(cc), units::legacy.angstrom_to_bohr(cc), z_asymptotic,
                  Constants::pi / 3.0, "C6v"} {}

// Fortran: getVH2Cu111Persson. (x, y) are skewed coordinates (bohr), returned V in eV.
double PerssonCu111Potential::value(double x, double y, double z, double r, double theta, double phi, int /*iz*/,
                                    int /*ir*/) const {
  // Bond vector (Cartesian) and centre of mass from skewed to Cartesian coordinates.
  // Fortran: transFromSphericalToCartesian, transFromSkewedToCartesian
  const double x_rel = r * std::sin(theta) * std::cos(phi);
  const double y_rel = r * std::sin(theta) * std::sin(phi);
  const double z_rel = r * std::cos(theta);
  const double x_cm = x + 0.5 * y;
  const double y_cm = 0.5 * std::sqrt(3.0) * y;

  // Atom positions, converted to Angstrom with Persson's factor.
  const double f = bohr_to_angstrom_persson;
  const Position h1{(x_cm - x_rel / 2.0) * f, (y_cm - y_rel / 2.0) * f, (z - z_rel / 2.0) * f};
  const Position h2{(x_cm + x_rel / 2.0) * f, (y_cm + y_rel / 2.0) * f, (z + z_rel / 2.0) * f};
  return leps(h1, h2);
}

}  // namespace spodvr

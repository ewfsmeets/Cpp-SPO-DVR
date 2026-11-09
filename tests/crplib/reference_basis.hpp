#pragma once
// The Fourier bases of the CRP interpolation and the site positions, written down independently of
// the library code (the tests compare the library against them). Sites in lattice units: skewed
// (u, v) for C3v/C6v (x = u + v/2, y = v sqrt(3)/2), Cartesian for C4v.

#include <array>
#include <cmath>
#include <numbers>
#include <vector>

namespace crplib_test {

struct Site {
  double u, v;
};

// The 12 sites of the hexagonal 6D interpolation, in Psite order.
inline const std::array<Site, 12> hexagonal_molecule_sites = {{{0, 0},
                                                               {1. / 2, 0},
                                                               {1. / 3, 1. / 3},
                                                               {1. / 6, 1. / 6},
                                                               {1. / 2, 1. / 2},
                                                               {0, 1. / 2},
                                                               {2. / 3, 1. / 6},
                                                               {1. / 6, 2. / 3},
                                                               {2. / 3, 2. / 3},
                                                               {5. / 6, 5. / 6},
                                                               {1. / 3, 5. / 6},
                                                               {5. / 6, 1. / 3}}};

// The atom-surface sites (skewed), in CRPData.3D order.
inline const std::vector<Site> c6v_atom_sites = {{0, 0},           {0, 1. / 2}, {1. / 3, 1. / 3}, {1. / 6, 1. / 6},
                                                 {1. / 6, 1. / 3}, {0, 1. / 6}, {0, 2. / 3}};
inline const std::vector<Site> c3v_atom_sites = {{0, 0},           {1. / 2, 1. / 2}, {1. / 3, 1. / 3}, {1. / 6, 2. / 3},
                                                 {1. / 6, 1. / 2}, {0, 1. / 6},      {0, 1. / 3},      {1. / 3, 5. / 6},
                                                 {1. / 2, 2. / 3}, {2. / 3, 2. / 3}};

// FOUR6II basis at skewed lattice coordinates (u, v).
inline std::array<double, 12> hexagonal_molecule_basis(double u, double v) {
  const double x = 2 * std::numbers::pi * u, y = 2 * std::numbers::pi * v;
  return {1.0,
          std::cos(x),
          std::cos(y),
          std::cos(x + y),
          std::sin(x),
          std::sin(y),
          std::sin(x + y),
          std::cos(2 * x + y),
          std::cos(x + 2 * y),
          std::cos(x - y),
          std::cos(2 * x) + std::cos(2 * y) + std::cos(2 * (x + y)),
          -std::sin(2 * x) - std::sin(2 * y) + std::sin(2 * (x + y))};
}

// FOUR3II basis at skewed lattice coordinates (u, v): the 7 cosine stars (C6v), with 3 sine
// stars inserted (C3v, 10 functions).
inline std::vector<double> hexagonal_atom_basis(double u, double v, bool c3v) {
  const double x = 2 * std::numbers::pi * u, y = 2 * std::numbers::pi * v;
  auto c = [](double a) { return std::cos(a); };
  auto s = [](double a) { return std::sin(a); };
  const double star1 = c(x) + c(y) + c(x + y);
  const double sine1 = -s(x) - s(y) + s(x + y);
  const double star2 = c(2 * x + y) + c(x + 2 * y) + c(x - y);
  const double star3 = c(2 * x) + c(2 * y) + c(2 * (x + y));
  const double sine3 = -s(2 * x) - s(2 * y) + s(2 * (x + y));
  const double star4 = c(3 * x + y) + c(x + 3 * y) + c(2 * x + 3 * y) + c(3 * x + 2 * y) + c(-x + 2 * y) + c(2 * x - y);
  const double sine4 =
      -s(3 * x + y) - s(x + 3 * y) + s(2 * x + 3 * y) + s(3 * x + 2 * y) - s(-x + 2 * y) - s(2 * x - y);
  const double star5 = c(3 * x) + c(3 * y) + c(3 * (x + y));
  const double star6 = c(2 * (2 * x + y)) + c(2 * (x + 2 * y)) + c(2 * (x - y));
  if (!c3v) return {1.0, star1, star2, star3, star4, star5, star6};
  return {1.0, star1, sine1, star2, star3, sine3, star4, sine4, star5, star6};
}

// Skewed lattice coordinates of a Cartesian point (bohr) for lattice constant a.
inline Site skewed(double x, double y, double a) {
  const double v = y / (a * std::sqrt(3.0) / 2.0);
  return {x / a - v / 2.0, v};
}

}  // namespace crplib_test

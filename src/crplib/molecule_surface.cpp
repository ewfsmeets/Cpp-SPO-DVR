// The molecule-surface PES: reading the molecular data into splines, the molecule geometry, the
// near-surface interpolation V = Va + Vb + I, the switch to the 2D gas-phase potential and the
// chain rule to atom gradients. Fortran: CRP_C3v/C4v/C6v_MoleculeSurface.f90 (Read_MoleculeSurface_
// Data, POT6D, RRINTII, RRINT), CRP_GasPhase.f90, CRP_General.f90 CAMBCORD. The C3v/C6v and C4v
// variants differ only in the site and Fourier step (hexagonal_ / square_interaction.cpp) and in
// the operation order of two chain-rule terms.
#include "crplib/molecule_surface.hpp"

#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <string>

#include "crplib/constants.hpp"
#include "crplib/errors.hpp"
#include "crplib/molecule_surface_detail.hpp"

namespace crplib {

using detail::BicubicValue;

namespace detail {

// Fortran: CAMBCORD, with ccdef = 0.5 (all symmetries)
MoleculeGeometry molecule_geometry(const Vec3& a, const Vec3& b) {
  constexpr double ccdef = 0.500;
  const double r = std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z - a.z) * (b.z - a.z));
  if (r == 0.0) throw EvaluationError("CRP molecule-surface PES: the two atoms coincide (Fortran CAMBCORD: r = 0)");
  MoleculeGeometry g{};
  g.r = r;
  g.x = a.x * (1.0 - ccdef) + b.x * ccdef;
  g.y = a.y * (1.0 - ccdef) + b.y * ccdef;
  g.z = a.z * (1.0 - ccdef) + b.z * ccdef;
  double cos_theta = (b.z - a.z) / r;
  if (cos_theta > 1.0) cos_theta = 1.0;  // clipped at +1 only, as the Fortran (no SPO-DVR grid reaches -1)
  g.theta = std::acos(cos_theta);
  const double dx = b.x - a.x;
  const double dy = b.y - a.y;
  g.phi = (dx == 0.0 && dy == 0.0) ? 0.0 : std::atan2(dy, dx);
  return g;
}

}  // namespace detail

namespace {

// Derivatives of r, theta and phi with respect to the coordinates of atom a (those of atom b are
// their negatives). Fortran: drXa ... dPhiYa of RRINTII / RRINT
struct GeometryDerivatives {
  double r_x, r_y, r_z;
  double theta_x, theta_y, theta_z;
  double phi_x, phi_y;
};

// Fortran: RRINTII (C3v, C6v)
GeometryDerivatives hexagonal_chain_rule(const Vec3& a, const Vec3& b, double r) {
  GeometryDerivatives d{};
  d.r_x = -(b.x - a.x) / r;
  d.r_y = -(b.y - a.y) / r;
  d.r_z = -(b.z - a.z) / r;
  const double rho2 = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
  const double rho = std::sqrt(rho2);
  if (rho2 > 0.00) {
    d.theta_x = -(b.z - a.z) * (b.x - a.x) / (r * r * rho);
    d.theta_y = -(b.z - a.z) * (b.y - a.y) / (r * r * rho);
    d.phi_x = (b.y - a.y) / rho2;
    d.phi_y = -(b.x - a.x) / rho2;
  }
  d.theta_z = rho / (r * r);
  return d;
}

// Fortran: RRINT (C4v); the same derivatives with another order of operations (r**2.0, SQRT inside).
GeometryDerivatives square_chain_rule(const Vec3& a, const Vec3& b, double r) {
  GeometryDerivatives d{};
  d.r_x = -(b.x - a.x) / r;
  d.r_y = -(b.y - a.y) / r;
  d.r_z = -(b.z - a.z) / r;
  const double rho2 = (b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y);
  if (rho2 > 0.0) {
    d.theta_x = -(b.z - a.z) * (b.x - a.x) / (r * r) / std::sqrt(rho2);
    d.theta_y = -(b.z - a.z) * (b.y - a.y) / (r * r) / std::sqrt(rho2);
    d.phi_x = (b.y - a.y) / rho2;
    d.phi_y = -(b.x - a.x) / rho2;
  }
  d.theta_z = std::sqrt(rho2) / (r * r);
  return d;
}

std::vector<double> converted_lengths(const std::vector<double>& x, const Units& from) {
  std::vector<double> y;
  for (double v : x) y.push_back(convert_length(v, from.length, LengthUnit::bohr));
  return y;
}
std::vector<double> converted_energies(const std::vector<double>& x, const Units& from) {
  std::vector<double> y;
  for (double v : x) y.push_back(convert_energy(v, from.energy, EnergyUnit::hartree));
  return y;
}

// The centre of a molecular cut's site in Cartesian bohr, as Read_MoleculeSurface_Data computes it
// for the reduction: skewed lattice units for the hexagonal symmetries, Cartesian for C4v.
std::pair<double, double> cut_site(Symmetry symmetry, const MoleculeCut& cut, double a) {
  if (symmetry == Symmetry::c4v) return {cut.x * a, cut.y * a};
  const double ucm = cut.x * a;
  const double vcm = cut.y * a;
  return {ucm + vcm / 2.0, std::sqrt(3.0) * vcm / 2.0};
}

// Bitwise equal grids (S4: equal comparisons alone would merge -0 and +0 knots).
bool same_grid(const std::vector<double>& a, const std::vector<double>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
  return true;
}

}  // namespace

// Fortran: Initialize_CRP_*_MoleculeSurface_PES and Read_MoleculeSurface_Data
MoleculeSurface::MoleculeSurface(const DataSet& data, EvaluationOptions options)
    : symmetry_(data.symmetry),
      lattice_constant_(convert_length(data.lattice_constant, data.units.length, LengthUnit::bohr)),
      share_grid_lookup_(options.share_grid_lookup) {
  validate_data_set(data);
  atom_ = make_atom_surface(data, options);
  const Units& in = data.units;
  zso_min_ = convert_length(data.zso_min, in.length, LengthUnit::bohr);
  zso_max_ = convert_length(data.zso_max, in.length, LengthUnit::bohr);
  // Zdep: not-a-knot at the start, zero slope at the end. Fortran: Compute_Spline1D(Zdep, 0, 0.00, 2, 0.00)
  z_dependence_ =
      detail::CubicSpline1D(converted_lengths(data.z_dependence.x, in), converted_energies(data.z_dependence.f, in),
                            detail::not_a_knot, {detail::EndCondition::first_derivative, 0.00});
  // Fortran: Compute_Spline1D(Rdep) (not-a-knot at both ends)
  h2_ = detail::CubicSpline1D(converted_lengths(data.h2_potential.x, in), converted_energies(data.h2_potential.f, in),
                              detail::not_a_knot, detail::not_a_knot);

  // The range common to all cuts (also those beyond the ones used). The Fortran starts from SAVEd
  // locals Rmin = 0, Rmax = 100, Zmin = -100.
  double r_min = 0.0, r_max = 100.0, z_min = -100.0;
  for (const MoleculeCut& cut : data.molecule_cuts) {
    const std::vector<double> r = converted_lengths(cut.r, in);
    const std::vector<double> z = converted_lengths(cut.z, in);
    std::vector<double> v = converted_energies(cut.v, in);
    if (!data.cuts_interpolation) {
      // The data are total energies: subtract both atom-surface potentials (corrugation reduction).
      // The Fortran evaluates cos(theta) and cos(phi) here with libmvec (a few ulp different).
      const auto [xcm, ycm] = cut_site(symmetry_, cut, lattice_constant_);
      const double theta = cut.theta * detail::pi / 180.0;
      const double phi = cut.phi * detail::pi / 180.0;
      for (std::size_t j = 0; j < z.size(); ++j) {
        for (std::size_t i = 0; i < r.size(); ++i) {
          const Vec3 a{xcm - 0.5 * r[i] * std::sin(theta) * std::cos(phi),
                       ycm - 0.5 * r[i] * std::sin(theta) * std::sin(phi), z[j] - 0.5 * r[i] * std::cos(theta)};
          const Vec3 b{xcm + 0.5 * r[i] * std::sin(theta) * std::cos(phi),
                       ycm + 0.5 * r[i] * std::sin(theta) * std::sin(phi), z[j] + 0.5 * r[i] * std::cos(theta)};
          double va = 0.0, vb = 0.0;
          try {
            va = atom_->evaluate(a).v;
            vb = atom_->evaluate(b).v;
          } catch (const EvaluationError& e) {
            throw DataError("Molecule/" + cut.file_name + ": " + e.what() +
                            " (CutsInterpolation is F: the atom-surface PES is needed at every cut point)");
          }
          double& value = v[i + r.size() * j];
          value = value - va - vb;
        }
      }
    }
    if (r_min < r.front()) r_min = r.front();
    if (r_max > r.back()) r_max = r.back();
    if (z_min < z.front()) z_min = z.front();
    cuts_.emplace_back(r, z, v);
  }
  r_min_ = r_min;
  r_max_ = r_max;
  z_min_ = z_min;
  // S4: group the used cuts by their grid.
  for (int k = 0; k < molecule_cut_count(symmetry_); ++k) {
    int group = -1;
    for (std::size_t g = 0; g < grid_cut_.size() && group < 0; ++g) {
      const detail::BicubicSpline& first = cuts_[grid_cut_[g]];
      if (same_grid(first.x(), cuts_[k].x()) && same_grid(first.y(), cuts_[k].y())) group = static_cast<int>(g);
    }
    if (group < 0) {
      group = static_cast<int>(grid_cut_.size());
      grid_cut_.push_back(k);
    }
    group_of_cut_.push_back(group);
  }
  warnings_ = data_set_warnings(data);
  for (const std::string& w : atom_->warnings()) warnings_.push_back(w);
}

std::vector<std::string> MoleculeSurface::warnings() const { return warnings_; }

// Fortran: Switch_Function
MoleculeSurface::Switch MoleculeSurface::switch_function(double z) const {
  if (z <= zso_min_) return {1.00, 0.00};
  if (z >= zso_max_) return {0.00, 0.00};
  const double f = 0.500 + 0.500 * std::cos((z - zso_min_) * detail::pi / (zso_max_ - zso_min_));
  const double df =
      -0.500 * detail::pi / (zso_max_ - zso_min_) * std::sin((z - zso_min_) * detail::pi / (zso_max_ - zso_min_));
  return {f, df};
}

// Fortran: Calculate_CRP_GasPhase_Potential
MoleculeSurface::GasPhaseValue MoleculeSurface::gas_phase(double z, double r) const {
  const auto h2 = h2_.value_and_slope_splint(r);
  double v1d = 0.0, d_z = 0.0;
  if (z >= zso_min_ && z <= z_dependence_.x().back()) {
    const auto zdep = z_dependence_.value_and_slope_splint(z);
    v1d = zdep.value;
    d_z = zdep.slope;
  } else if (z > z_dependence_.x().back()) {
    v1d = z_dependence_.y().back();
    d_z = 0.00;
  } else {
    throw EvaluationError("CRP gas-phase potential: Z = " + std::to_string(z) +
                          " below ZsoMin (Fortran: ERROR IN V2DINF)");
  }
  return {h2.value + v1d, d_z, h2.slope};
}

bool MoleculeSurface::covers(double r, double z) const {
  if (z >= zso_max_) return !(r < h2_.x().front() || r > h2_.x().back());
  const int n = molecule_cut_count(symmetry_);
  for (int k = 0; k < n; ++k)
    if (!cuts_[k].locate(r, z)) return false;
  return true;
}

// The cuts' values and derivatives at (r, Z); false if (r, Z) is outside a cut's grid (the
// Fortran's switch, whichever cut it is). S4: the cell is located once per group of cuts with the
// same grid (the same cell and weights, so the same values).
bool MoleculeSurface::evaluate_cuts(double r, double z, std::span<BicubicValue> values) const {
  const int n = static_cast<int>(values.size());
  if (!share_grid_lookup_) {
    for (int k = 0; k < n; ++k) {
      const std::optional<BicubicValue> v = cuts_[k].evaluate(r, z);
      if (!v) return false;
      values[k] = *v;
    }
    return true;
  }
  std::array<detail::GridCell, 29> cells{};
  for (std::size_t g = 0; g < grid_cut_.size(); ++g) {
    const std::optional<detail::GridCell> cell = cuts_[grid_cut_[g]].locate(r, z);
    if (!cell) return false;
    cells[g] = *cell;
  }
  for (int k = 0; k < n; ++k) values[k] = cuts_[k].evaluate_in_cell(cells[group_of_cut_[k]]);
  return true;
}

// Fortran: POT6D -> CAMBCORD -> RRINTII (C3v, C6v) / RRINT (C4v)
std::optional<MoleculeSurfaceValue> MoleculeSurface::evaluate(const Vec3& a, const Vec3& b) const {
  const detail::MoleculeGeometry g = detail::molecule_geometry(a, b);

  if (!(g.z < zso_max_)) {
    // Far from the surface: the 2D gas-phase potential only
    if (g.r < h2_.x().front() || g.r > h2_.x().back()) return std::nullopt;  // Fortran: switch = .TRUE.
    const GasPhaseValue gp = gas_phase(g.z, g.r);
    const double r_xa = -(b.x - a.x) / g.r, r_xb = -r_xa;
    const double r_ya = -(b.y - a.y) / g.r, r_yb = -r_ya;
    const double r_za = -(b.z - a.z) / g.r, r_zb = -r_za;
    return MoleculeSurfaceValue{gp.v,
                                {r_xa * gp.d_r, r_ya * gp.d_r, 0.500 * gp.d_z + r_za * gp.d_r},
                                {r_xb * gp.d_r, r_yb * gp.d_r, 0.500 * gp.d_z + r_zb * gp.d_r}};
  }

  // The cuts at the molecule's (r, Z); outside any cut's grid the Fortran gives up (switch).
  const int n = molecule_cut_count(symmetry_);
  std::array<BicubicValue, 29> values{};
  if (!evaluate_cuts(g.r, g.z, std::span<BicubicValue>(values).first(n))) return std::nullopt;
  const std::span<const BicubicValue> cut_values = std::span<const BicubicValue>(values).first(n);
  const bool square = symmetry_ == Symmetry::c4v;
  const detail::Interaction vi =
      square ? detail::square_interaction(cut_values, g.theta, g.phi, g.x, g.y, lattice_constant_)
             : detail::hexagonal_interaction(symmetry_, cut_values, g.theta, g.phi, g.x, g.y, lattice_constant_);
  const GeometryDerivatives d = square ? square_chain_rule(a, b, g.r) : hexagonal_chain_rule(a, b, g.r);
  const AtomSurfaceValue va = atom_->evaluate(a);
  const AtomSurfaceValue vb = atom_->evaluate(b);

  // dV/d(atom coordinate) of Va + Vb + I: the centre moves by half the atom displacement (0.5 *
  // dI/dX, dI/dY, dI/dZ), r, theta and phi by the derivatives above (those of atom b negated).
  const double i_v = vi[0], i_z = vi[1], i_r = vi[2], i_phi = vi[3], i_theta = vi[4], i_x = vi[5], i_y = vi[6];
  const Vec3 full_a{va.gradient.x + 0.500 * i_x + d.r_x * i_r + d.theta_x * i_theta + d.phi_x * i_phi,
                    va.gradient.y + 0.500 * i_y + d.r_y * i_r + d.theta_y * i_theta + d.phi_y * i_phi,
                    va.gradient.z + 0.500 * i_z + d.r_z * i_r + d.theta_z * i_theta};
  const Vec3 full_b{vb.gradient.x + 0.500 * i_x + (-d.r_x) * i_r + (-d.theta_x) * i_theta + (-d.phi_x) * i_phi,
                    vb.gradient.y + 0.500 * i_y + (-d.r_y) * i_r + (-d.theta_y) * i_theta + (-d.phi_y) * i_phi,
                    vb.gradient.z + 0.500 * i_z + (-d.r_z) * i_r + (-d.theta_z) * i_theta};
  if (g.z <= zso_min_) return MoleculeSurfaceValue{va.v + vb.v + i_v, full_a, full_b};

  // Switching region: V = (Va + Vb + I) F + (1 - F) V2D
  const Switch sw = switch_function(g.z);
  const GasPhaseValue gp = gas_phase(g.z, g.r);
  const double df_a = 0.500 * sw.df_dz;  // Fortran: dFsodZa (= dFsodZb)
  const double df_b = 0.500 * sw.df_dz;
  const double v = (va.v + vb.v + i_v) * sw.f + (1.00 - sw.f) * gp.v;
  const Vec3 grad_a{
      full_a.x * sw.f + (1.00 - sw.f) * gp.d_r * d.r_x, full_a.y * sw.f + (1.00 - sw.f) * gp.d_r * d.r_y,
      full_a.z * sw.f + (va.v + vb.v + i_v) * df_a - df_a * gp.v + (1.00 - sw.f) * (gp.d_z * 0.500 + d.r_z * gp.d_r)};
  const Vec3 grad_b{full_b.x * sw.f + (1.00 - sw.f) * gp.d_r * (-d.r_x),
                    full_b.y * sw.f + (1.00 - sw.f) * gp.d_r * (-d.r_y),
                    full_b.z * sw.f + (va.v + vb.v + i_v) * df_b - df_b * gp.v +
                        (1.00 - sw.f) * (gp.d_z * 0.500 + (-d.r_z) * gp.d_r)};
  return MoleculeSurfaceValue{v, grad_a, grad_b};
}

}  // namespace crplib

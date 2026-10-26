// Analysis basis: diffraction channels, analysed rotational states, rotation-diffraction states
// and their channel energies. Fortran: RotDifBasis.f90, GridEnergies.f90 (edif, enymo)
#include "analysis/analysis_basis.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"

namespace spodvr {

int diffraction_order(const std::string& point_group, int n, int m) {
  const int abs_n = std::abs(n);
  const int abs_m = std::abs(m);
  if (fortran_compare::seq(point_group, "C2v") || fortran_compare::seq(point_group, "C4v")) return abs_n + abs_m;
  if (fortran_compare::seq(point_group, "C3") || fortran_compare::seq(point_group, "C3v") || fortran_compare::seq(point_group, "C6v")) {
    if (n * m > 0) return std::max(abs_n, abs_m);
    return abs_n + abs_m;
  }
  throw NotSupportedError("diffraction order: point group " + point_group + " not implemented");
}

namespace {

// Quantum numbers -(n/2 - 1) .. n/2 along one lattice direction (only 0 without that DoF).
struct QuantumRange {
  int min;
  int max;
};

QuantumRange quantum_range(bool dof, int n_points) {
  if (!dof) return {0, 0};
  return {-(n_points / 2 - 1), n_points / 2};
}

// FFT slot of quantum number q on an n-point grid: q >= 0 at q, q < 0 at n + q.
int fft_slot(int q, int n_points) { return q >= 0 ? q : n_points + q; }

}  // namespace

int max_open_diffraction_order(const Config& config, const GasPhase& gas, int mx, int my) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  const QuantumRange qx = quantum_range(p.x_dof, mx);
  const QuantumRange qy = quantum_range(p.y_dof, my);
  const double mass = d.total_mass;
  const double cos_skew = std::cos(d.system.skewing_angle);

  // Parallel kinetic energy of channel (i, j), evaluated as in Fortran.
  auto parallel_energy = [&](int i, int j) {
    const double kx = d.kx + i * d.delta_kx;
    const double ky = d.ky + j * d.delta_ky;
    const double ekinx = 0.5 * (kx * kx) / mass;
    const double ekiny = 0.5 * (ky * ky) / mass;
    const double ekxky = -cos_skew * (kx * ky) / mass;
    return ekinx + ekxky + ekiny;
  };

  // Largest energy available for the parallel motion: EMaxWP + E(vin, jin) - E(0, 0) + E_par.
  const double e00 = gas.states.energy(0, 0);
  const double e_initial = gas.states.energy(p.v_in, p.j_in);
  const double e_total = units::legacy.ev_to_hartree(p.e_max_wave_packet) + (e_initial - e00) + parallel_energy(0, 0);

  int max_order = 0;
  for (int j = qy.min; j <= qy.max; ++j)
    for (int i = qx.min; i <= qx.max; ++i)
      if (parallel_energy(i, j) <= e_total)
        max_order = std::max(max_order, diffraction_order(d.system.point_group, i, j));
  return max_order;
}

AnalysisBasis::AnalysisBasis(const Config& config, const AngularGrid& angular, const RotationalBasis& basis,
                             const GasPhase& gas, const std::vector<double>& xy_kinetic_energies)
    : basis_(basis), j_max_(gas.j_max_analysis), n_v_(gas.n_v_analysis) {
  const InputParameters& p = config.input;
  const std::string& point_group = config.derived.system.point_group;

  // Diffraction channels n outer, m inner, then sorted by order (stable). Fortran: DiffractionBasis
  const QuantumRange qx = quantum_range(p.x_dof, p.nx);
  const QuantumRange qy = quantum_range(p.y_dof, p.ny);
  for (int n = qx.min; n <= qx.max; ++n)
    for (int m = qy.min; m <= qy.max; ++m)
      channels_.push_back({n, m, diffraction_order(point_group, n, m), fft_slot(n, p.nx), fft_slot(m, p.ny)});
  std::stable_sort(channels_.begin(), channels_.end(),
                   [](const DiffractionChannel& a, const DiffractionChannel& b) { return a.order < b.order; });

  highest_order_ = channels_.back().order;
  for (int order = 0; order <= highest_order_; ++order) {
    const auto count = std::count_if(channels_.begin(), channels_.end(),
                                     [order](const DiffractionChannel& c) { return c.order == order; });
    max_per_order_ = std::max(max_per_order_, static_cast<int>(count));
  }

  // Analysed orders: DifOrdAna, or (absent) all energetically open orders.
  max_order_ = p.diffraction_order_analysis ? *p.diffraction_order_analysis
                                            : max_open_diffraction_order(config, gas, p.nx, p.ny);
  too_small_ = max_open_diffraction_order(config, gas, p.nx + 2, p.ny + 2) > max_order_;
  for (const DiffractionChannel& c : channels_)
    if (c.order <= max_order_) n_diffraction_ = n_diffraction_ + 1;

  for (const DiffractionChannel& c : channels_)
    diffraction_energies_.push_back(xy_kinetic_energies[c.iy + p.ny * c.ix]);

  // Analysed rotational states: those of the propagation basis with J <= jana, in basis order.
  // Fortran: RotationBasis
  for (int irot = 0; irot < basis.size(); ++irot)
    if (basis.j(irot) <= j_max_) rotations_.push_back(irot);

  // The initial state (JIn, MJIn) in the specular diffraction channel (the first one).
  // Fortran: RotationDiffractionBasis (inchna)
  for (int irotana = 0; irotana < n_rotation(); ++irotana) {
    const int irot = rotations_[irotana];
    if (basis.j(irot) == p.j_in && basis.mj(irot) == p.mj_in) initial_ = irotana * n_diffraction_;
  }

  // Channel energies E_dif + E(v, j). Fortran: SetupGridEnergies (enymo)
  channel_energies_.resize(n_states() * n_v_);
  for (int v = 0; v < n_v_; ++v)
    for (int is = 0; is < n_states(); ++is)
      channel_energies_[is + n_states() * v] = diffraction_energy(state_diffraction(is)) + gas.states.energy(v, j(is));

  // Regular order: J = JMin .. jana (step dJ), mJ = -MjMax .. MjMax. Fortran: mapToScatteringRotationalBasis
  for (int jj = angular.j_min(); jj <= j_max_; jj += angular.delta_j())
    for (int mm = -basis.mj_max(); mm <= basis.mj_max(); ++mm) {
      if (std::abs(mm) > jj) continue;
      for (int irotana = 0; irotana < n_rotation(); ++irotana)
        if (basis.j(rotations_[irotana]) == jj && basis.mj(rotations_[irotana]) == mm)
          regular_rotations_.push_back(irotana);
    }
}

int AnalysisBasis::j(int is) const { return basis_.j(rotations_[state_rotation(is)]); }

int AnalysisBasis::mj(int is) const { return basis_.mj(rotations_[state_rotation(is)]); }

std::vector<int> AnalysisBasis::regular_state_order() const {
  std::vector<int> order;
  for (const int irotana : regular_rotations_)
    for (int idif = 0; idif < n_diffraction_; ++idif) order.push_back(irotana * n_diffraction_ + idif);
  return order;
}

}  // namespace spodvr

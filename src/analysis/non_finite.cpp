// The first non-finite number (NaN or inf) of the scattering results, for the run's exit code 5.
// C++ only: the Fortran has no such check.
#include "analysis/non_finite.hpp"

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace spodvr {

namespace {

bool is_finite(double x) { return std::isfinite(x); }
bool is_finite(std::complex<double> x) { return std::isfinite(x.real()) && std::isfinite(x.imag()); }

// The index of the first non-finite element of `values`, nothing if all are finite.
template <typename T>
std::optional<int> first_non_finite(const std::vector<T>& values) {
  for (std::size_t i = 0; i < values.size(); ++i)
    if (!is_finite(values[i])) return static_cast<int>(i);
  return std::nullopt;
}

// The first non-finite number of one energy's results, in ResultQuantity order.
std::optional<NonFiniteValue> find_in_energy(const EnergyResult& e, int energy_index) {
  struct Scalar {
    ResultQuantity quantity;
    double value;
  };
  struct Vector {
    ResultQuantity quantity;
    const std::vector<double>& values;
  };
  using enum ResultQuantity;

  for (const Scalar& s : {Scalar{total_energy, e.total_energy}, Scalar{collision_energy, e.collision_energy},
                          Scalar{z_energy, e.z_energy}, Scalar{total_momentum, e.total_momentum}})
    if (!is_finite(s.value)) return NonFiniteValue{energy_index, s.quantity, 0};
  if (const std::optional<int> i = first_non_finite(e.s_matrix)) return NonFiniteValue{energy_index, s_matrix, *i};
  for (const Vector& v :
       {Vector{probability, e.probability}, Vector{phase, e.phase}, Vector{p_vibration, e.p_vibration},
        Vector{p_rotation, e.p_rotation}, Vector{p_order, e.p_order}, Vector{e_rotation_v, e.e_rotation_v},
        Vector{e_diffraction_v, e.e_diffraction_v}})
    if (const std::optional<int> i = first_non_finite(v.values)) return NonFiniteValue{energy_index, v.quantity, *i};
  for (const Scalar& s :
       {Scalar{p_total, e.p_total}, Scalar{p_reaction, e.p_reaction}, Scalar{e_vibration, e.e_vibration},
        Scalar{e_rotation, e.e_rotation}, Scalar{e_diffraction, e.e_diffraction}})
    if (!is_finite(s.value)) return NonFiniteValue{energy_index, s.quantity, 0};
  return std::nullopt;
}

}  // namespace

std::optional<NonFiniteValue> find_non_finite(const ScatteringResults& results) {
  for (std::size_t ien = 0; ien < results.energies.size(); ++ien)
    if (const std::optional<NonFiniteValue> found = find_in_energy(results.energies[ien], static_cast<int>(ien)))
      return found;
  return std::nullopt;
}

}  // namespace spodvr

#pragma once

#include <optional>

#include "analysis/s_matrix.hpp"

namespace spodvr {

// The numbers of an EnergyResult, in the order find_non_finite scans them (the order of their
// declaration in EnergyResult; `open` holds no floating-point number and is not scanned).
enum class ResultQuantity {
  total_energy,
  collision_energy,
  z_energy,
  total_momentum,
  s_matrix,         // (is, v)
  probability,      // (is, v)
  phase,            // (is, v)
  p_vibration,      // (v)
  p_rotation,       // (irotana, v)
  p_order,          // (order, irotana, v)
  e_rotation_v,     // (v)
  e_diffraction_v,  // (v)
  p_total,
  p_reaction,
  e_vibration,
  e_rotation,
  e_diffraction,
};

// Where a non-finite number (NaN or +-inf) sits in a ScatteringResults: the energy, the quantity
// and the flat index into that quantity's vector (layouts in EnergyResult; 0 for a scalar).
struct NonFiniteValue {
  int energy_index;
  ResultQuantity quantity;
  int index;
};

// The first NaN or infinite number of `results`: energies in order, within an energy the
// quantities in ResultQuantity order, within a vector ascending index (a complex S counts when
// either part is not finite). Only finiteness is checked: P > 1, a negative reaction
// probability, S = 0 and empty lists are legitimate results. Nothing if every number is finite.
// C++ only (the Fortran writes such results without notice).
[[nodiscard]] std::optional<NonFiniteValue> find_non_finite(const ScatteringResults& results);

}  // namespace spodvr

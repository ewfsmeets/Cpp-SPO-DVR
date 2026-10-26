#pragma once

#include <string>
#include <vector>

#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "input/config.hpp"

namespace spodvr {

// Diffraction order of the channel (n, m) for the surface point group: |n| + |m| for C2v and
// C4v; for C3, C3v and C6v max(|n|, |m|) if n m > 0, else |n| + |m|.
// Fortran: RotDifBasis.f90 getDiffractionOrder
[[nodiscard]] int diffraction_order(const std::string& point_group, int n, int m);

// One diffraction channel (n, m): the reciprocal lattice vector n b1 + m b2 of the parallel
// momentum, its order and its slot (ix, iy) in the (x, y) FFT output. 0-based slots.
struct DiffractionChannel {
  int n;
  int m;
  int order;
  int ix;  // Fortran: PntToNQ - 1
  int iy;  // Fortran: PntToMQ - 1
};

// Highest diffraction order that is energetically open for the largest energy in the wave
// packet (EMaxWP) on an mx x my momentum grid. Fortran: RotDifBasis.f90 getMaxDiffOrderAnalysis
[[nodiscard]] int max_open_diffraction_order(const Config& config, const GasPhase& gas, int mx, int my);

// The states in which the scattered wave packet is analysed: all diffraction channels of order
// <= mdifa combined with all rotational states of the propagation basis with J <= jana, for the
// vibrational levels v < nviba. A rotation-diffraction state index is
// is = irotana * ndifa + idif (rotation outer, diffraction inner, as in Fortran). 0-based.
// Fortran: RotDifBasis.f90 (DiffractionBasis, RotationBasis, RotationDiffractionBasis),
// GridEnergies.f90 (edif, enymo)
class AnalysisBasis {
 public:
  // xy_kinetic_energies: ekinxy at ixy = iy + ny*ix (PropagationTables::xy_kinetic_energies).
  AnalysisBasis(const Config& config, const AngularGrid& angular, const RotationalBasis& basis, const GasPhase& gas,
                const std::vector<double>& xy_kinetic_energies);

  // --- diffraction --------------------------------------------------------------------------
  // All nx*ny channels, sorted by increasing order (stable). Fortran: DifBas, ndftow
  [[nodiscard]] const std::vector<DiffractionChannel>& diffraction_channels() const { return channels_; }
  [[nodiscard]] int n_diffraction() const { return n_diffraction_; }   // Fortran: ndifa
  [[nodiscard]] int max_order() const { return max_order_; }           // Fortran: mdifa
  [[nodiscard]] int highest_order() const { return highest_order_; }   // Fortran: maxord
  [[nodiscard]] int max_channels_per_order() const { return max_per_order_; }  // Fortran: maxnum
  // True if the order open on an (nx+2, ny+2) grid exceeds mdifa (Fortran prints a warning).
  [[nodiscard]] bool diffraction_basis_too_small() const { return too_small_; }

  // --- rotation -----------------------------------------------------------------------------
  [[nodiscard]] int n_rotation() const { return static_cast<int>(rotations_.size()); }  // Fortran: nana
  // Index in the propagation basis of analysed rotational state irotana. Fortran: RBIndex
  [[nodiscard]] int rotation_index(int irotana) const { return rotations_[irotana]; }
  [[nodiscard]] int j_max() const { return j_max_; }   // Fortran: jana

  // --- rotation-diffraction states ----------------------------------------------------------
  [[nodiscard]] int n_states() const { return n_rotation() * n_diffraction_; }  // Fortran: nstata
  [[nodiscard]] int n_v() const { return n_v_; }                                // Fortran: nviba
  [[nodiscard]] int state_rotation(int is) const { return is / n_diffraction_; }   // Fortran: PntToRot - 1
  [[nodiscard]] int state_diffraction(int is) const { return is % n_diffraction_; }  // Fortran: PntToDif - 1
  [[nodiscard]] int j(int is) const;
  [[nodiscard]] int mj(int is) const;
  [[nodiscard]] const DiffractionChannel& channel(int is) const { return channels_[state_diffraction(is)]; }
  [[nodiscard]] int initial_state() const { return initial_; }  // (JIn, MJIn, n = m = 0). Fortran: inchna - 1

  // Parallel kinetic energy of diffraction channel idif (hartree). Fortran: edif
  [[nodiscard]] double diffraction_energy(int idif) const { return diffraction_energies_[idif]; }
  // Channel energy E_dif + E(v, j) of state (is, v) (hartree). Fortran: enymo(is, v)
  [[nodiscard]] double channel_energy(int is, int v) const { return channel_energies_[is + n_states() * v]; }

  // Analysed rotational states (irotana) in the "regular" order J ascending, then
  // mJ = -MjMax .. MjMax, used for printing. Fortran: mapToScatteringRotationalBasis
  [[nodiscard]] const std::vector<int>& regular_rotation_order() const { return regular_rotations_; }
  // All states in regular order: rotations as above, diffraction channels inner.
  // Fortran: mapToScatteringRotDifBasis
  [[nodiscard]] std::vector<int> regular_state_order() const;

 private:
  const RotationalBasis& basis_;
  std::vector<DiffractionChannel> channels_;
  int n_diffraction_ = 0;
  int max_order_ = 0;
  int highest_order_ = 0;
  int max_per_order_ = 0;
  bool too_small_ = false;
  std::vector<int> rotations_;
  int j_max_ = 0;
  int n_v_ = 0;
  int initial_ = -1;
  std::vector<double> diffraction_energies_;
  std::vector<double> channel_energies_;
  std::vector<int> regular_rotations_;
};

}  // namespace spodvr

#pragma once
// The optional datasets of a frame (FrameExtras, C++ only): reductions of psi_full that need a transform per Z plane.
//   angular           prob_angle_z [nZ, ntheta, nphi]: the probability per DVR node (theta_t, phi_k)
//   diffraction       pop_nm_z [nZ, nX, nY]: the probability per parallel-momentum point (FFT order), the initial
//                     parallel momentum shifted out as in the BK analysis
//   vibration         pop_vj_z [nZ, nR, nJ]: the probability per gas-phase state phi_vJ(r) (all nR of each J)
//   rotation_moments  p2_numerator_z, rot_numerator_z [nZ]: <P2(cos theta)> and <J(J+1) / (2 mu r^2)> times the
//                     probability at each Z
//   momentum          prob_kr [nR]: the probability per k_r point (FFT order); prob_kz_spec [nZSp]: per k_Z point
//                     of the specular wave
// Each is summed in a fixed order whatever the number of threads. At each Z the first three sum to the frame's
// probability there; prob_kr sums to the scattering grid's probability, prob_kz_spec to the norm of chi.

#include <complex>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "fft/fft.hpp"
#include "frames/frame_reducers.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "propagation/angular_transform.hpp"
#include "propagation/propagation_tables.hpp"

namespace spodvr {

// Which extras a run computes (FrameExtras).
struct FrameExtraSwitches {
  bool angular = false;
  bool diffraction = false;
  bool vibration = false;
  bool rotation_moments = false;
  bool momentum = false;

  [[nodiscard]] static FrameExtraSwitches of(const InputParameters& p);
  [[nodiscard]] bool any() const { return angular || diffraction || vibration || rotation_moments || momentum; }
};

// The sizes of the extras' datasets (zero for an extra that is off).
struct FrameExtraSizes {
  std::size_t angle_z = 0;
  std::size_t nm_z = 0;
  std::size_t vj_z = 0;
  std::size_t moments_z = 0;
  std::size_t kr = 0;
  std::size_t kz_spec = 0;
};
[[nodiscard]] FrameExtraSizes frame_extra_sizes(const Config& config, const SpatialGrids& grids,
                                                const AngularGrid& angular, const RotationalBasis& basis);

// <J m| P2(cos theta) |J' m> of normalised spherical harmonics (the program's phase: Condon-Shortley for m > 0):
// non-zero for J' = J and J' = J +- 2 only; zero for J or J' below |m|.
[[nodiscard]] double p2_element(int j, int jp, int m);

// The distinct J values of a basis, ascending (the last axis of pop_vj_z).
[[nodiscard]] std::vector<int> basis_j_values(const RotationalBasis& basis);

// The extras of a run: tables and transforms made once, then fed every Z plane of psi_full of a frame.
class FrameExtras {
 public:
  // The tables must outlive the extras (they keep views of shift_forward, rotational_energies and, through the
  // layout, specular_to_scattering); the other arguments are used here only.
  FrameExtras(const Config& config, const SpatialGrids& grids, const AngularGrid& angular, const RotationalBasis& basis,
              const GasPhase& gas, const PropagationTables& tables, const Fft& fft, const FrameLayout& layout);
  ~FrameExtras();
  FrameExtras(const FrameExtras&) = delete;
  FrameExtras& operator=(const FrameExtras&) = delete;

  // The record's extras of Z index iz < nZ from that plane of psi_full, layout (irot, ixy, ir); prob_kr sums over
  // the planes and starts again at iz = 0.
  void reduce_plane(int iz, std::span<const std::complex<double>> full_plane, FrameRecord& record);
  // The record's prob_kz_spec from the specular wave (momentum).
  void reduce_specular(const SpecularWave& specular, FrameRecord& record);

  [[nodiscard]] const FrameExtraSwitches& switches() const { return on_; }
  // The wall-clock time of each extra so far (angular, diffraction, vibration, rotation_moments, momentum).
  [[nodiscard]] const std::vector<double>& seconds() const { return seconds_; }
  // The doubles the extras hold during a frame: their datasets in the record, the plane of psi_full, their tables
  // and their scratch (per thread where each thread has its own), for the memory estimate.
  [[nodiscard]] static std::size_t memory_values(const FrameLayout& layout, const AngularGrid& angular,
                                                 const FrameExtraSwitches& on, int n_j, int n_threads);

 private:
  struct Tables;
  void reduce_angular(int iz, std::span<const std::complex<double>> plane, FrameRecord& record);
  void reduce_diffraction(int iz, std::span<const std::complex<double>> plane, FrameRecord& record);
  void reduce_vibration(int iz, std::span<const std::complex<double>> plane, FrameRecord& record);
  void reduce_rotation_moments(int iz, std::span<const std::complex<double>> plane, FrameRecord& record);
  void reduce_momentum(std::span<const std::complex<double>> plane, FrameRecord& record);

  FrameLayout layout_;
  FrameExtraSwitches on_;
  std::unique_ptr<Tables> tables_;
  std::vector<double> seconds_;
};

}  // namespace spodvr

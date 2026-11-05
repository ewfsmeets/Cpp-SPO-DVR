#pragma once

#include <complex>
#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "analysis/analysis_basis.hpp"
#include "analysis/balint_kurti.hpp"
#include "analysis/s_matrix.hpp"
#include "core/aligned_vector.hpp"
#include "fft/fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/config.hpp"

namespace spodvr {

// Fine analysis (Fortran nread = 2): the S-matrix on the full spectral grid E_k = k pi / T of
// the BK time series of a finished run, and its linear interpolation to the NrOfEs analysis
// energies. Energies E_k are measured from the initial channel energy enymo(inchna, nvin), so
// E_k is the translational energy in Z of the elastic channel.
// Fortran: AsympInfo.f90 ReadOrWriteAsympInfo (nread = 2), Analyzor.f90 gtaife, getbkc,
// hlftrf, shbkar, probaf, StateToStateProbs and smatrix_calc (nread = 2)

// The spectral grid. Grid point k (1-based, as in Fortran) has E_z = k * spres; the grid keeps
// the points first_index .. last_index whose E_z lies in the 95% confidence interval.
// Fortran: ReadOrWriteAsympInfo (ndatol, ndused, notraf, spres, icomin, icomax, nen95)
struct FineAnalysisGrid {
  int n_records_in_file;     // BK records in the file
  int n_records_available;   // = nint(tEnd / deltat) (at most the file's). Fortran: ndatol
  int n_records_used;        // largest FFT length <= that. Fortran: ndused
  int transform_length;      // 2 * ndused. Fortran: notraf
  double time_step;          // nanspo * dt. Fortran: deltat
  double effective_time;     // ndused * deltat. Fortran: ttime
  double spectral_resolution;  // pi / ttime (eV). Fortran: spres
  double meaningful_min;     // eV. Fortran: emintr
  double meaningful_max;     // eV. Fortran: emaxtr
  int first_index;           // Fortran: icomin
  int last_index;            // Fortran: icomax

  [[nodiscard]] int size() const { return last_index >= first_index ? last_index - first_index + 1 : 0; }  // nen95
  // E_z of grid point i = 0 .. size()-1 (eV). Fortran: fetrz(i+1) = ftrans(i+1)
  [[nodiscard]] double z_energy(int i) const { return (first_index + i) * spectral_resolution; }
  [[nodiscard]] std::vector<double> z_energies() const;
};

// Throws InputError if tEnd needs more records than the file has (Fortran: stop) or if fewer
// than 2 records are available (GetLowestFFTSize would return 1).
[[nodiscard]] FineAnalysisGrid make_fine_analysis_grid(const Config& config, int n_records_in_file);

// The BK coefficients of the first ndused records, frequency shifted by the initial channel
// energy: c(j) = balcf(t_j) exp(i ev0jl t_j), t_j = j deltat, c(0) = 0. One contiguous series of
// ndused + 1 values per state k = is + nstata * v (the Fortran bkar layout).
// Fortran: getbkc
class BkSeries {
 public:
  BkSeries(const FineAnalysisGrid& grid, const AnalysisBasis& basis, double initial_channel_energy);

  // Stores record j = 1, 2, ... (records beyond ndused are ignored). `time` is checked:
  // it must be j * deltat to within rounding.
  void add_record(double time, const BkCoefficients& bk);
  // Throws CheckpointError unless all ndused records were added.
  void check_complete() const;

  [[nodiscard]] int n_series() const { return n_series_; }   // nstata * nviba
  [[nodiscard]] int length() const { return n_used_ + 1; }   // ndused + 1 (j = 0 .. ndused)
  [[nodiscard]] std::span<const std::complex<double>> series(int k) const;

 private:
  int n_used_;
  int n_series_;
  double time_step_;
  double energy_;  // ev0jl (hartree)
  int n_added_ = 0;
  AlignedVector<std::complex<double>> values_;
};

// The time-to-energy transform of every series at the grid points of `grid`:
//   A_k(E_i) = deltat * sum_{j=0}^{ndused} w_j c_k(j) exp(i E_i t_j),  w_0 = w_ndused = 1/2, else 1,
// i.e. the trapezoid rule for the integral over [0, T]. Computed with one zero-padded FFT of
// length 2 ndused (sign +1) per series. Fortran hlftrf does the same with two FFTs of the
// Hermitian and anti-Hermitian extensions of c (equal up to rounding; tested).
// Result layout: k + n_series * i. Fortran: hlftrf, shbkar
[[nodiscard]] std::vector<std::complex<double>> half_fourier_transform(const BkSeries& series,
                                                                       const FineAnalysisGrid& grid, const Fft& fft);

// S-matrix, probabilities, sums and averages at the energies z_energies (E_z in eV; one
// EnergyResult each), from the transformed series A_k(E) at those energies (layout
// k + nstata*nviba * i). Fortran: probaf (prefactors; the sums are those of smatrix_combine)
[[nodiscard]] ScatteringResults fine_spectrum(const Config& config, std::span<const double> z_energies,
                                              const AnalysisBasis& basis, const GasPhase& gas,
                                              std::span<const std::complex<double>> amplitudes);

// S at E_z between two grid points (fraction f from point 1 to point 2): |S|^2 is interpolated
// linearly (as Fortran) and the phase is that of the linearly interpolated S.
// The Fortran interpolation of the phase is wrong (ATAN, Re S in place of Im S, signs taken
// separately); see docs/fortran_fine_analysis_bugs.md.
[[nodiscard]] std::complex<double> interpolate_s_matrix(std::complex<double> s1, std::complex<double> s2, double f);

// The grid interval of E_z for the interpolation: upper neighbour (0-based, in 1 .. size-1) and
// fraction f = (E_z - E(upper-1)) / (E(upper) - E(upper-1)). The upper neighbour is the first
// grid point with E >= E_z, clamped (Fortran: idat in 2 .. nen95, plain >=), so f can be
// slightly larger than 1 just above the last point. Needs size() >= 2.
struct GridInterval {
  int upper;
  double f;
};
[[nodiscard]] GridInterval locate_on_grid(const FineAnalysisGrid& grid, double z_energy);

// An analysis energy for which no interpolated result exists.
struct SkippedEnergy {
  double z_energy;      // eV
  std::string reason;
};

struct FineAnalysisResults {
  FineAnalysisGrid grid;
  ScatteringResults spectrum;      // at the grid points
  ScatteringResults interpolated;  // at the analysis energies inside the grid
  std::vector<SkippedEnergy> skipped;
  // Populated channels whose phase changes by more than pi/2 between neighbouring grid points
  // (the interpolation is then unreliable: propagate longer).
  int large_phase_steps = 0;
};

// The interpolated results at the NrOfEs analysis energies (`energies`). Energies outside the
// confidence interval or the grid are skipped (Fortran prints a message). With fewer than 2
// grid points nothing is interpolated. Fortran: StateToStateProbs and smatrix_calc (nread = 2)
[[nodiscard]] FineAnalysisResults interpolate_to_analysis_energies(const Config& config, const FineAnalysisGrid& grid,
                                                                   const AnalysisBasis& basis, const GasPhase& gas,
                                                                   const AnalysisEnergies& energies,
                                                                   ScatteringResults spectrum);

// Counts the populated channels (|S|^2 > 1e-8 of the largest |S|^2) whose phase changes by more
// than pi/2 between neighbouring grid points.
[[nodiscard]] int count_large_phase_steps(const ScatteringResults& spectrum);

}  // namespace spodvr

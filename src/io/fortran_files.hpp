#pragma once

#include <filesystem>
#include <functional>

#include "analysis/analysis_basis.hpp"
#include "analysis/balint_kurti.hpp"
#include "gasphase/gas_phase.hpp"
#include "input/config.hpp"
#include "io/fortran_unformatted.hpp"
#include "propagation/analysis_surface.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr {

// Restart files in the layout of the Fortran code (gfortran unformatted), so that runs can be
// continued across the two programs.

// Wave function file (Fortran writes fort.9, reads fort.8): nproj; SpecWave;
// (FixedZAnalysisSpec, izanlySpecChan); then ScatWave(:, :, izr) for izr = 1 .. nz*nr.
// Fortran: WaveFunction.f90 SaveWaveFunc / ReadWaveFunc
void write_fortran_wave_function(const std::filesystem::path& path, const WavePacket& psi,
                                 const AnalysisSurfaceState& surface);

// The Fortran file holds only the surface position and whether it is fixed. The other state
// (anlyTol, oldValZanly) is not saved by the Fortran code (it is left uninitialised on a
// restart); it is taken from `fresh`, the state of a new analysis surface for the same input,
// with oldValZanly = izanlySpecChan.
[[nodiscard]] AnalysisSurfaceState read_fortran_wave_function(const std::filesystem::path& path, WavePacket& psi,
                                                              const AnalysisSurfaceState& fresh);

// BK coefficient file (Fortran writes fort.21, reads fort.22): ten header records
// (sizes; deltat, rkz, z0, alfwv, ZAnalysis, eviroi; enrgy, thtai, phii, jin, mjin, ichn;
// the six label arrays of the rotation-diffraction basis; enymo), then one record balcf per
// analysis step. Fortran: AsympInfo.f90 ReadOrWriteAsympInfo, Analyzor.f90 BalintKurtiCoeffs
class FortranBkFileWriter {
 public:
  // n_records: the number of BK records the file will hold (Fortran: ndatnw), written first.
  FortranBkFileWriter(const std::filesystem::path& path, const Config& config, const RotationalBasis& rotations,
                      const AnalysisBasis& basis, const GasPhase& gas, int n_records);
  void append(const BkCoefficients& bk);

 private:
  FortranRecordWriter file_;
};

// Calls on_record(time, coefficients) for every record, with time = istep * deltat as in
// Fortran GetSMatrixElemsOfOldRun. The sizes, deltat and basis labels must match the input.
// Returns the number of records.
int read_fortran_bk_file(const std::filesystem::path& path, const Config& config, const RotationalBasis& rotations,
                         const AnalysisBasis& basis,
                         const std::function<void(double, const BkCoefficients&)>& on_record);

}  // namespace spodvr

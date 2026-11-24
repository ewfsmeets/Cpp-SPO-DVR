// Restart files in the Fortran layout. Fortran: WaveFunction.f90 SaveWaveFunc / ReadWaveFunc,
// AsympInfo.f90 ReadOrWriteAsympInfo, RotDifBasis.f90 WriteToFile, Analyzor.f90 BalintKurtiCoeffs
#include "io/fortran_files.hpp"

#include <array>
#include <complex>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace spodvr {

namespace {

using cplx = std::complex<double>;
using int32 = std::int32_t;

void check_value(int32 stored, int expected, const std::string& what, const std::filesystem::path& path) {
  if (stored != expected)
    throw CheckpointError(what + " is " + std::to_string(expected) + " but " + std::to_string(stored) + " in '" +
                          path.string() + "'");
}

// The six label arrays of the rotation-diffraction basis (1-based pointers), in file order.
// Fortran: RotDifBasis.f90 WriteToFile (PntToRot, PntToDif, NQuantum, MQuantum, JQuantum, MjQuantum)
std::array<std::vector<int32>, 6> state_labels(const AnalysisBasis& basis) {
  std::array<std::vector<int32>, 6> labels;
  for (int is = 0; is < basis.n_states(); ++is) {
    const DiffractionChannel& c = basis.channel(is);
    labels[0].push_back(basis.state_rotation(is) + 1);
    labels[1].push_back(basis.state_diffraction(is) + 1);
    labels[2].push_back(c.n);
    labels[3].push_back(c.m);
    labels[4].push_back(basis.j(is));
    labels[5].push_back(basis.mj(is));
  }
  return labels;
}

// Record 3 of a wave-function file, (FixedZAnalysisSpec, izanlySpecChan), on top of `state`
// (oldValZanly = izanlySpecChan). izanlySpecChan must lie on the specular grid (1 .. n_z_specular):
// a corrupt file is refused here, before the scattering wave is read.
AnalysisSurfaceState read_analysis_record(FortranRecordReader& file, const std::filesystem::path& path,
                                          int n_z_specular, AnalysisSurfaceState state) {
  FortranRecordParser analysis = file.read("record 3 (FixedZAnalysisSpec, izanlySpecChan)");
  state.fixed = analysis.read<int32>() != 0;
  const int32 iz_analysis = analysis.read<int32>();  // 1-based; checked before any arithmetic (a corrupt file)
  analysis.expect_end();
  if (iz_analysis < 1 || iz_analysis > n_z_specular)
    throw CheckpointError("'" + path.string() + "': izanlySpecChan " + std::to_string(iz_analysis) +
                          " is outside the specular grid, 1 .. " + std::to_string(n_z_specular));
  state.index = iz_analysis - 1;
  state.previous_index = state.index;
  return state;
}

double analysis_time_step(const Config& config) {
  return config.input.time_step * config.derived.steps_per_analysis;  // Fortran: deltat = dtlan * nanspo
}

}  // namespace

void write_fortran_wave_function(const std::filesystem::path& path, const WavePacket& psi,
                                 const AnalysisSurfaceState& surface) {
  FortranRecordWriter file(path);
  file.write(FortranRecord().add<int32>(2));  // nproj
  file.write_array(std::span<const cplx>(psi.specular.values));
  file.write(FortranRecord().add<int32>(surface.fixed ? 1 : 0).add<int32>(surface.index + 1));
  for (int izr = 0; izr < psi.scattering.n_blocks(); ++izr) file.write_array(psi.scattering.block(izr));
  file.close();
}

AnalysisSurfaceState read_fortran_wave_function(const std::filesystem::path& path, WavePacket& psi,
                                                const AnalysisSurfaceState& fresh) {
  FortranRecordReader file(path);
  FortranRecordParser nproj = file.read("record 1 (nproj)");
  check_value(nproj.read<int32>(), 2, "ProjOp", path);
  nproj.expect_end();

  FortranRecordParser specular = file.read("record 2 (SpecWave)");
  specular.read_array(std::span<cplx>(psi.specular.values));
  specular.expect_end();

  const AnalysisSurfaceState state = read_analysis_record(file, path, psi.specular.n_z(), fresh);

  for (int izr = 0; izr < psi.scattering.n_blocks(); ++izr) {
    FortranRecordParser block = file.read("ScatWave record " + std::to_string(izr + 1));
    block.read_array(psi.scattering.block(izr));
    block.expect_end();
  }
  return state;
}

AnalysisSurfaceState read_fortran_analysis_surface(const std::filesystem::path& path, int n_r, int n_z_specular) {
  FortranRecordReader file(path);
  FortranRecordParser nproj = file.read("record 1 (nproj)");
  check_value(nproj.read<int32>(), 2, "ProjOp", path);
  nproj.expect_end();

  std::vector<cplx> specular_wave(std::size_t(n_r) * n_z_specular);
  FortranRecordParser specular = file.read("record 2 (SpecWave)");
  specular.read_array(std::span<cplx>(specular_wave));
  specular.expect_end();

  return read_analysis_record(file, path, n_z_specular, AnalysisSurfaceState{});
}

FortranBkFileWriter::FortranBkFileWriter(const std::filesystem::path& path, const Config& config,
                                         const RotationalBasis& rotations, const AnalysisBasis& basis,
                                         const GasPhase& gas, int n_records)
    : file_(path) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  const int n_dif_total = static_cast<int>(basis.diffraction_channels().size());
  file_.write(FortranRecord()
                  .add<int32>(n_records)
                  .add<int32>(n_dif_total * rotations.size())  // nstatw
                  .add<int32>(n_dif_total)                     // ndftow
                  .add<int32>(rotations.size()));              // nrotw
  file_.write(FortranRecord()
                  .add(analysis_time_step(config))
                  .add(d.kz)
                  .add(p.z0)
                  .add(d.wave_packet_alpha)
                  .add(p.z_analysis)
                  .add(gas.initial.energy));
  file_.write(FortranRecord()
                  .add(d.translational_energy)
                  .add(d.theta_incidence)
                  .add(d.phi_incidence)
                  .add<int32>(p.j_in)
                  .add<int32>(p.mj_in)
                  .add<int32>(rotations.initial_index() + 1));
  for (const std::vector<int32>& labels : state_labels(basis)) file_.write_array(std::span<const int32>(labels));
  std::vector<double> channel_energies;  // Fortran: enymo(is, v)
  for (int v = 0; v < basis.n_v(); ++v)
    for (int is = 0; is < basis.n_states(); ++is) channel_energies.push_back(basis.channel_energy(is, v));
  file_.write_array(std::span<const double>(channel_energies));
  file_.flush();
}

void FortranBkFileWriter::append(const BkCoefficients& bk) {
  file_.write_array(std::span<const cplx>(bk.values));
  file_.flush();
}

int read_fortran_bk_file(const std::filesystem::path& path, const Config& config, const RotationalBasis& rotations,
                         const AnalysisBasis& basis,
                         const std::function<void(double, const BkCoefficients&)>& on_record) {
  FortranRecordReader file(path);
  const int n_dif_total = static_cast<int>(basis.diffraction_channels().size());

  FortranRecordParser sizes = file.read("record 1 (sizes)");
  const int n_records = sizes.read<int32>();
  check_value(sizes.read<int32>(), n_dif_total * rotations.size(), "nstatw", path);
  check_value(sizes.read<int32>(), n_dif_total, "ndftow", path);
  check_value(sizes.read<int32>(), rotations.size(), "nrotw", path);
  sizes.expect_end();

  FortranRecordParser times = file.read("record 2 (deltat, ...)");
  const double deltat = times.read<double>();
  if (deltat != analysis_time_step(config))
    throw CheckpointError("the analysis time step (dtAna) differs from the one in '" + path.string() + "'");
  (void)file.read("record 3");

  const auto expected = state_labels(basis);
  const std::array<std::string_view, 6> names = {"PntToRot", "PntToDif", "n", "m", "J", "mJ"};
  for (int k = 0; k < 6; ++k) {
    FortranRecordParser labels = file.read(std::string("record ") + std::to_string(4 + k));
    if (labels.size() != expected[k].size() * sizeof(int32))
      throw CheckpointError("the analysis basis (nstata) differs from the one in '" + path.string() + "'");
    for (const int32 label : expected[k])
      check_value(labels.read<int32>(), label, "a state label (" + std::string(names[k]) + ")", path);
  }
  (void)file.read("record 10 (enymo)");

  BkCoefficients bk(basis.n_states(), basis.n_v());
  for (int istep = 1; istep <= n_records; ++istep) {
    FortranRecordParser record = file.read("BK record " + std::to_string(istep));
    record.read_array(std::span<cplx>(bk.values));
    record.expect_end();
    on_record(istep * deltat, bk);  // Fortran: GetSMatrixElemsOfOldRun
  }
  return n_records;
}

}  // namespace spodvr

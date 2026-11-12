// Restart files in the Fortran layout. Fortran: WaveFunction.f90 SaveWaveFunc / ReadWaveFunc,
// AsympInfo.f90 ReadOrWriteAsympInfo, RotDifBasis.f90 WriteToFile, Analyzor.f90 BalintKurtiCoeffs,
// FluxAnalyzor.f90 writeFluxCoeffs / readFluxCoeffs
#include "io/fortran_files.hpp"

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <format>
#include <span>
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

  FortranRecordParser analysis = file.read("record 3 (FixedZAnalysisSpec, izanlySpecChan)");
  AnalysisSurfaceState state = fresh;
  state.fixed = analysis.read<int32>() != 0;
  state.index = analysis.read<int32>() - 1;
  state.previous_index = state.index;
  analysis.expect_end();

  for (int izr = 0; izr < psi.scattering.n_blocks(); ++izr) {
    FortranRecordParser block = file.read("ScatWave record " + std::to_string(izr + 1));
    block.read_array(psi.scattering.block(izr));
    block.expect_end();
  }
  return state;
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

namespace {

// The record order of the flux files, (iz, n) with n outer, is the memory order of (block, iz, n):
// record k holds the n_blocks values from k n_blocks on.
std::size_t flux_record_count(const FluxCoefficients& c) { return c.index.n1 * c.index.n2; }

// Fortran: FluxAnalyzor.f90 writeFluxCoeffs, write(funit) phi(:, iz, ien) for ien, then iz
void write_flux_array(const std::filesystem::path& path, const FluxCoefficients& c, std::span<const cplx> values) {
  FortranRecordWriter file(path);
  for (std::size_t k = 0; k < flux_record_count(c); ++k) file.write_array(values.subspan(k * c.index.n0, c.index.n0));
  file.close();
}

std::string flux_count_error(const std::filesystem::path& path, const std::string& found, const FluxCoefficients& c) {
  return std::format(
      "'{}' holds {} flux coefficient records, expected {} ({} flux energies x {} planes Z < zasymp): "
      "EMinFluxAna, EMaxFluxAna, DeltaEFluxAna or the Z grid differ from the run that wrote it",
      path.string(), found, flux_record_count(c), c.n_energies(), c.n_planes());
}

// Fortran: FluxAnalyzor.f90 readFluxCoeffs, read(funit) phi(:, iz, ien) for ien, then iz. The
// Fortran accepts longer records and extra records; here the file must fit the input exactly (flux Q8).
void read_flux_array(const std::filesystem::path& path, const FluxCoefficients& c, std::span<cplx> values) {
  FortranRecordReader file(path);
  const std::size_t record_bytes = c.index.n0 * sizeof(cplx);
  for (std::size_t k = 0; k < flux_record_count(c); ++k) {
    if (file.at_end()) throw CheckpointError(flux_count_error(path, std::to_string(k), c));
    FortranRecordParser record = file.read(std::format("flux coefficient record {}", k + 1));
    if (record.size() != record_bytes)
      throw CheckpointError(std::format(
          "record {} of '{}' has {} bytes, expected {} ({} complex values: the rotational functions x nX x nY "
          "of the input)",
          k + 1, path.string(), record.size(), record_bytes, c.n_blocks()));
    record.read_array(values.subspan(k * c.index.n0, c.index.n0));
  }
  if (!file.at_end())
    throw CheckpointError(flux_count_error(path, "more than " + std::to_string(flux_record_count(c)), c));
}

}  // namespace

void write_fortran_flux_coefficients(const std::filesystem::path& phi_path,
                                     const std::filesystem::path& derivative_path,
                                     const FluxCoefficients& coefficients) {
  write_flux_array(phi_path, coefficients, coefficients.phi);
  write_flux_array(derivative_path, coefficients, coefficients.phi_derivative);
}

void read_fortran_flux_coefficients(const std::filesystem::path& phi_path, const std::filesystem::path& derivative_path,
                                    FluxCoefficients& coefficients) {
  read_flux_array(phi_path, coefficients, coefficients.phi);
  read_flux_array(derivative_path, coefficients, coefficients.phi_derivative);
}

}  // namespace spodvr

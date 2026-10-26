// Native checkpoint (wave function) and BK coefficient files.
// Fortran: WaveFunction.f90 SaveWaveFunc / ReadWaveFunc, AsympInfo.f90 ReadOrWriteAsympInfo
#include "io/checkpoint.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "report/format.hpp"

namespace spodvr {

namespace {

using cplx = std::complex<double>;
using Magic = std::array<char, 16>;

constexpr Magic checkpoint_magic = {'S', 'P', 'O', 'D', 'V', 'R', '-', 'C', 'H', 'E', 'C', 'K', 'P', 'N', 'T', '\0'};
constexpr Magic bk_magic = {'S', 'P', 'O', 'D', 'V', 'R', '-', 'B', 'K', 'C', 'O', 'E', 'F', 'F', 'S', '\0'};
constexpr std::uint32_t format_version = 1;
constexpr std::uint64_t end_marker = 0x454e44'4f46'4350ULL;  // marks a completely written checkpoint

// Input parameters that a restart must not change, with their JSON key names.
std::vector<std::pair<std::string, double>> restart_parameters(const Config& config) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  return {{"ZStart", p.z_start},
          {"dZ", p.dz},
          {"nZ", p.nz},
          {"nZSpecular", p.nz_specular},
          {"RStart", p.r_start},
          {"dR", p.dr},
          {"nR", p.nr},
          {"nX", p.nx},
          {"nY", p.ny},
          {"JMax", p.j_max},
          {"MJMax", d.mj_max_basis},
          {"Mass1", p.mass1},
          {"Mass2", p.mass2},
          {"HomoNuc", p.homonuclear ? 1.0 : 0.0},
          {"dtStep", p.time_step},
          {"dtAna", p.dt_analysis},
          {"Zana", p.z_analysis},
          {"Z0", p.z0},
          {"EMinWP", p.e_min_wave_packet},
          {"EMaxWP", p.e_max_wave_packet},
          {"EParallel", p.e_parallel},
          {"PhiParallel", p.phi_parallel},
          {"VIn", p.v_in},
          {"JIn", p.j_in},
          {"MJIn", p.mj_in},
          {"ZOptStart", p.z_opt_start},
          {"ZOptEnd", p.z_opt_end},
          {"ROptStart", p.r_opt_start},
          {"ROptEnd", p.r_opt_end},
          {"ZSpOptStart", p.zsp_opt_start},
          {"ZSpOptEnd", p.zsp_opt_end},
          {"ZOptFactor (effective)", d.z_opt_factor},
          {"ROptFactor (effective)", d.r_opt_factor},
          {"ZSpOptFactor (effective)", d.zsp_opt_factor}};
}

void write_header(BinaryWriter& file, const Magic& magic, const Config& config) {
  file.write(magic);
  file.write(format_version);
  file.write_string(config.input.system_id);
  const auto parameters = restart_parameters(config);
  file.write(static_cast<std::uint64_t>(parameters.size()));
  for (const auto& [name, value] : parameters) {
    file.write_string(name);
    file.write(value);
  }
}

// Reads and checks the header: file type, version and the restart parameters.
void read_header(BinaryReader& file, const Magic& magic, const std::string& kind, const Config& config) {
  const std::string where = " in " + kind + " '" + file.path().string() + "'";
  if (file.read<Magic>() != magic) throw CheckpointError("not a " + kind + ": '" + file.path().string() + "'");
  if (const auto version = file.read<std::uint32_t>(); version != format_version)
    throw CheckpointError("unsupported format version " + std::to_string(version) + where);
  if (const std::string id = file.read_string(); id != config.input.system_id)
    throw CheckpointError("SystemID is '" + config.input.system_id + "' in the input but '" + id + "'" + where);

  const auto expected = restart_parameters(config);
  if (file.read<std::uint64_t>() != expected.size()) throw CheckpointError("unexpected parameter list" + where);
  for (const auto& [name, value] : expected) {
    const std::string stored_name = file.read_string();
    const double stored = file.read<double>();
    if (stored_name != name) throw CheckpointError("unexpected parameter list" + where);
    if (stored != value)
      throw CheckpointError(name + " is " + full_precision(value) + " in the input but " + full_precision(stored) +
                            where + "; a restart must use the same parameters");
  }
}

void check_size(BinaryReader& file, const std::string& what, std::int32_t expected) {
  if (const auto stored = file.read<std::int32_t>(); stored != expected)
    throw CheckpointError(what + " is " + std::to_string(expected) + " but " + std::to_string(stored) + " in '" +
                          file.path().string() + "'");
}

}  // namespace

void write_checkpoint(const std::filesystem::path& path, const Config& config, double time, const WavePacket& psi,
                      const AnalysisSurfaceState& surface) {
  std::filesystem::path temporary = path;
  temporary += ".tmp";
  {
    BinaryWriter file(temporary);
    write_header(file, checkpoint_magic, config);
    file.write(time);
    for (const std::size_t n : {psi.scattering.index.n0, psi.scattering.index.n1, psi.scattering.index.n2,
                                psi.specular.index.n0, psi.specular.index.n1})
      file.write(static_cast<std::int32_t>(n));
    file.write(static_cast<std::int32_t>(surface.index));
    file.write(static_cast<std::int32_t>(surface.previous_index));
    file.write(static_cast<std::uint8_t>(surface.fixed ? 1 : 0));
    file.write(surface.norm_tolerance);
    file.write_array(std::span<const cplx>(psi.specular.values));
    file.write_array(std::span<const cplx>(psi.scattering.values));
    file.write(end_marker);
    file.close();
  }
  std::filesystem::rename(temporary, path);
}

CheckpointInfo read_checkpoint(const std::filesystem::path& path, const Config& config, WavePacket& psi) {
  BinaryReader file(path);
  read_header(file, checkpoint_magic, "checkpoint", config);
  CheckpointInfo info{};
  info.time = file.read<double>();
  check_size(file, "the number of rotational states", static_cast<std::int32_t>(psi.scattering.index.n0));
  check_size(file, "the number of (x, y) points", static_cast<std::int32_t>(psi.scattering.index.n1));
  check_size(file, "the number of (Z, r) points", static_cast<std::int32_t>(psi.scattering.index.n2));
  check_size(file, "nR", static_cast<std::int32_t>(psi.specular.index.n0));
  check_size(file, "nZSpecular", static_cast<std::int32_t>(psi.specular.index.n1));
  info.surface.index = file.read<std::int32_t>();
  info.surface.previous_index = file.read<std::int32_t>();
  info.surface.fixed = file.read<std::uint8_t>() != 0;
  info.surface.norm_tolerance = file.read<double>();
  file.read_array(std::span<cplx>(psi.specular.values));
  file.read_array(std::span<cplx>(psi.scattering.values));
  if (file.read<std::uint64_t>() != end_marker || !file.at_end())
    throw CheckpointError("checkpoint '" + path.string() + "' is corrupt (bad end marker)");
  return info;
}

BkFileWriter::BkFileWriter(const std::filesystem::path& path, const Config& config, const AnalysisBasis& basis)
    : file_(path) {
  write_header(file_, bk_magic, config);
  file_.write(static_cast<std::int32_t>(basis.n_states()));
  file_.write(static_cast<std::int32_t>(basis.n_v()));
  for (int is = 0; is < basis.n_states(); ++is) {
    const DiffractionChannel& c = basis.channel(is);
    for (const int label : {basis.j(is), basis.mj(is), c.n, c.m, c.order})
      file_.write(static_cast<std::int32_t>(label));
  }
  for (int v = 0; v < basis.n_v(); ++v)
    for (int is = 0; is < basis.n_states(); ++is) file_.write(basis.channel_energy(is, v));
  file_.flush();
}

void BkFileWriter::append(double time, const BkCoefficients& bk) {
  file_.write(time);
  file_.write_array(std::span<const cplx>(bk.values));
  file_.flush();
}

std::optional<double> read_bk_file(const std::filesystem::path& path, const Config& config, const AnalysisBasis& basis,
                                   const std::function<void(double, const BkCoefficients&)>& on_record) {
  BinaryReader file(path);
  read_header(file, bk_magic, "BK coefficient file", config);
  check_size(file, "the number of analysed states (nstata)", basis.n_states());
  check_size(file, "the number of analysed vibrational levels (nviba)", basis.n_v());
  for (int is = 0; is < basis.n_states(); ++is) {
    const DiffractionChannel& c = basis.channel(is);
    for (const int label : {basis.j(is), basis.mj(is), c.n, c.m, c.order})
      check_size(file, "a state label of the analysis basis", label);
  }
  for (int k = 0; k < basis.n_states() * basis.n_v(); ++k) (void)file.read<double>();  // channel energies

  std::optional<double> last_time;
  BkCoefficients bk(basis.n_states(), basis.n_v());
  while (!file.at_end()) {
    const double time = file.read<double>();
    file.read_array(std::span<cplx>(bk.values));
    on_record(time, bk);
    last_time = time;
  }
  return last_time;
}

}  // namespace spodvr

// Printed sections and full-precision data files of the setup results.
#include "report/sections.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

#include "core/unit_conversion.hpp"

namespace spodvr {

namespace {

std::string yes_no(bool b) { return b ? "yes" : "no"; }

// Fortran nread: 0 new, 1 restart, 2 analysis only, 6 potential only.
std::string run_mode_name(int run_mode) {
  switch (run_mode) {
    case 0: return "new calculation";
    case 1: return "restart";
    case 2: return "analysis only";
    case 6: return "potential only";
    default: return "unknown (" + std::to_string(run_mode) + ")";
  }
}

std::ofstream open_output(const std::filesystem::path& path) {
  std::ofstream file(path);
  if (!file) throw std::runtime_error("cannot open output file '" + path.string() + "'");
  return file;
}

}  // namespace

void print_input_summary(std::ostream& out, const Style& style, const Config& config) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  print_section(out, style, "INPUT");
  KeyValueList list;
  list.add("Title", p.title);
  list.add("System", p.system_id);
  list.add("Run mode", run_mode_name(d.run_mode));
  list.add("Masses", fixed(p.mass1, 6) + ", " + fixed(p.mass2, 6), "amu");
  list.add("Homonuclear", yes_no(p.homonuclear));
  list.add("Total mass", full_precision(d.total_mass), "m_e");
  list.add("Reduced mass", full_precision(d.reduced_mass), "m_e");
  list.add("Lattice constants (x, y)",
           fixed(d.system.lattice_constant_x, 6) + ", " + fixed(d.system.lattice_constant_y, 6), "bohr");
  list.add("Skewing angle", fixed(rad_to_deg(d.system.skewing_angle), 6), "deg");
  list.add("Point group", d.system.point_group);
  list.add("Z asymptotic", fixed(d.z_asymptotic, 6), "bohr");
  list.add("Switching to gas phase",
           p.use_switching ? fixed(p.z_min_switch, 4) + " .. " + fixed(p.z_max_switch, 4) + " bohr" : "no");
  list.add("Time step", fixed(p.time_step, 6), "a.u.");
  list.add("Steps (total / per analysis)", std::to_string(p.n_steps) + " / " + std::to_string(d.steps_per_analysis));
  list.add("Translational energy", full_precision(d.translational_energy), "eV");
  list.add("Wave packet alpha", full_precision(d.wave_packet_alpha), "a.u.");
  list.add("Initial kz", full_precision(d.kz), "a.u.");
  list.add("95% confidence interval", fixed(d.confidence_e_min, 6) + " .. " + fixed(d.confidence_e_max, 6), "eV");
  list.add("Max. analysis energy", fixed(d.max_analysis_energy, 6), "eV");
  list.add(
      "CAP strengths (Z, r, Zsp)",
      scientific(d.z_opt_factor, 8) + ", " + scientific(d.r_opt_factor, 8) + ", " + scientific(d.zsp_opt_factor, 8),
      "hartree");
  list.print(out, style);
}

void print_grids(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& g,
                 const AngularGrid& angular) {
  const auto& units = units::legacy;
  print_section(out, style, "GRIDS");
  KeyValueList list;
  auto grid_text = [](const UniformGrid& grid) {
    return std::to_string(grid.size()) + " points, " + fixed(grid[0], 6) + " .. " + fixed(grid.back(), 6) + ", step " +
           fixed(grid.spacing(), 6);
  };
  list.add("Z (scattering)", grid_text(g.z_scattering), "bohr");
  list.add("Z (specular)", grid_text(g.z_specular), "bohr");
  list.add("r", grid_text(g.r), "bohr");
  list.add("x", grid_text(g.x), "bohr");
  list.add("y", grid_text(g.y), "bohr");
  list.add("dkz (scattering / specular)",
           full_precision(g.kz_scattering.spacing()) + " / " + full_precision(g.kz_specular.spacing()), "1/bohr");
  list.add("dkr", full_precision(g.kr.spacing()), "1/bohr");
  list.add("dkx, dky", full_precision(config.derived.delta_kx) + ", " + full_precision(config.derived.delta_ky),
           "1/bohr");
  list.add("Max. kinetic energy in r", fixed(units.hartree_to_ev(g.kr.max_kinetic_energy()), 6), "eV");
  list.add("Analysis at Z (index)", fixed(g.z_analysis, 6) + " (" + std::to_string(g.iz_analysis) + ")", "bohr");
  list.add("Z points in interaction region", std::to_string(g.nz_potential));
  list.add("Optical potential Z (indices)",
           std::to_string(g.z_optical.first) + " .. " + std::to_string(g.z_optical.last));
  if (g.z_specular_optical)
    list.add("Optical potential Z specular (indices)",
             std::to_string(g.z_specular_optical->first) + " .. " + std::to_string(g.z_specular_optical->last));
  list.add("Optical potential r (indices)",
           std::to_string(g.r_optical.first) + " .. " + std::to_string(g.r_optical.last));
  list.add("Theta points", std::to_string(angular.n_theta()));
  list.add("Phi points", std::to_string(angular.n_phi()));
  list.add("J range (step)", std::to_string(angular.j_min()) + " .. " + std::to_string(angular.j_max()) + " (" +
                                 std::to_string(angular.delta_j()) + ")");
  list.add("MJMax (basis)", std::to_string(angular.mj_max()));
  list.print(out, style);
  out << "\n  Indices are 0-based.\n\n";

  Table table({"itheta", "theta (rad)", "weight"});
  for (int t = 0; t < angular.n_theta(); ++t)
    table.add_row({std::to_string(t), full_precision(angular.thetas()[t]), full_precision(angular.theta_weights()[t])});
  table.print(out, style);
}

void print_memory_estimate(std::ostream& out, const Style& style, const MemoryEstimate& estimate,
                           std::optional<std::size_t> physical_bytes) {
  print_section(out, style, "MEMORY ESTIMATE");
  Table table({"array", "values", "size", "note"});
  for (const MemoryItem& item : estimate.items)
    table.add_row({item.name, std::to_string(item.values), format_bytes(item.bytes), item.note});
  table.print(out, style);
  out << "\n";
  KeyValueList list;
  if (estimate.propagation_bytes == 0) {  // potential-only run
    list.add("Total", format_bytes(estimate.setup_bytes));
  } else {
    list.add("Setting up the propagation", format_bytes(estimate.setup_bytes));
    list.add("During the propagation", format_bytes(estimate.propagation_bytes));
  }
  list.add("Physical memory", physical_bytes ? format_bytes(*physical_bytes) : "unknown");
  list.print(out, style);
  out << "\n  Small tables and FFT plans are not included.\n";
  if (physical_bytes && estimate.peak_bytes > *physical_bytes)
    out << "  " << style.warning("Warning") << ": the run needs " << format_bytes(estimate.peak_bytes)
        << ", more than the physical memory.\n";
}

void print_asymptotic_potential(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                                const PotentialGrid& potential) {
  const auto& units = units::legacy;
  print_section(out, style, "GAS-PHASE POTENTIAL");
  out << "  V(r) at Z = " << fixed(config.derived.z_asymptotic, 6) << " bohr, theta = pi/2\n\n";
  Table table({"ir", "r (bohr)", "V (eV)"});
  for (int ir = 0; ir < grids.r.size(); ++ir)
    table.add_row({std::to_string(ir), full_precision(grids.r[ir]),
                   full_precision(units.hartree_to_ev(potential.asymptotic[ir]))});
  table.print(out, style);
}

void print_gas_phase(std::ostream& out, const Style& style, const Config& config, const GasPhase& gas) {
  const auto& units = units::legacy;
  print_section(out, style, "GAS-PHASE ROVIBRATIONAL STATES");
  Table table({"v", "j", "E (eV)", "E - E(vin,jin) (eV)", "psi(r_first)", "psi(r_last)", "norm", "<r> (bohr)"});
  for (const auto& s : gas.summary)
    table.add_row({std::to_string(s.v), std::to_string(s.j), full_precision(units.hartree_to_ev(s.energy)),
                   full_precision(units.hartree_to_ev(s.relative_energy)), scientific(s.psi_first, 8),
                   scientific(s.psi_last, 8), fixed(s.norm, 12), fixed(s.r_mean, 12)});
  table.print(out, style);
  out << '\n';
  KeyValueList list;
  list.add("Open vibrational levels", std::to_string(gas.v_max_open));
  list.add("Highest open j", std::to_string(gas.j_max_open));
  list.add("Vibrational levels analysed", std::to_string(gas.n_v_analysis));
  list.add("Highest j analysed", std::to_string(gas.j_max_analysis));
  list.add("Initial state (v, j)", std::to_string(config.input.v_in) + ", " + std::to_string(config.input.j_in));
  list.add("  energy", full_precision(units.hartree_to_ev(gas.initial.energy)), "eV");
  list.add("  potential energy", full_precision(units.hartree_to_ev(gas.initial.potential_energy)), "eV");
  list.add("  rotational energy", full_precision(units.hartree_to_ev(gas.initial.rotational_energy)), "eV");
  list.add("  vibrational energy", full_precision(units.hartree_to_ev(gas.initial.vibrational_energy)), "eV");
  list.add("  rovibrational kinetic energy",
           full_precision(units.hartree_to_ev(gas.initial.rovibrational_kinetic_energy)), "eV");
  list.add("  <r>", full_precision(gas.initial.r_mean), "bohr");
  list.print(out, style);
  for (const auto& w : gas.warnings) out << "  " << style.warning("WARNING: " + w) << '\n';
}

void write_asymptotic_potential_file(const std::filesystem::path& path, const SpatialGrids& grids,
                                     const PotentialGrid& potential) {
  const auto& units = units::legacy;
  auto file = open_output(path);
  file << "# Gas-phase potential V(r) (Z = zasymp, theta = pi/2)\n# ir  r(bohr)  V(hartree)  V(eV) = "
          "V(hartree)*27.2114\n";
  for (int ir = 0; ir < grids.r.size(); ++ir) {
    const double v = potential.asymptotic[ir];
    file << ir << ' ' << full_precision(grids.r[ir]) << ' ' << full_precision(v) << ' '
         << full_precision(units.hartree_to_ev(v)) << '\n';
  }
}

void write_gas_phase_file(const std::filesystem::path& path, const SpatialGrids& grids, const GasPhase& gas) {
  auto file = open_output(path);
  const auto& states = gas.states;
  file << "# Gas-phase rovibrational states (atomic units)\n# r grid:";
  for (int ir = 0; ir < grids.r.size(); ++ir) file << ' ' << full_precision(grids.r[ir]);
  file << "\n# v j E(hartree) psi(r_0) ... psi(r_{nr-1})\n";
  for (int v = 0; v < states.n_v(); ++v) {
    for (int j = 0; j <= states.j_max(); ++j) {
      file << v << ' ' << j << ' ' << full_precision(states.energy(v, j));
      for (const double y : states.wave_function(v, j)) file << ' ' << full_precision(y);
      file << '\n';
    }
  }
}

}  // namespace spodvr

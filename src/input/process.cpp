#include "input/process.hpp"

#include <cmath>

#include "core/constants.hpp"
#include "core/unit_conversion.hpp"
#include "input/validate.hpp"
#include "input/wave_packet_params.hpp"
#include "potential/potential_factory.hpp"

namespace spodvr {

namespace {

// Fortran: InputNew.f90 ProcessInput
void process_input(const InputParameters& p, DerivedParameters& d) {
  const auto& units = units::legacy;

  // Run mode. Fortran: nread
  d.run_mode = p.new_calculation ? 0 : 1;
  if (p.only_analysis) d.run_mode = 2;
  if (p.only_potential) d.run_mode = 6;

  d.v_in_index = p.v_in + 1;
  d.mj_max_basis = p.mj_max < 1 ? p.j_max : p.mj_max;

  const double total_mass = units.amu_to_me(p.mass1 + p.mass2);
  const double reduced_mass = units.amu_to_me((p.mass1 * p.mass2) / (p.mass1 + p.mass2));

  d.optical_potential_power = 2;
  d.z_opt_factor = p.z_opt_factor.value_or(0.0);
  d.r_opt_factor = p.r_opt_factor.value_or(0.0);
  d.zsp_opt_factor = p.zsp_opt_factor.value_or(0.0);
  if (p.z_opt_energy && *p.z_opt_energy > 0.0)
    d.z_opt_factor = optical_potential_strength(units.ev_to_hartree(*p.z_opt_energy),
                                                p.z_opt_end - p.z_opt_start, total_mass);
  if (p.r_opt_energy && *p.r_opt_energy > 0.0)
    d.r_opt_factor = optical_potential_strength(units.ev_to_hartree(*p.r_opt_energy),
                                                p.r_opt_end - p.r_opt_start, reduced_mass);
  if (p.zsp_opt_energy && *p.zsp_opt_energy > 0.0)
    d.zsp_opt_factor = optical_potential_strength(units.ev_to_hartree(*p.zsp_opt_energy),
                                                  p.zsp_opt_end - p.zsp_opt_start, total_mass);

  d.rotational_constant = units.wavenumber_to_hartree(59.304);

  d.dx = d.system.lattice_constant_x / static_cast<double>(p.nx);
  d.dy = d.system.lattice_constant_y / static_cast<double>(p.ny);

  d.translational_energy =
      p.e_parallel + units.hartree_to_ev(wave_packet_average_energy(
                         units.ev_to_hartree(p.e_min_wave_packet),
                         units.ev_to_hartree(p.e_max_wave_packet), total_mass));
  d.wave_packet_alpha = wave_packet_alpha(units.ev_to_hartree(p.e_min_wave_packet),
                                          units.ev_to_hartree(p.e_max_wave_packet), total_mass);
  // In degrees here; converted to radians in the Setup part below (as in the Fortran code).
  d.theta_incidence = rad_to_deg(std::asin(std::sqrt(p.e_parallel / d.translational_energy)));
  d.phi_incidence = p.phi_parallel;

  d.analysis_energy_start = p.e_min_wave_packet + p.e_parallel;
  d.analysis_energy_step =
      (p.e_max_wave_packet - p.e_min_wave_packet) / static_cast<double>(p.n_energies - 1);

  d.steps_per_analysis = static_cast<int>(p.dt_analysis / p.time_step);
  d.analyses_per_average = static_cast<int>(p.dt_averages / p.dt_analysis);

  d.fixed_z_analysis_specular = !p.dynamical_zsp_analysis;
}

// Fortran: Setup.f90 MiscellaneousSetup (everything except the grids)
void miscellaneous_setup(const InputParameters& p, DerivedParameters& d) {
  const auto& units = units::legacy;
  constexpr double pi = Constants::pi;
  const double skew = d.system.skewing_angle;
  const double alat = d.system.lattice_constant_x;
  const double blat = d.system.lattice_constant_y;

  d.total_mass = units.amu_to_me(p.mass1 + p.mass2);
  d.reduced_mass = units.amu_to_me(p.mass1 * p.mass2 / (p.mass1 + p.mass2));

  d.delta_kx = 2.0 * pi / (alat * std::sin(skew));
  d.delta_ky = 2.0 * pi / (blat * std::sin(skew));

  d.theta_incidence = pi / 180.0 * d.theta_incidence;
  d.phi_incidence = pi / 180.0 * d.phi_incidence;
  d.k_total = std::sqrt(units.ev_to_hartree(2.0 * d.total_mass * d.translational_energy));

  const double kx_cartesian = d.k_total * std::sin(d.theta_incidence) * std::cos(d.phi_incidence);
  const double ky_cartesian = d.k_total * std::sin(d.theta_incidence) * std::sin(d.phi_incidence);
  d.kx = kx_cartesian / std::sin(skew);
  d.ky = (std::cos(skew) * kx_cartesian + std::sin(skew) * ky_cartesian) / std::sin(skew);
  d.kz = d.k_total * std::cos(d.theta_incidence);

  // 6D calculation (nmoder = 4).
  d.unit_cell_area = alat * blat * std::sin(skew);
  d.unit_cell_norm = std::sqrt(1.0 / d.unit_cell_area);
  d.xy_volume_element = d.dx * d.dy * std::sin(skew);
  d.xyz_volume_element = d.xy_volume_element * p.dz;

  d.v_min = units.ev_to_hartree(p.v_min);
  d.v_max = units.ev_to_hartree(p.v_max);

  d.analysis_time_step = static_cast<double>(d.steps_per_analysis) * p.time_step;
  d.max_analysis_energy = units.hartree_to_ev(pi / d.analysis_time_step);

  const double k_max = d.kz + 2.0 * std::sqrt(d.wave_packet_alpha);
  const double k_min = d.kz - 2.0 * std::sqrt(d.wave_packet_alpha);
  d.confidence_e_max = 0.5 / d.total_mass * (k_max * k_max) * units.hartree_to_ev(1.0);
  d.confidence_e_min = 0.5 / d.total_mass * (k_min * k_min) * units.hartree_to_ev(1.0);
}

}  // namespace

Config build_config(const InputParameters& input) {
  Config config{input, {}};
  DerivedParameters& d = config.derived;

  // Fortran: InputNew.f90 SetupInput
  d.system = make_potential(input.system_id)->system_parameters();
  d.z_asymptotic_unswitched = d.system.z_asymptotic;
  d.z_asymptotic = input.use_switching ? input.z_max_switch : d.system.z_asymptotic;

  validate_input(config.input, d.system, d.z_asymptotic);
  process_input(config.input, d);
  miscellaneous_setup(config.input, d);
  return config;
}

}  // namespace spodvr

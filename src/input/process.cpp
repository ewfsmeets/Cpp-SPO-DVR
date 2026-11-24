// Derived parameters from the validated input. Fortran: InputNew.f90 ProcessInput,
// Setup.f90 MiscellaneousSetup (scalar part)
#include "input/process.hpp"

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "core/constants.hpp"
#include "core/unit_conversion.hpp"
#include "input/validate.hpp"
#include "input/validate_outputs.hpp"
#include "input/wave_packet_params.hpp"
#include "potential/potential_factory.hpp"
#include "potential/scm/scm_crp_potential.hpp"
#include "symmetry/point_group.hpp"

namespace spodvr {

namespace {

// Strength (hartree) of an absorber of the given length tuned for energy_ev (eV). An absorber
// whose L/λ is not inside the table is added to `warnings` (outside it the Fortran warns).
// Fortran: MyFunctions.f90 getACoeffOptPotential, called from InputNew.f90 ProcessInput
double strength_for_energy(const std::string& energy_key, double energy_ev, double length, double mass,
                           std::vector<AbsorberTableWarning>& warnings) {
  const double ekin = units::legacy.ev_to_hartree(energy_ev);
  const AbsorberReducedLength reduced = absorber_reduced_length(ekin, length, mass);
  if (reduced.range != AbsorberTableRange::inside) warnings.push_back({energy_key, reduced.value, reduced.range});
  return optical_potential_strength(ekin, reduced);
}

// Strength of a low-Z absorber (hartree): from the energy if it is > 0, else the factor; 0 when
// the absorber is off or its energy is 0 (as in the Fortran, whose default strength is 0 here).
// Fortran: InputNew.f90 ProcessInput (vabszblow, vabssblow)
double low_z_absorber_strength(const std::string& energy_key, const LowZAbsorberInput& absorber, double total_mass,
                               std::vector<AbsorberTableWarning>& warnings) {
  if (!absorber.enabled()) return 0.0;
  if (absorber.energy && *absorber.energy > 0.0)
    return strength_for_energy(energy_key, *absorber.energy, *absorber.end - *absorber.start, total_mass, warnings);
  return absorber.factor.value_or(0.0);
}

// Fortran: InputNew.f90 ProcessInput
void process_input(const InputParameters& p, DerivedParameters& d) {
  const auto& units = units::legacy;

  // Run mode. Fortran: nread
  d.run_mode = p.new_calculation ? 0 : 1;
  if (p.only_analysis) d.run_mode = 2;
  if (p.only_potential) d.run_mode = 6;

  d.v_in_index = p.v_in + 1;
  d.mj_max_basis = p.mj_max < 1 ? p.j_max : p.mj_max;
  d.phi_points_multiple =
      p.use_symmetry && p.phis_optimal ? phi_points_multiple(PointGroup(d.system.point_group).operations()) : 1;

  const double total_mass = units.amu_to_me(p.mass1 + p.mass2);
  const double reduced_mass = units.amu_to_me((p.mass1 * p.mass2) / (p.mass1 + p.mass2));

  d.optical_potential_power = 2;
  // Without an energy the factor is the strength (0 when absent). The high-end energies are > 0 by the
  // input check (r1718); only ZOptLowEnergy may be 0, which gives strength 0 (absorber D1).
  d.z_opt_factor = p.z_opt_factor.value_or(0.0);
  d.r_opt_factor = p.r_opt_factor.value_or(0.0);
  d.zsp_opt_factor = p.zsp_opt_factor.value_or(0.0);
  // In the order of the Fortran's ProcessInput, which is the order of its warnings.
  std::vector<AbsorberTableWarning>& warnings = d.absorber_table_warnings;
  if (p.z_opt_energy && *p.z_opt_energy > 0.0)
    d.z_opt_factor =
        strength_for_energy("ZOptEnergy", *p.z_opt_energy, p.z_opt_end - p.z_opt_start, total_mass, warnings);
  d.z_opt_low_factor = low_z_absorber_strength("ZOptLowEnergy", p.z_opt_low, total_mass, warnings);
  if (p.r_opt_energy && *p.r_opt_energy > 0.0)
    d.r_opt_factor =
        strength_for_energy("ROptEnergy", *p.r_opt_energy, p.r_opt_end - p.r_opt_start, reduced_mass, warnings);
  if (p.zsp_opt_energy && *p.zsp_opt_energy > 0.0)
    d.zsp_opt_factor =
        strength_for_energy("ZSpOptEnergy", *p.zsp_opt_energy, p.zsp_opt_end - p.zsp_opt_start, total_mass, warnings);
  d.zsp_opt_low_factor = low_z_absorber_strength("ZSpOptLowEnergy", p.zsp_opt_low, total_mass, warnings);

  d.rotational_constant = units.wavenumber_to_hartree(59.304);

  d.dx = d.system.lattice_constant_x / p.nx;
  d.dy = d.system.lattice_constant_y / p.ny;

  d.translational_energy = p.e_parallel + units.hartree_to_ev(wave_packet_average_energy(
                                              units.ev_to_hartree(p.e_min_wave_packet),
                                              units.ev_to_hartree(p.e_max_wave_packet), total_mass));
  d.wave_packet_alpha =
      wave_packet_alpha(units.ev_to_hartree(p.e_min_wave_packet), units.ev_to_hartree(p.e_max_wave_packet), total_mass);
  // In degrees here; converted to radians in the Setup part below (as in the Fortran code).
  d.theta_incidence = rad_to_deg(std::asin(std::sqrt(p.e_parallel / d.translational_energy)));
  d.phi_incidence = p.phi_parallel;

  d.analysis_energy_start = p.e_min_wave_packet + p.e_parallel;
  d.analysis_energy_step = (p.e_max_wave_packet - p.e_min_wave_packet) / (p.n_energies - 1);

  d.steps_per_analysis = steps_per_analysis(p);
  d.analyses_per_average = analyses_per_average(p);
  // Rounded, not truncated as the counts above: dtFrames must be a whole multiple of dtAna (validate_outputs.cpp).
  d.analyses_per_frame = p.dt_frames ? static_cast<int>(std::lround(*p.dt_frames / p.dt_analysis)) : 0;

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

  d.analysis_time_step = d.steps_per_analysis * p.time_step;
  d.max_analysis_energy = units.hartree_to_ev(pi / d.analysis_time_step);

  const double k_max = d.kz + 2.0 * std::sqrt(d.wave_packet_alpha);
  const double k_min = d.kz - 2.0 * std::sqrt(d.wave_packet_alpha);
  d.confidence_e_max = 0.5 / d.total_mass * (k_max * k_max) * units.hartree_to_ev(1.0);
  d.confidence_e_min = 0.5 / d.total_mass * (k_min * k_min) * units.hartree_to_ev(1.0);
}

}  // namespace

Config build_config(const InputParameters& input) { return build_config(input, *make_potential(input)); }

Config build_config(const InputParameters& input, const Potential& potential) {
  Config config{input, {}};
  DerivedParameters& d = config.derived;

  // Fortran: InputNew.f90 SetupInput
  d.system = potential.system_parameters();
  // What identifies an SCM run in its restart files and analysis.json.
  if (const auto scm = potential.scm_potential()) d.scm = scm->get().identity();
  d.z_asymptotic_unswitched = d.system.z_asymptotic;
  d.z_asymptotic = input.use_switching ? input.z_max_switch : d.system.z_asymptotic;

  validate_input(config.input, d.system, d.z_asymptotic);
  process_input(config.input, d);
  miscellaneous_setup(config.input, d);
  validate_output_keys(config);
  return config;
}

// Fortran: InputNew.f90 ProcessInput, nanspo = INT(dtAnalysis/dtlan)
int steps_per_analysis(const InputParameters& p) { return static_cast<int>(p.dt_analysis / p.time_step); }

// Fortran: InputNew.f90 ProcessInput, nave = INT(dtAverages/dtAnalysis)
int analyses_per_average(const InputParameters& p) { return static_cast<int>(p.dt_averages / p.dt_analysis); }

}  // namespace spodvr

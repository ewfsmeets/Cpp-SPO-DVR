#pragma once

#include <optional>
#include <string>

#include "potential/potential.hpp"

namespace spodvr {

// The input parameters as read from the JSON input file, in input units (bohr, eV, amu, a.u. of
// time, degrees). Defaults are those of the Fortran code (globaldata.f90 SetupGlobalDataModule);
// keys that the Fortran code requires have no meaningful default here and must be present.
// The comment after each member gives the JSON key and the Fortran variable.
struct InputParameters {
  // --- general --------------------------------------------------------------------------
  std::string title = "SPO-DVR REACTIVE SCATTERING CALCULATION.";  // Title / title
  std::string system_id;                                           // SystemID / systemID
  bool new_calculation = true;                                     // NewCalc / newcalc
  bool z_dof = false;                                              // ZDoF
  bool r_dof = false;                                              // RDoF
  bool x_dof = false;                                              // XDoF
  bool y_dof = false;                                              // YDoF
  bool theta_dof = false;                                          // ThetaDoF
  bool phi_dof = false;                                            // PhiDoF

  // --- molecule ---------------------------------------------------------------------------
  double mass1 = 1.0078;      // Mass1 / amass1 (amu)
  double mass2 = 1.0078;      // Mass2 / amass2 (amu)
  bool homonuclear = true;    // HomoNuc / lhomo

  // --- grids ------------------------------------------------------------------------------
  double z_start = -1.0;      // ZStart / zmin
  double dz = 0.25;           // dZ / delz
  int nz = 56;                // nZ / nz
  int nz_specular = 80;       // nZSpecular / nzsp
  double r_start = 0.4;       // RStart / rmin
  double dr = 0.25;           // dR / delr
  int nr = 32;                // nR / nr
  int nx = 4;                 // nX / nx
  int ny = 4;                 // nY / ny
  int j_max = 6;              // JMax / jmaxw
  int mj_max = -1;            // MJMax / mjmaxw (a negative value is replaced by JMax on reading)

  int projection_operator = 2;  // ProjOp / nproj

  // --- propagation ------------------------------------------------------------------------
  double t_start = 0.0;       // tStart / starttime
  double t_end = 1000.0;      // tEnd / tottym
  double time_step = 10.0;    // dtStep / dtlan
  int n_steps = 100;          // nStep / nstep

  // --- analysis ---------------------------------------------------------------------------
  bool only_analysis = false;           // onlyAnalysis
  bool also_analysis = false;           // alsoAnalysis
  double dt_analysis = 50.0;            // dtAna / dtAnalysis
  bool dynamical_zsp_analysis = true;   // DynamZSpAna / dynamicalZSpAna
  double dt_averages = 50.0;            // dtAve / dtAverages
  double z_analysis = 7.0;              // Zana / ZAnalysis
  std::optional<int> j_analysis;        // JAna / jana (absent: determined by the program)
  std::optional<int> v_analysis;        // VAna / nviba (absent: determined by the program)
  std::optional<int> diffraction_order_analysis;  // DifOrdAna / mdifa (absent: all open orders)
  int n_energies = 41;                  // NrOfEs / nenrgy
  std::optional<double> z_analysis_specular;  // ZAnaSp / ZAnalysisSp

  // --- flux analysis ----------------------------------------------------------------------
  bool do_flux_analysis = false;        // doFluxAnalysis
  double r_flux_analysis = 0.0;         // rFluxAnalysis (required if doFluxAnalysis)
  double e_min_flux_analysis = 0.0;     // EMinFluxAna (required if doFluxAnalysis)
  double e_max_flux_analysis = 0.0;     // EMaxFluxAna (required if doFluxAnalysis)
  double delta_e_flux_analysis = 0.0;   // DeltaEFluxAna (required if doFluxAnalysis)

  // --- initial wave packet and state ------------------------------------------------------
  double e_min_wave_packet = 0.15;      // EMinWP (eV)
  double e_max_wave_packet = 0.55;      // EMaxWP (eV)
  double z0 = 10.0;                     // Z0 / z0
  int v_in = 0;                         // VIn / vin
  int j_in = 0;                         // JIn / jin
  int mj_in = 0;                        // MJIn / mjin
  double e_parallel = 0.0;              // EParallel / epar (eV)
  double phi_parallel = 0.0;            // PhiParallel / phipar (degrees)

  // --- potential --------------------------------------------------------------------------
  bool only_potential = false;          // Potential / onlyPotential
  bool impose_v_min = false;            // ImposeVMin / lvmin
  double v_min = -5.5;                  // VMin / vmin (eV)
  bool impose_v_max = false;            // ImposeVMax / lvmax
  double v_max = 10.0;                  // VMax / vmax (eV)
  bool use_switching = false;           // UseSwitching / lswitch
  double z_min_switch = 6.0;            // ZminSwitch / zmin_switch
  double z_max_switch = 8.0;            // ZmaxSwitch / zmax_switch

  // --- complex absorbing potentials (exactly one of energy / factor per CAP) --------------
  double z_opt_start = 7.0;             // ZOptStart / zsplmi
  double z_opt_end = 11.25;             // ZOptEnd / zsplma
  std::optional<double> z_opt_energy;   // ZOptEnergy (eV)
  std::optional<double> z_opt_factor;   // ZOptFactor / vabszb (hartree)
  double r_opt_start = 4.15;            // ROptStart / rsplmi
  double r_opt_end = 8.15;              // ROptEnd / rsplma
  std::optional<double> r_opt_energy;   // ROptEnergy (eV)
  std::optional<double> r_opt_factor;   // ROptFactor / vabsrb (hartree)
  double zsp_opt_start = 13.0;          // ZSpOptStart / ssplmi
  double zsp_opt_end = 18.0;            // ZSpOptEnd / ssplma
  std::optional<double> zsp_opt_energy;  // ZSpOptEnergy (eV)
  std::optional<double> zsp_opt_factor;  // ZSpOptFactor / vabssb (hartree)

  // --- output -----------------------------------------------------------------------------
  bool write_s_matrix = false;          // writeSMatrix / lwrit3
  bool write_probabilities = true;      // writeProbs / lwri25
  bool write_rd_probabilities = true;   // writeRDProbs / lwdiro
  bool do_dump = false;                 // DoDump
  bool full_s_matrix = false;           // FullSMatrix

  // --- tolerances -------------------------------------------------------------------------
  double tol_specular = 0.0;            // TolSp / tolsp (required: no Fortran default)
};

// Quantities derived from the input. Fortran: InputNew.f90 ProcessInput and the scalar part of
// Setup.f90 MiscellaneousSetup. Atomic units unless stated otherwise.
struct DerivedParameters {
  SystemParameters system;          // from the PES (Fortran: setupParametersH2SurfaceSystem)
  double z_asymptotic_unswitched;   // Fortran: zasymp_unswitched
  double z_asymptotic;              // = ZmaxSwitch if switching is used. Fortran: zasymp

  int run_mode;                     // 0 new, 1 restart, 2 analysis only, 6 potential only. Fortran: nread
  // MJMax used for the angular grid and rotational basis: MJMax < 1 means JMax (so MJMax = 0
  // also means JMax; a Fortran quirk that is kept). Fortran: setupdvr.f90 / transform.f90
  int mj_max_basis;
  int v_in_index;                   // vin + 1. Fortran: nvin

  double total_mass;                // Fortran: rmass / totmass (m_e)
  double reduced_mass;              // Fortran: rdmass (m_e)
  double rotational_constant;       // Fortran: rotcon (hartree; 59.304 cm^-1)

  double z_opt_factor;              // Fortran: vabszb
  double r_opt_factor;              // Fortran: vabsrb
  double zsp_opt_factor;            // Fortran: vabssb
  int optical_potential_power;      // Fortran: nabpow

  double dx;                        // Fortran: delx
  double dy;                        // Fortran: dely

  double translational_energy;      // eV. Fortran: enrgy
  double wave_packet_alpha;         // Fortran: alfwv
  double theta_incidence;           // rad. Fortran: thtai (after conversion in Setup)
  double phi_incidence;             // rad. Fortran: phii (after conversion in Setup)
  double analysis_energy_start;     // eV. Fortran: enrgi
  double analysis_energy_step;      // eV. Fortran: denrgy

  int steps_per_analysis;           // Fortran: nanspo
  int analyses_per_average;         // Fortran: nave
  bool fixed_z_analysis_specular;   // Fortran: FixedZAnalysisSpec

  // Setup.f90 MiscellaneousSetup
  double delta_kx;                  // reciprocal lattice spacing. Fortran: delkx
  double delta_ky;                  // Fortran: delky
  double k_total;                   // Fortran: rkt
  double kx;                        // initial momentum, reciprocal coordinates. Fortran: rkx
  double ky;                        // Fortran: rky
  double kz;                        // Fortran: rkz
  double unit_cell_area;            // Fortran: surfun
  double unit_cell_norm;            // Fortran: uclnrm
  double xy_volume_element;         // Fortran: facun
  double xyz_volume_element;        // Fortran: facn3d
  double v_min;                     // hartree. Fortran: vmin (after conversion)
  double v_max;                     // hartree. Fortran: vmax (after conversion)
  double analysis_time_step;        // Fortran: deltat
  double max_analysis_energy;       // eV. Fortran: etomax
  double confidence_e_min;          // eV. Fortran: ectmin
  double confidence_e_max;          // eV. Fortran: ectmax
};

// The complete, validated and processed configuration (immutable once built).
struct Config {
  InputParameters input;
  DerivedParameters derived;
};

}  // namespace spodvr

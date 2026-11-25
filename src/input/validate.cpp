// Input consistency checks. Fortran: InputNew.f90 ValidateInput
#include "input/validate.hpp"

#include <cmath>
#include <cstdlib>
#include <format>
#include <string>
#include <string_view>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fft_sizes.hpp"
#include "input/validation_helpers.hpp"
#include "potential/potential_factory.hpp"

namespace spodvr {

using validation::check;
using validation::nint;

// Fortran: InputNew.f90 ValidateInput (upstream r1727). The checks appear in the same order and with
// the Fortran's messages, with typos corrected and some made longer (e.g. NrOfEs, EMaxFluxAna,
// DeltaEFluxAna). C++-only checks lie in between (RestartDir, RestartFormat, CheckpointFormat,
// analysis D3/D4, SymmetryTolerance).
void validate_input(const InputParameters& p, const SystemParameters& system, double zasymp) {
  using namespace fortran_compare;
  const auto& units = units::legacy;

  // --- system ID and degrees of freedom -------------------------------------------------
  check(!potential_is_implemented(p.system_id), "SystemID", "Unknown system ID");
  check(!p.z_dof, "ZDoF", "Currently only 6D calculations are supported.");
  check(!p.r_dof, "RDoF", "Currently only 6D calculations are supported.");
  check(!p.x_dof, "XDoF", "Currently only 6D calculations are supported.");
  check(!p.y_dof, "YDoF", "Currently only 6D calculations are supported.");
  check(!p.theta_dof, "ThetaDoF", "Currently only 6D calculations are supported.");
  check(!p.phi_dof, "PhiDoF", "Currently only 6D calculations are supported.");

  // --- masses -----------------------------------------------------------------------------
  check(rle(p.mass1, 0.0), "Mass1", "Mass of atom 1 <= 0.");
  check(rle(p.mass2, 0.0), "Mass2", "Mass of atom 2 <= 0.");
  if (p.homonuclear) check(rne(p.mass1, p.mass2), "HomoNuc", "masses are different, but molecule is homo-nuclear.");
  if (!p.homonuclear) check(req(p.mass1, p.mass2), "HomoNuc", "masses are the same, but molecule is hetero-nuclear.");

  // --- scattering and specular grid in Z --------------------------------------------------
  check(p.nz < 1, "nZ", "nZ < 1.");
  check(!is_valid_fft_size(p.nz), "nZ", "nZ is not a valid FFT size.");
  check(rle(p.dz, 0.0), "dZ", "dZ <= 0.");
  // nproj /= 0 always holds here (ProjOp = 0 is rejected while reading).
  check(p.nz_specular < 1, "nZSpecular", "nZSpecular < 1.");
  check(!is_valid_fft_size(p.nz_specular), "nZSpecular", "nZSpecular is not a valid FFT size.");

  // --- grid in r --------------------------------------------------------------------------
  // The 5.24 rule (5.19 rejected only RStart < 0): the r grid must start above 0.
  check(rle(p.r_start, 0.0), "RStart", "RStart <= 0.");
  check(p.nr < 1, "nR", "nR < 1.");
  check(!is_valid_fft_size(p.nr), "nR", "nR is not a valid FFT size.");
  check(rle(p.dr, 0.0), "dR", "dR <= 0.");

  // --- grid in x and y --------------------------------------------------------------------
  check(p.nx < 1, "nX", "nX < 1.");
  check(p.ny < 1, "nY", "nY < 1.");
  check(!is_valid_fft_size(p.nx), "nX", "nX is not a valid FFT size.");
  check(!is_valid_fft_size(p.ny), "nY", "nY is not a valid FFT size.");
  if (seq(system.point_group, "C4v") || seq(system.point_group, "C6v") || seq(system.point_group, "C3v"))
    check(p.x_dof && p.y_dof && p.nx != p.ny, "nY", "nX /= nY.");

  // --- rotational basis -------------------------------------------------------------------
  check(p.j_in < 0, "JIn", "JIn < 0.");
  check(p.j_max < 0, "JMax", "JMax < 0.");
  check(std::abs(p.mj_max) > p.j_max, "MJMax", "|MJMax| > JMax.");
  check(std::abs(p.mj_in) > p.j_in, "MJIn", "|MJIn| > JIn.");
  check(p.j_in >= p.j_max, "JIn", "JIn >= JMax.");  // upstream r1726 (was >)
  check(std::abs(p.mj_in) > std::abs(p.mj_max), "MJIn", "|MJIn| > |MJMax|.");
  if (p.homonuclear)
    check(p.j_max % 2 != p.j_in % 2, "JMax",
          "the molecule is homo-nuclear but JMax is even and JIn is odd, or vice versa.");
  if (p.homonuclear)
    check(p.j_max % 2 != std::abs(p.mj_max) % 2, "MJMax",
          "the molecule is homo-nuclear, but JMax is even and |MJMax| is odd, or vice versa.");

  // --- projection operator ----------------------------------------------------------------
  check(p.projection_operator != 0 && p.projection_operator != 1 && p.projection_operator != 2, "ProjOp",
        "ProjOp must be 0, 1 or 2.");

  // --- split operator propagation ---------------------------------------------------------
  const double n_steps = p.n_steps;
  check(rlt(p.t_start, 0.0), "tStart", "tStart < 0.");
  check(p.new_calculation && rne(p.t_start, 0.0), "tStart", "This is a new calculation but tStart /= 0.");
  check(!p.new_calculation && req(p.t_start, 0.0), "tStart", "This is a restart but tStart = 0.");
  check(!p.new_calculation && p.restart_dir.empty(), "RestartDir",
        "This is a restart (NewCalc = false) but RestartDir is not given.");
  check(p.only_analysis && p.restart_dir.empty(), "RestartDir",
        "An analysis-only run (onlyAnalysis) needs RestartDir: the directory with the BK coefficient file.");
  check(p.restart_format != "native" && p.restart_format != "fortran", "RestartFormat",
        "RestartFormat must be \"native\" or \"fortran\".");
  check(p.checkpoint_format != "native" && p.checkpoint_format != "fortran" && p.checkpoint_format != "both",
        "CheckpointFormat", "CheckpointFormat must be \"native\", \"fortran\" or \"both\".");
  check(rlt(p.t_end, 0.0), "tEnd", "tEnd < 0.");
  check(rgt(p.t_start, p.t_end), "tStart", "tStart > tEnd.");
  check(rne(p.t_end, p.t_start + n_steps * p.time_step), "tEnd", "tEnd /= tStart + nStep*dtStep.");
  check(rle(p.time_step, 0.0), "dtStep", "dtStep <= 0.");
  check(p.n_steps < 0, "nStep", "nStep < 0.");

  // --- analysis ---------------------------------------------------------------------------
  check(rle(p.dt_analysis, 0.0), "dtAna", "dtAna <= 0.");
  check(rlt(p.dt_analysis, p.time_step), "dtAna", "dtAna < dtStep.");
  // analysis D4 (C++ only): the tolerant test above lets dtAna lie up to 1e-12 below dtStep, where the program's
  // dtAna / dtStep truncates to 0 time steps per analysis (the time loop would never end). The ratio itself is
  // compared, not its int: a huge ratio must not overflow the conversion.
  check(p.dt_analysis / p.time_step < 1.0, "dtAna",
        std::format("dtAna = {} lies within the 1e-12 tolerance below dtStep = {}, so dtAna / dtStep truncates to 0 "
                    "time steps per analysis; use dtAna >= dtStep.",
                    p.dt_analysis, p.time_step));
  check(rle(p.dt_averages, 0.0), "dtAve", "dtAve <= 0.");
  check(rlt(p.dt_averages, p.time_step), "dtAve", "dtAve < dtStep.");

  if (p.j_analysis) check(*p.j_analysis < 0, "JAna", "JAna < 0.");
  if (p.v_analysis) check(*p.v_analysis < 0, "VAna", "VAna < 0.");
  if (p.diffraction_order_analysis) check(*p.diffraction_order_analysis < 0, "DifOrdAna", "DifOrdAna < 0.");
  // NrOfEs < 2, as the Fortran since r1719 (analysis D1 is resolved; the -95 energy choice is gone).
  check(p.n_energies < 2, "NrOfEs",
        std::format("NrOfEs = {}: the analysis needs at least 2 energies, (EMaxWP-EMinWP)/(NrOfEs-1) apart; use "
                    "NrOfEs >= 2.",
                    p.n_energies));

  check(rne(nint(n_steps * p.time_step / p.dt_analysis) * p.dt_analysis, n_steps * p.time_step), "dtAna",
        "Propagation time must be an integer times dtAna.");
  // analysis D4: dtAna / dtStep truncates (a ratio just below k, a decimal such as 1.2 / 0.1 = 11.999999999999998, or
  // a ratio that is not whole) and nStep is not a whole number of the truncated steps: the time loop would stop
  // before tEnd. A truncated count above nStep leaves no whole analysis step (then only nStep = 0 is a multiple); at
  // most nStep, the conversion to int cannot overflow. (The count, not the ratio: dtAna 5e-13 above 5 dtStep with
  // nStep 5 is one whole analysis step.)
  const double analysis_ratio = p.dt_analysis / p.time_step;
  check(std::floor(analysis_ratio) > p.n_steps ? p.n_steps != 0 : p.n_steps % static_cast<int>(analysis_ratio) != 0,
        "dtAna",
        std::format("nStep = {} is not a multiple of the {} time steps per analysis that dtAna / dtStep = {} "
                    "truncates to; the run would stop before tEnd. Use a dtStep and dtAna whose ratio is exact in "
                    "binary, e.g. dtStep 0.125 or 0.25.",
                    p.n_steps, std::floor(analysis_ratio), analysis_ratio));
  check(rlt(p.dt_averages, p.dt_analysis), "dtAve", "dtAve < dtAna.");
  check(rne(nint(p.dt_averages / p.dt_analysis) * p.dt_analysis, p.dt_averages), "dtAve",
        "dtAve must be an integer times dtAna.");
  // analysis D4: dtAve up to 1e-12 below dtAna truncates to 0 analyses per average (a division by zero in the time
  // loop). dtAve just below m dtAna (m - 1 analyses per average) stays allowed: the averages carry their real times.
  check(p.dt_averages / p.dt_analysis < 1.0, "dtAve",
        std::format("dtAve = {} lies within the 1e-12 tolerance below dtAna = {}, so dtAve / dtAna truncates to 0 "
                    "analyses per average; use dtAve >= dtAna.",
                    p.dt_averages, p.dt_analysis));

  check(rlt(p.z_analysis, zasymp), "Zana", "ZAna < ZAsymptotic.");
  check(rne(nint((p.z_analysis - p.z_start) / p.dz) * p.dz, p.z_analysis - p.z_start), "Zana",
        "ZAna not on a grid point.");
  check(rlt(p.z_analysis, p.z_start), "Zana", "ZAna < ZStart.");
  check(rgt(p.z_analysis, p.z_start + (p.nz - 1) * p.dz), "Zana", "ZAnalysis is too big for the scattering grid in Z.");
  check(rgt(p.z_analysis, p.z_opt_start), "Zana", "ZAnalysis > ZOptStart.");

  if (p.j_analysis) check(*p.j_analysis > p.j_max, "JAna", "JAna > JMax.");
  if (p.v_analysis) check(*p.v_analysis > p.nr - 1, "VAna", "VAna > nR-1.");
  if (p.diffraction_order_analysis)
    check(*p.diffraction_order_analysis > p.nx + p.ny - 2, "DifOrdAna",
          "DifOrdAna > nX+nY-2 (the maximum diffraction order possible).");

  check(rlt(Constants::pi / p.dt_analysis, units.ev_to_hartree(p.e_max_wave_packet)), "dtAna",
        "dtAna is too big for fine analysis.");

  // --- flux analysis ----------------------------------------------------------------------
  if (p.do_flux_analysis) {
    check(rlt(p.r_flux_analysis, p.r_start), "rFluxAnalysis", "rFluxAnalysis < rmin");
    check(rgt(p.r_flux_analysis, p.r_start + (p.nr - 1) * p.dr), "rFluxAnalysis", "rFluxAnalysis > rmax");
    check(rne(nint((p.r_flux_analysis - p.r_start) / p.dr) * p.dr, p.r_flux_analysis - p.r_start), "rFluxAnalysis",
          "rFluxAnalysis not on a grid point.");
    check(rlt(p.e_min_flux_analysis, p.e_min_wave_packet), "EMinFluxAna", "EMinFluxAna < EMinWP");
    check(rgt(p.e_min_flux_analysis, p.e_max_wave_packet), "EMinFluxAna", "EMinFluxAna > EMaxWP");
    check(rlt(p.e_max_flux_analysis, p.e_min_wave_packet), "EMaxFluxAna", "EMaxFluxAna < EMinWP");
    // flux Q1 (fixed upstream in r1714, which checks the range as here): the range must end inside the
    // wave packet.
    check(rgt(p.e_max_flux_analysis, p.e_max_wave_packet), "EMaxFluxAna",
          std::format("EMaxFluxAna > EMaxWP ({} eV > {} eV): the flux analysis energies must lie inside the wave "
                      "packet's energy range; use EMaxFluxAna <= EMaxWP",
                      p.e_max_flux_analysis, p.e_max_wave_packet));
    check(rge(p.e_min_flux_analysis, p.e_max_flux_analysis), "EMinFluxAna", "EMinFluxAna >= EMaxFluxAna");
    // flux Q2 (fixed upstream in r1719): an exact <= 0, as the Fortran now checks it.
    check(p.delta_e_flux_analysis <= 0.0, "DeltaEFluxAna",
          std::format("DeltaEFluxAna <= 0 ({} eV): the flux analysis energy grid needs a positive spacing; use "
                      "DeltaEFluxAna > 0",
                      p.delta_e_flux_analysis));
    check(rgt(p.delta_e_flux_analysis, p.e_max_flux_analysis - p.e_min_flux_analysis), "DeltaEFluxAna",
          "DeltaEFluxAna > EMaxFluxAna-EMinFluxAna");
    check(rgt(p.r_flux_analysis, p.r_opt_start), "rFluxAnalysis", "rFluxAnalysis > ROptStart");
  }

  // --- initial wave packet ----------------------------------------------------------------
  check(rle(p.e_min_wave_packet, 0.0), "EMinWP", "EMinWP <= 0.");
  check(rle(p.e_max_wave_packet, 0.0), "EMaxWP", "EMaxWP <= 0.");
  check(rge(p.e_min_wave_packet, p.e_max_wave_packet), "EMinWP", "EMinWP >= EMaxWP.");
  check(rlt(p.z0, zasymp), "Z0", "Z0 < ZAsymptotic.");
  check(rgt(p.z0, p.zsp_opt_start), "Z0", "Z0 > ZSpOptStart.");

  // --- initial quantum state --------------------------------------------------------------
  check(p.v_in < 0, "VIn", "VIn < 0.");
  check(p.j_in < 0, "JIn", "JIn < 0.");
  check(std::abs(p.mj_in) > p.j_in, "MJIn", "|MJIn| > JIn.");
  check(p.v_in > p.nr - 1, "VIn", "VIn > nR-1.");
  // Unreachable since r1726: the rotational-basis block above already rejects JIn >= JMax. Kept
  // because the Fortran keeps it (InputNew.f90 ValidateInput, the initial-state block).
  check(p.j_in > p.j_max, "JIn", "JIn > JMax.");
  check(std::abs(p.mj_in) > p.mj_max, "MJIn", "MJIn > MJMax.");
  // Upstream r1726: VAna counts levels (v = 0 to VAna-1) and JAna is checked the same way.
  if (p.v_analysis) check(p.v_in >= *p.v_analysis, "VIn", "VIn >= VAna.");
  if (p.j_analysis) check(p.j_in >= *p.j_analysis, "JIn", "JIn >= JAna.");
  // analysis D3 (C++ only): the gas phase analyses at most nR-1 vibrational levels (GasPhase.f90
  // SetupAsympWF stops the open-level search there; with VAna given, VIn < VAna <= nR-1 already), so
  // VIn = nR-1 would leave the initial level out of the analysis, in the Fortran too.
  check(p.v_in == p.nr - 1, "VIn",
        std::format("VIn = nR-1 = {}: the analysis takes at most nR-1 vibrational levels (v = 0 to nR-2), so the "
                    "initial level would not be analysed; use VIn < nR-1.",
                    p.v_in));

  // --- potential --------------------------------------------------------------------------
  if (p.impose_v_min && p.impose_v_max) check(rgt(p.v_min, p.v_max), "VMin", "VMin > VMax.");
  check(p.symmetry_tolerance <= 0.0, "SymmetryTolerance", "SymmetryTolerance <= 0.");

  // --- optical potential in Z, scattering grid --------------------------------------------
  const double z_end = p.z_start + (p.nz - 1) * p.dz;
  check(rgt(p.z_opt_start, p.z_opt_end), "ZOptStart", "ZOptStart > ZOptEnd.");
  check(rgt(p.z_opt_end, z_end), "ZOptEnd", "ZOptEnd is too big for the scattering grid in Z.");
  // Duplicates of the ZAna block's "ZAnalysis > ZOptStart." (rgt), kept because the Fortran keeps them (InputNew.f90
  // ValidateInput, the optical-potential block). The second is unreachable (the same test). The first fires before
  // it only where ZOptStart - ZAna is exactly +-1e-12: rge is false there, so rlt (= !rge) is true while rgt is false
  // (at +1e-12 although ZOptStart lies above ZAna; MyOperators.f90 behaves the same).
  check(rlt(p.z_opt_start, p.z_analysis), "ZOptStart", "ZOptStart < ZAna.");
  check(rgt(p.z_analysis, p.z_opt_start), "Zana", "ZAna > ZOptStart.");

  // --- optical potential in r -------------------------------------------------------------
  check(rgt(p.r_opt_start, p.r_opt_end), "ROptStart", "ROptStart > ROptEnd.");
  check(rgt(p.r_opt_end, p.r_start + (p.nr - 1) * p.dr), "ROptEnd", "ROptEnd is too big for the grid in r.");

  // --- optical potential in Z, specular grid (used throughout for ProjOp = 2) -------------
  if (p.projection_operator == 2) {
    check(rgt(p.zsp_opt_start, p.zsp_opt_end), "ZSpOptStart", "ZSpOptStart > ZSpOptEnd.");
    check(rgt(p.zsp_opt_end, p.z_start + (p.nz_specular - 1) * p.dz), "ZSpOptEnd",
          "ZSpOptEnd is too big for the specular grid in Z.");
    check(rlt(p.zsp_opt_start, zasymp), "ZSpOptStart", "ZSpOptStart < ZAsymptotic.");
  }

  // --- parallel motion --------------------------------------------------------------------
  check(rlt(p.e_parallel, 0.0), "EParallel", "EParallel < 0.");
  check(rlt(p.phi_parallel, 0.0), "PhiParallel", "PhiParallel < 0.");
  check(rgt(p.phi_parallel, 360.0), "PhiParallel", "PhiParallel > 360 degrees.");

  // --- switching function -----------------------------------------------------------------
  if (p.use_switching) check(rge(p.z_min_switch, p.z_max_switch), "ZminSwitch", "ZminSwitch >= ZmaxSwitch.");

  // --- 2D cut of the PES (write2DPES) ------------------------------------------------------
  if (p.write_2d_pes) {
    check(rlt(p.x_2d_pes, 0.0), "x2DPES", "x2DPES < 0.");
    check(rgt(p.x_2d_pes, system.lattice_constant_x), "x2DPES", "x2DPES > alat.");
    check(rlt(p.y_2d_pes, 0.0), "y2DPES", "y2DPES < 0.");
    check(rgt(p.y_2d_pes, system.lattice_constant_y), "y2DPES", "y2DPES > blat.");
    check(rlt(p.theta_2d_pes, 0.0), "theta2DPES", "theta2DPES < 0.");
    check(rgt(p.theta_2d_pes, 360.0), "theta2DPES", "theta2DPES > 360.");
    check(rlt(p.phi_2d_pes, 0.0), "phi2DPES", "phi2DPES < 0.");
    check(rgt(p.phi_2d_pes, 360.0), "phi2DPES", "phi2DPES > 360.");
  }
}

}  // namespace spodvr

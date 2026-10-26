// Input consistency checks. Fortran: InputNew.f90 ValidateInput
#include "input/validate.hpp"

#include <cmath>
#include <cstdlib>
#include <string>

#include "core/constants.hpp"
#include "core/error.hpp"
#include "core/real_compare.hpp"
#include "core/unit_conversion.hpp"
#include "fft/fft_sizes.hpp"
#include "potential/potential_factory.hpp"

namespace spodvr {

namespace {

void check(bool violated, const char* parameter, const std::string& message) {
  if (violated) throw InputError(parameter, message);
}

// Fortran NINT: round half away from zero.
double nint(double x) { return static_cast<double>(std::lround(x)); }

}  // namespace

// Fortran: InputNew.f90 ValidateInput. The checks appear in the same order and with the same
// messages (a few message typos are corrected, see CLAUDE.md).
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
  if (p.homonuclear)
    check(rne(p.mass1, p.mass2), "HomoNuc", "masses are different, but molecule is homo-nuclear.");
  if (!p.homonuclear)
    check(req(p.mass1, p.mass2), "HomoNuc", "masses are the same, but molecule is hetero-nuclear.");

  // --- scattering and specular grid in Z --------------------------------------------------
  check(p.nz < 1, "nZ", "nZ < 1.");
  check(!is_valid_fft_size(p.nz), "nZ", "nZ is not a valid FFT size.");
  check(rle(p.dz, 0.0), "dZ", "dZ <= 0.");
  // nproj /= 0 always holds here (ProjOp = 0 is rejected while reading).
  check(p.nz_specular < 1, "nZSpecular", "nZSpecular < 1.");
  check(!is_valid_fft_size(p.nz_specular), "nZSpecular", "nZSpecular is not a valid FFT size.");

  // --- grid in r --------------------------------------------------------------------------
  check(rlt(p.r_start, 0.0), "RStart", "RStart < 0.");
  check(p.nr < 1, "nR", "nR < 1.");
  check(!is_valid_fft_size(p.nr), "nR", "nR is not a valid FFT size.");
  check(rle(p.dr, 0.0), "dR", "dR <= 0.");

  // --- grid in x and y --------------------------------------------------------------------
  check(p.nx < 1, "nX", "nX < 1.");
  check(p.ny < 1, "nY", "nY < 1.");
  check(!is_valid_fft_size(p.nx), "nX", "nX is not a valid FFT size.");
  check(!is_valid_fft_size(p.ny), "nY", "nY is not a valid FFT size.");
  if (seq(system.point_group, "C4v") || seq(system.point_group, "C6v") ||
      seq(system.point_group, "C3v"))
    check(p.x_dof && p.y_dof && p.nx != p.ny, "nY", "nX /= nY.");

  // --- rotational basis -------------------------------------------------------------------
  check(p.j_in < 0, "JIn", "JIn < 0.");
  check(p.j_max < 0, "JMax", "JMax < 0.");
  check(std::abs(p.mj_max) > p.j_max, "MJMax", "|MJMax| > JMax.");
  check(std::abs(p.mj_in) > p.j_in, "MJIn", "|MJIn| > JIn.");
  check(p.j_in > p.j_max, "JIn", "JIn > JMax.");
  check(std::abs(p.mj_in) > std::abs(p.mj_max), "MJIn", "|MJIn| > |MJMax|.");
  if (p.homonuclear)
    check(p.j_max % 2 != p.j_in % 2, "JMax",
          "the molecule is homo-nuclear but JMax is even and JIn is odd, or vice versa.");
  if (p.homonuclear)
    check(p.j_max % 2 != std::abs(p.mj_max) % 2, "MJMax",
          "the molecule is homo-nuclear, but JMax is even and |MJMax| is odd, or vice versa.");

  // --- projection operator ----------------------------------------------------------------
  check(p.projection_operator != 0 && p.projection_operator != 1 && p.projection_operator != 2,
        "ProjOp", "ProjOp must be 0, 1 or 2.");

  // --- split operator propagation ---------------------------------------------------------
  const double n_steps = p.n_steps;
  check(rlt(p.t_start, 0.0), "tStart", "tStart < 0.");
  check(p.new_calculation && rne(p.t_start, 0.0), "tStart",
        "This is a new calculation but tStart /= 0.");
  check(!p.new_calculation && req(p.t_start, 0.0), "tStart", "This is a restart but tStart = 0.");
  check(rlt(p.t_end, 0.0), "tEnd", "tEnd < 0.");
  check(rgt(p.t_start, p.t_end), "tStart", "tStart > tEnd.");
  check(rne(p.t_end, p.t_start + n_steps * p.time_step), "tEnd", "tEnd /= tStart + nStep*dtStep.");
  check(rle(p.time_step, 0.0), "dtStep", "dtStep <= 0.");
  check(p.n_steps < 0, "nStep", "nStep < 0.");

  // --- analysis ---------------------------------------------------------------------------
  check(rle(p.dt_analysis, 0.0), "dtAna", "dtAna <= 0.");
  check(rlt(p.dt_analysis, p.time_step), "dtAna", "dtAna < dtStep.");
  check(rle(p.dt_averages, 0.0), "dtAve", "dtAve <= 0.");
  check(rlt(p.dt_averages, p.time_step), "dtAve", "dtAve < dtStep.");

  if (p.j_analysis) check(*p.j_analysis < 0, "JAna", "JAna < 0.");
  if (p.v_analysis) check(*p.v_analysis < 0, "VAna", "VAna < 0.");
  if (p.diffraction_order_analysis)
    check(*p.diffraction_order_analysis < 0, "DifOrdAna", "DifOrdAna < 0.");
  if (p.n_energies != Constants::nen95_sentinel) check(p.n_energies < 0, "NrOfEs", "NrOfEs < 0.");

  check(rne(nint(n_steps * p.time_step / p.dt_analysis) * p.dt_analysis, n_steps * p.time_step),
        "dtAna", "Propagation time must be an integer times dtAna.");
  check(rlt(p.dt_averages, p.dt_analysis), "dtAve", "dtAve < dtAna.");
  check(rne(nint(p.dt_averages / p.dt_analysis) * p.dt_analysis, p.dt_averages), "dtAve",
        "dtAve must be an integer times dtAna.");

  check(rlt(p.z_analysis, zasymp), "Zana", "ZAna < ZAsymptotic.");
  check(rne(nint((p.z_analysis - p.z_start) / p.dz) * p.dz, p.z_analysis - p.z_start), "Zana",
        "ZAna not on a grid point.");
  check(rlt(p.z_analysis, p.z_start), "Zana", "ZAna < ZStart.");
  check(rgt(p.z_analysis, p.z_start + (p.nz - 1) * p.dz), "Zana",
        "ZAnalysis is too big for the scattering grid in Z.");
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
    check(rgt(p.r_flux_analysis, p.r_start + (p.nr - 1) * p.dr),
          "rFluxAnalysis", "rFluxAnalysis > rmax");
    check(rne(nint((p.r_flux_analysis - p.r_start) / p.dr) * p.dr, p.r_flux_analysis - p.r_start),
          "rFluxAnalysis", "rFluxAnalysis not on a grid point.");
    check(rlt(p.e_min_flux_analysis, p.e_min_wave_packet), "EMinFluxAna", "EMinFluxAna < EMinWP");
    // Fortran message says 'EMaxFluxAna < EMinWP' for this check.
    check(rlt(p.e_max_flux_analysis, p.e_max_wave_packet), "EMaxFluxAna", "EMaxFluxAna < EMaxWP");
    check(rge(p.e_min_flux_analysis, p.e_max_flux_analysis), "EMinFluxAna",
          "EMinFluxAna >= EMaxFluxAna");
    check(p.delta_e_flux_analysis < 0.0, "DeltaEFluxAna", "DeltaEFluxAna < 0");
    check(rgt(p.delta_e_flux_analysis, p.e_max_flux_analysis - p.e_min_flux_analysis),
          "DeltaEFluxAna", "DeltaEFluxAna > EMaxFluxAna-EMinFluxAna");
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
  check(p.j_in > p.j_max, "JIn", "JIn > JMax.");
  check(std::abs(p.mj_in) > p.mj_max, "MJIn", "MJIn > MJMax.");
  if (p.v_analysis) check(p.v_in > *p.v_analysis, "VIn", "VIn > VAna.");
  if (p.j_analysis) check(p.j_in > *p.j_analysis, "JIn", "JIn > JAna.");

  // --- potential --------------------------------------------------------------------------
  if (p.impose_v_min && p.impose_v_max) check(rgt(p.v_min, p.v_max), "VMin", "VMin > VMax.");

  // --- optical potential in Z, scattering grid --------------------------------------------
  const double z_end = p.z_start + (p.nz - 1) * p.dz;
  check(rgt(p.z_opt_start, p.z_opt_end), "ZOptStart", "ZOptStart > ZOptEnd.");
  check(rgt(p.z_opt_end, z_end), "ZOptEnd", "ZOptEnd is too big for the scattering grid in Z.");
  check(rlt(p.z_opt_start, p.z_analysis), "ZOptStart", "ZOptStart < ZAna.");
  check(rgt(p.z_analysis, p.z_opt_start), "Zana", "ZAna > ZOptStart.");

  // --- optical potential in r -------------------------------------------------------------
  check(rgt(p.r_opt_start, p.r_opt_end), "ROptStart", "ROptStart > ROptEnd.");
  check(rgt(p.r_opt_end, p.r_start + (p.nr - 1) * p.dr), "ROptEnd",
        "ROptEnd is too big for the grid in r.");

  // --- optical potential in Z, specular grid (used throughout for ProjOp = 2) -------------
  if (p.projection_operator == 2) {
    check(rgt(p.zsp_opt_start, p.zsp_opt_end), "ZSpOptStart", "ZSpOptStart > ZSpOptEnd.");
    check(rgt(p.zsp_opt_end, p.z_start + (p.nz_specular - 1) * p.dz),
          "ZSpOptEnd", "ZSpOptEnd is too big for the specular grid in Z.");
    check(rlt(p.zsp_opt_start, zasymp), "ZSpOptStart", "ZSpOptStart < ZAsymptotic.");
  }

  // --- parallel motion --------------------------------------------------------------------
  check(rlt(p.e_parallel, 0.0), "EParallel", "EParallel < 0.");
  check(rlt(p.phi_parallel, 0.0), "PhiParallel", "PhiParallel < 0.");
  check(rgt(p.phi_parallel, 360.0), "PhiParallel", "PhiParallel > 360 degrees.");

  // --- switching function -----------------------------------------------------------------
  if (p.use_switching)
    check(rge(p.z_min_switch, p.z_max_switch), "ZminSwitch", "ZminSwitch >= ZmaxSwitch.");
}

}  // namespace spodvr

// The SCM distorted-surface PES (SystemID CRPDistrtSurface): the CRP potential plus ΔV, its set-up
// checks and its system parameters. Fortran: DistortedSurfacePES.dir/DistortedSurfaceWrapper.f90
// (GetCRPDistSurfParams, GetCRPDistSurfPotential), SCM_Wrapper.f90 (SetSCMCRPParameters)
#include "potential/scm/scm_crp_potential.hpp"

#include <format>

#include "core/error.hpp"
#include "core/unit_conversion.hpp"

namespace spodvr {

namespace {

// s = a_T/a₀, the slab's lattice constant over the CRP's, if the slab is rescaled to the CRP
// lattice. a_T must be positive then (SCM D6; the Fortran does not check it).
// Fortran: SCM_Wrapper.f90 SetSCMCRPParameters (LatExpansion)
std::optional<double> rescaling_factor_of(const ScmInput& input, const ScmSlab& slab, const SystemParameters& crp) {
  if (!input.rescale) return std::nullopt;
  const double slab_lattice_constant = slab.header.lattice_constant;
  if (!(slab_lattice_constant > 0.0))
    throw InputError("PosFilepath", std::format("{}: the lattice constant a_T on line 3 is {} bohr; with "
                                                "EnableSurfExpan it rescales the lattice (s = a_T/a₀) and must be "
                                                "positive",
                                                input.position_file, slab_lattice_constant));
  return slab_lattice_constant / crp.lattice_constant_x;
}

// SCM D14: the slab file holds an FCC111 surface (the reader requires it), so the CRP data set must
// describe a hexagonal one. The Fortran does not check this.
void require_hexagonal_data_set(const CrpPotential& crp) {
  const std::string& symmetry = crp.summary().symmetry;
  if (symmetry == "C3v" || symmetry == "C6v") return;
  throw InputError("CRPlibDataSet",
                   std::format("the data set is {}, but the SCM slab (PosFilepath) is an FCC111 surface: SystemID "
                               "CRPDistrtSurface needs a C3v or C6v data set",
                               symmetry));
}

// SCM D14 and D9, before ΔV is set up.
ScmLatticeCheck checked_lattice(const CrpPotential& crp, const ScmSlab& slab) {
  require_hexagonal_data_set(crp);
  ScmLatticeCheck check = check_slab_lattice(slab, crp.system_parameters(), crp.summary());
  require_slab_lattice(check, slab, crp.system_parameters());
  return check;
}

// The CRP's system parameters, with the slab's lattice constant a_T for both lattice vectors when
// rescaling. Fortran: DistortedSurfaceWrapper.f90 GetCRPDistSurfParams (xlat_const, ylat_const)
SystemParameters system_parameters_of(const SystemParameters& crp, const ScmSlab& slab, bool rescaling) {
  SystemParameters parameters = crp;
  if (rescaling) {
    parameters.lattice_constant_x = slab.header.lattice_constant;
    parameters.lattice_constant_y = slab.header.lattice_constant;
  }
  return parameters;
}

}  // namespace

ScmCrpPotential::ScmCrpPotential(const std::filesystem::path& data_set, const ScmInput& input)
    : crp_(data_set),
      slab_(read_scm_slab(input.position_file, input.configuration, input.cutoff)),
      rescaling_factor_(rescaling_factor_of(input, slab_, crp_.system_parameters())),
      crp_scale_(rescaling_factor_.value_or(1.0)),
      lattice_check_(checked_lattice(crp_, slab_)),
      correction_(read_scm_parameters(input.parameter_file), slab_, rescaling_factor_,
                  crp_.system_parameters().skewing_angle),
      parameters_(system_parameters_of(crp_.system_parameters(), slab_, rescaling_factor_.has_value())) {}

double ScmCrpPotential::value(double u, double v, double z, double r, double theta, double phi, int iz, int ir) const {
  const CartesianVector xy = correction_.to_cartesian(u, v);
  const double v_crp_ev = crp_.value(u / crp_scale_, v / crp_scale_, z, r, theta, phi, iz, ir);
  // ΔV in eV with Constants_eV2Ha (a division), added to the CRP's eV (SCM D4).
  return v_crp_ev + units::precise.hartree_to_ev(correction_.value_hartree(xy.x, xy.y, z, r, theta, phi));
}

}  // namespace spodvr

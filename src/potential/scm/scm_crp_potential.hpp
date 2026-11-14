#pragma once
// The SCM distorted-surface PES (SystemID CRPDistrtSurface): a CRP potential plus the SCM correction
// ΔV of one slab snapshot. Fortran: DistortedSurfacePES.dir/DistortedSurfaceWrapper.f90
// (GetCRPDistSurfParams, GetCRPDistSurfPotential), H2SurfaceSystem.f90 (CRPDistrtSurface branch)

#include <filesystem>
#include <optional>
#include <string>

#include "input/config.hpp"
#include "potential/crp_potential.hpp"
#include "potential/potential.hpp"
#include "potential/scm/scm_correction.hpp"
#include "potential/scm/scm_diagnostics.hpp"
#include "potential/scm/scm_files.hpp"

namespace spodvr {

// V(u, v, …) = CRP(u/s, v/s, …) + ΔV(x, y, …)/Constants_eV2Ha in eV, (x, y) the Cartesian position of
// (u, v); s = a_T/a₀ with rescaling (EnableSurfExpan), else 1 (SCM D4). Point group, z_asymptotic
// and skew angle are the CRP's. Evaluation is const and thread-safe.
class ScmCrpPotential final : public Potential {
 public:
  // Reads the data set and the SCM files (relative to the working directory) and checks them (SCM
  // D6, D9, D14); throws InputError naming the key. Fortran: InitialiseDistPESModule,
  // GetCRPDistSurfParams (DistortedSurfaceWrapper.f90), SetSCMCRPParameters, SetupSCMModule
  ScmCrpPotential(const std::filesystem::path& data_set, const ScmInput& input);

  // V in eV in the Fortran's order: (x, y), the CRP at (u/s, v/s), then ΔV added; r is an r grid
  // point (see scm_pair_terms). Fortran: DistortedSurfaceWrapper.f90 GetCRPDistSurfPotential
  [[nodiscard]] double value(double u, double v, double z, double r, double theta, double phi, int iz,
                             int ir) const override;
  // The CRP's, with both lattice constants a_T when rescaling. Fortran: GetCRPDistSurfParams
  [[nodiscard]] const SystemParameters& system_parameters() const override { return parameters_; }
  [[nodiscard]] std::string name() const override { return "CRPDistrtSurface"; }

  // The CRP potential of the data set (its own, unrescaled system parameters and summary).
  [[nodiscard]] const CrpPotential& crp() const { return crp_; }
  // The slab snapshot as read (layers 1 and 2 inside the cut-off).
  [[nodiscard]] const ScmSlab& slab() const { return slab_; }
  // s = a_T/a₀ with rescaling, else empty. Fortran: LatExpansion (thermExpansion)
  [[nodiscard]] std::optional<double> rescaling_factor() const { return rescaling_factor_; }
  // SCM D9: how far the slab's ideal atoms lie from the CRP lattice.
  [[nodiscard]] const ScmLatticeCheck& lattice_check() const { return lattice_check_; }

 private:
  CrpPotential crp_;
  ScmSlab slab_;
  std::optional<double> rescaling_factor_;
  double crp_scale_;  // s with rescaling, else 1 (an exact division). Fortran: thermExpansion
  ScmLatticeCheck lattice_check_;
  ScmCorrection correction_;
  SystemParameters parameters_;
};

}  // namespace spodvr

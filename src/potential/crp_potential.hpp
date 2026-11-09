#pragma once

#include <filesystem>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "potential/potential.hpp"

namespace crplib {
class QdWrapper;
}

namespace spodvr {

// The CRP data set as the setup report shows it (lengths in bohr, energies in eV).
struct CrpDataSetSummary {
  std::string directory;
  std::optional<std::string> description;  // first line of CRPData.info
  std::string symmetry;                    // C3v, C4v or C6v
  std::string file_units;                  // units of the data-set files, e.g. "angstrom, eV"
  int surface_layers = 0;
  int atom_cuts = 0;
  int molecule_cuts = 0;                   // in CRPData.6D
  int molecule_cuts_used = 0;              // the number the interpolation uses
  bool cuts_hold_interpolation = false;    // CutsInterpolation: the cuts hold I, not V
  double r_min = 0.0, r_max = 0.0;         // range covered by every molecular cut. Fortran: Rmin, Rmax
  double z_min = 0.0;                      // Fortran: Zmin
  double zso_min = 0.0, zso_max = 0.0;     // switching to the 2D gas-phase potential. Fortran: ZsoMin, ZsoMax
  double z_asymptotic = 0.0;               // Fortran: Zasy
  double energy_shift = 0.0;               // V is shifted by -VH2(r_last) (eV). Fortran: VH2%F(NX)
  int cluster_atoms = 0;                   // atoms of the repulsive sum. Fortran: MURII, MURO
  int cluster_atoms_within_reach = 0;      // of these, the ones the corrugation function can reach
  std::vector<std::string> warnings;       // conditions the Fortran accepts that are probably not intended
};

// How the grid points relate to the range of the data (see crplib::RangeReport).
struct CrpRangeReport {
  int r_points = 0;
  int r_clamped = 0;            // evaluated at Rmin or Rmax
  int z_points = 0;
  int z_clamped = 0;            // evaluated at Zmin
  long long zr_points = 0;
  long long zr_uncovered = 0;   // V = -VH2(r_last): the Fortran's silent V = 0 branch
};

// The CRP PES of a data set (SystemID "CRPlibrary"), evaluated by the port of the CRP library
// (src/crplib/, docs/crp_library.md). x, y are skewed lattice coordinates (bohr), as for every
// Potential. r and Z are clamped into the range of the data. This is the only file that sees
// the library. Fortran: H2SurfaceSystem.f90 (CRPlibrary branch), CRP_QD_Wrapper.f90
class CrpPotential final : public Potential {
 public:
  // Reads and checks the data set (directory relative to the working directory). Throws
  // InputError("CRPlibDataSet", ...) if it cannot be used.
  explicit CrpPotential(const std::filesystem::path& data_set);
  ~CrpPotential() override;
  CrpPotential(const CrpPotential&) = delete;
  CrpPotential& operator=(const CrpPotential&) = delete;

  // V in eV. Throws InputError("CRPlibDataSet", ...) where the Fortran library stops (a geometry
  // outside the data, e.g. an atom below the atom-surface cuts). Fortran: Get_CRP_QD_Wrapper_Potential_6D
  [[nodiscard]] double value(double x, double y, double z, double r, double theta, double phi, int iz,
                             int ir) const override;
  // Fortran: Get_CRP_QD_Wrapper_Params
  [[nodiscard]] const SystemParameters& system_parameters() const override { return parameters_; }
  [[nodiscard]] std::string name() const override { return "CRPlibrary"; }

  [[nodiscard]] const CrpDataSetSummary& summary() const { return summary_; }
  // The clamped grid points and the (Z, r) points without data (z: the Z points where the 6D
  // potential is evaluated; r: the whole r grid).
  [[nodiscard]] CrpRangeReport range_report(std::span<const double> z, std::span<const double> r) const;

 private:
  std::unique_ptr<const crplib::QdWrapper> wrapper_;
  SystemParameters parameters_;
  CrpDataSetSummary summary_;
};

}  // namespace spodvr

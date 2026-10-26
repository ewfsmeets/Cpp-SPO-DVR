#pragma once

#include <optional>

#include "grids/momentum_grid.hpp"
#include "grids/uniform_grid.hpp"
#include "input/config.hpp"

namespace spodvr {

// Inclusive 0-based index range [first, last] of grid points inside an optical potential.
struct IndexRange {
  int first;
  int last;
  [[nodiscard]] int size() const { return last - first + 1; }
};

// The coordinate and momentum grids in Z, r, x and y plus the grid indices derived from them.
// All indices are 0-based (Fortran index - 1).
// Fortran: Setup.f90 MiscellaneousSetup
struct SpatialGrids {
  UniformGrid z_scattering;   // nZ points. Fortran: zgrid(1:nz)
  UniformGrid z_specular;     // nZSpecular points. Fortran: zgrid(1:nzsp)
  UniformGrid r;              // Fortran: rgrid
  UniformGrid x;              // Fortran: xgrid
  UniformGrid y;              // Fortran: ygrid

  MomentumGrid kz_scattering;  // Fortran: delkz, zk, ekinz
  MomentumGrid kz_specular;    // Fortran: delksp, zksp, ensp (= ezsp)
  MomentumGrid kr;             // Fortran: delkr, rkg, ekinr

  int iz_analysis;             // Fortran: izanly - 1
  double z_analysis;           // ZAna moved onto the grid. Fortran: ZAnalysis after Setup
  int nz_potential;            // number of Z points with Z < zasymp. Fortran: nzpot
  IndexRange z_optical;        // Fortran: nzchsp .. nzchmx
  std::optional<IndexRange> z_specular_optical;  // ProjOp = 2 only. Fortran: nschsp .. nschmx
  IndexRange r_optical;        // Fortran: nrchsp .. nrchmx
  std::optional<int> ir_flux_analysis;  // Fortran: irFluxAnalysis - 1
};

// Builds the grids and performs the grid-related input checks of Setup.f90 (throws InputError).
[[nodiscard]] SpatialGrids make_spatial_grids(const Config& config);

}  // namespace spodvr

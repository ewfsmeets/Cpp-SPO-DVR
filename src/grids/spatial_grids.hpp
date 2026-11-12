#pragma once

#include <optional>
#include <string>
#include <vector>

#include "core/index_range.hpp"
#include "grids/momentum_grid.hpp"
#include "grids/uniform_grid.hpp"
#include "input/config.hpp"

namespace spodvr {

// The end of a grid an absorbing potential sits at. High: zero at its start, strongest at its
// end (ZOpt, ROpt, ZSpOpt). Low: the mirror image, strongest at its start (ZOptLow, ZSpOptLow).
enum class AbsorberSide { High, Low };

// A grid point that this program and the Fortran put on different sides of an absorber's strongest
// end (the grid points differ by an ulp: docs/fortran_deviations.md). There one code multiplies
// it by exp(-1.5 A dt) every step and the other by 1; the set-up prints a warning.
struct SplitGridPoint {
  std::string key;       // input key of the strongest end (ZOptEnd, ROptEnd, ZSpOptEnd, …LowStart)
  double end;            // its value (bohr)
  int index;             // 0-based grid index
  double point;          // the point here (plain double)
  double fortran_point;  // the Fortran's point (quad precision, rounded once)
  bool inside_here;      // inside the absorber here (and outside in the Fortran, or vice versa)
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
  IndexRange z_optical;        // grid points inside the optical potential. Fortran: nzchsp .. nzchmx
  std::optional<IndexRange> z_specular_optical;  // ProjOp = 2 only. Fortran: nschsp .. nschmx
  IndexRange r_optical;        // Fortran: nrchsp .. nrchmx
  std::optional<int> ir_flux_analysis;  // Fortran: irFluxAnalysis - 1
  // Points inside the low-Z absorbers, if on (Fortran: nzchsplow .. nzchmxlow, nschsplow ..
  // nschmxlow; the Fortran prints nzchmilow = nzchsplow - 1, the last point below Start).
  std::optional<IndexRange> z_optical_low;
  std::optional<IndexRange> z_specular_optical_low;
  std::vector<SplitGridPoint> split_points;  // of every absorber that is on
};

// Builds the grids and performs the grid-related input checks of Setup.f90 (throws InputError).
[[nodiscard]] SpatialGrids make_spatial_grids(const Config& config);

}  // namespace spodvr

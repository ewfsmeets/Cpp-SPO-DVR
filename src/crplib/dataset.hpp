#pragma once
// A CRP data set: the content of the data-set directory, exactly as written in its files (values in
// the data set's own units; the PES constructors convert them as the Fortran does while reading).
// Fortran: the reading parts of CRP_Wrapper.f90 (SetUp_CRP_*_PES), CRP_C*v_AtomSurface.f90 and
// CRP_C*v_MoleculeSurface.f90 (Read_PES_Info, Read_*_Data). Layout of the files:
// docs/crp_library.md "Data-set format".

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "crplib/units.hpp"

namespace crplib {

// Surface symmetry of a data set (first line of CRPData).
enum class Symmetry { c3v, c4v, c6v };

// "C3v", "C4v", "C6v"
const char* symmetry_name(Symmetry symmetry);

// The number of cuts the interpolation uses: atom-surface sites (exact for C3v/C6v) and molecular
// configurations (the minimum; further cuts are read but unused).
int atom_cut_count(Symmetry symmetry);
int molecule_cut_count(Symmetry symmetry);

// A tabulated function of one variable (a two-column file). Fortran: Spline1D X, F
struct Table1D {
  std::vector<double> x;
  std::vector<double> f;
};

// A 1D cut through the atom-surface PES above one site. Fortran: PESCut3D
struct AtomCut {
  double x = 0.0;         // site, in units of the lattice constant, Cartesian (CRPData.3D)
  double y = 0.0;
  std::string file_name;  // in Atom/
  double slope_first = 0.0;  // dV/dz at the first and last z (energy/length). Fortran: dz1, dz2
  double slope_last = 0.0;
  Table1D potential;  // z -> V
};

// A 2D cut (r, Z) through the molecule-surface PES for one site and orientation. Fortran: PESCut6D
struct MoleculeCut {
  double x = 0.0;         // centre of mass, lattice units: skewed (u, v) for C3v/C6v, Cartesian for C4v
  double y = 0.0;
  double theta = 0.0;     // orientation (degrees); used only when CutsInterpolation is false
  double phi = 0.0;
  std::string file_name;  // in Molecule/
  std::vector<double> r;  // grid in r (fastest in the file) and in Z
  std::vector<double> z;
  std::vector<double> v;  // v[ir + r.size() * iz]
};

// Optional description of the data set (CRPData.info, three lines).
struct DataSetInfo {
  std::string short_info, long_info, reference;
};

struct DataSet {
  // CRPData
  Symmetry symmetry = Symmetry::c6v;
  Units units = internal_units;               // units of every file of the data set
  double lattice_constant = 0.0;              // nearest-neighbour distance on the surface
  std::vector<double> interlayer_spacings;    // NumLayers - 1 values
  double displacement_1 = 0.0;                // optional; Fortran: Disp1 (0 when absent)
  double displacement_2 = 0.0;                // Fortran: Disp2
  std::optional<DataSetInfo> info;            // CRPData.info

  // Atom-surface PES: Atom/intrep.dat (corrugation function and its end slopes) and the cuts.
  Table1D corrugation;
  double corrugation_slope_first = 0.0;  // Fortran: dz1, dz2 of intrep.dat
  double corrugation_slope_last = 0.0;
  std::vector<AtomCut> atom_cuts;

  // Molecule-surface PES: CRPData.6D, Molecule/extrapol.dat, Molecule/VH2.dat and the cuts.
  bool cuts_interpolation = false;  // true: the cuts hold the interpolation function, not V
  double zso_min = 0.0;             // switching to the 2D gas-phase potential between these Z
  double zso_max = 0.0;
  Table1D z_dependence;   // Z -> V of the gas-phase potential (extrapol.dat). Fortran: Zdep
  Table1D h2_potential;   // r -> V of the free molecule (VH2.dat). Fortran: Rdep
  std::vector<MoleculeCut> molecule_cuts;

  [[nodiscard]] int num_layers() const { return static_cast<int>(interlayer_spacings.size()) + 1; }
};

// Reads and validates the data set in `directory` (docs/crp_library.md "Validation"); throws
// DataError naming the file and the problem.
DataSet read_data_set(const std::filesystem::path& directory);

// Writes a data set in the format read_data_set (and the Fortran library) reads; every value is
// written with enough digits to be read back exactly. Creates the directory.
void write_data_set(const DataSet& data, const std::filesystem::path& directory);

// Throws DataError if the data set cannot be used (same checks as read_data_set).
void validate_data_set(const DataSet& data);

// Conditions the Fortran accepts but that are probably not intended (e.g. unused extra cuts).
std::vector<std::string> data_set_warnings(const DataSet& data);

}  // namespace crplib

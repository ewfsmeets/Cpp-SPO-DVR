#pragma once
// The two data files of the SCM distorted-surface correction: the parameter file (params.dat) and
// the slab position file (snapshots of a distorted metal slab). Fortran:
// DistortedSurfacePES.dir/SCM_IO.f90 (Read_SCM_Parameters, ReadSlabFileParams, ReadSCMSlab).
// Example files: potentials/H2Cu111SRP_C3v.dir/params.dat and examplePosfile.dat.

#include <array>
#include <filesystem>
#include <string>
#include <vector>

namespace spodvr {

// The pair-function parameters of one surface layer: one list-directed record of 18 values in
// params.dat (version 2), as written (atomic units). At bond length r the eight parameters are
// slope·r_c + intercept, with r_c = r clamped to [P17, P18].
// Fortran: SCM_Types.f90 SCM_Parameters_Rydberg%P(1:18)
struct ScmLayerParameters {
  // P1–P8: reference distance, decay, the four polynomial coefficients, switch start and end.
  std::array<double, 8> intercept{};
  std::array<double, 8> slope{};  // P9–P16: the change of P1–P8 per unit r
  double r_clamp_min = 0.0;       // P17
  double r_clamp_max = 0.0;       // P18 (larger than P17)
};

// The content of params.dat. Fortran: SCM_Types.f90 SCM_Data%ParameterVersion, %RydbergParams(1:2)
struct ScmParameters {
  int version = 0;                             // n of "VERSION n" on line 1; only 2 exists
  std::array<ScmLayerParameters, 2> layers{};  // [0] the top layer, [1] the layer below it
};

// A position in bohr: Cartesian, the ideal top layer at z = 0, the surface normal along +z.
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

// One slab atom: its position in the snapshot and in the ideal lattice, and its layer.
// Fortran: SCM_Types.f90 SCM_Data%VectorsDisplaced, %VectorsIdeal, %AtomType
struct ScmAtom {
  Vec3 distorted;
  Vec3 ideal;
  int layer = 0;  // 1 = the top layer, 2 = the layer below it
};

// The three header lines of a position file. Fortran: the SCM_Types.f90 SCM_Data fields named
// below, read by SCM_IO.f90 ReadSlabFileParams
struct ScmSlabHeader {
  std::string tag;                  // free text naming the file. Fortran: SlabFileTag
  std::string structure;            // the surface; only FCC111 exists. Fortran: SlabSurfStruct
  int nx = 0, ny = 0, nz = 0;       // atoms per row, rows per layer, layers. Fortran: SlabAtomSize
  int configuration_count = 0;      // snapshots in the file. Fortran: NumSlabConfigs
  int version = 0;                  // 3 or 2 (1 is not supported). Fortran: PosFileVersion
  int max_surface_atoms = 0;        // atoms per snapshot in version 2. Fortran: MaxNumSurfAtoms
  double lattice_constant = 0.0;    // a_T (bohr), the slab's own one. Fortran: DistLatticeConstant
};

// One snapshot (configuration) of a position file, reduced to the atoms the correction uses.
// Fortran: SCM_Types.f90 SCM_Data after SCM_IO.f90 ReadSCMSlab
struct ScmSlab {
  ScmSlabHeader header;
  int configuration = 0;  // the snapshot read, 1-based. Fortran: UsedSlabNum (input PosFileSlabNum)
  double cutoff = 0.0;    // bohr; applied in version 3 only. Fortran: CutOff (input SCMCutoff)
  // Layers 1 and 2 only, in file order. Version 3: the atoms whose distorted position lies inside
  // the cut-off sphere around the origin; version 2: the first max_surface_atoms of the snapshot.
  std::vector<ScmAtom> atoms;
};

// Reads params.dat: line 1 "VERSION 2" (tag case-sensitive), then one list-directed READ of the 18
// values per layer (a record may span lines; the rest of its last line is ignored). Throws
// InputError("SCMParFilepath", ...) naming the file, the line or the parameter, if the file cannot
// be read, has another version, or has P17 ≥ P18. Fortran: SCM_IO.f90 Read_SCM_Parameters
ScmParameters read_scm_parameters(const std::filesystem::path& file);

// Reads the header and snapshot `configuration` (1-based) of a position file of version 3 or 2:
// the earlier snapshots are skipped record by record, as the Fortran does. Throws InputError,
// naming the input key (PosFilepath, PosFileSlabNum or SCMCutoff), the file and the line or
// parameter, if the file cannot be read or does not hold a usable snapshot (structure not FCC111,
// fewer than 2 layers, a layer index ≤ 0, no atom of layers 1–2 kept, ...), if `configuration` is
// out of range, or if `cutoff` (bohr) is not positive. Fortran: SCM_IO.f90 GetSurfaceParameters,
// ReadSCMSlab
ScmSlab read_scm_slab(const std::filesystem::path& file, int configuration, double cutoff);

}  // namespace spodvr

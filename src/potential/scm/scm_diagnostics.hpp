#pragma once
// Checks of an SCM slab snapshot against its CRP data set, none of which the Fortran has (SCM D9):
// check_* returns the result as data for the set-up report, require_* throws an InputError naming
// the input key when the result is unacceptable.

#include "potential/crp_potential.hpp"
#include "potential/potential.hpp"
#include "potential/scm/scm_files.hpp"

namespace spodvr {

// Ideal atoms farther than this from the CRP lattice give a warning (bohr). The SRP48 example slab
// lies within 1.2e-6 bohr: its bohr constant differs from the CRP library's by 7.2e-8 relative, and
// its farthest kept atom is 15.8 bohr from the origin (SCM D9).
inline constexpr double scm_lattice_warning_distance = 1e-5;
// Ideal atoms farther than this from the CRP lattice are an error (bohr): ΔV would subtract the
// pair energies of another lattice than the one the CRP potential describes.
inline constexpr double scm_lattice_error_distance = 1e-3;

// How far the slab's ideal atom positions lie from the sites of the CRP lattice (SCM D9). The
// deviation of an atom is its distance (bohr) from the nearest site of its layer.
struct ScmLatticeCheck {
  int top_layer_atoms = 0;              // layer-1 atoms checked
  int second_layer_atoms = 0;           // layer-2 atoms checked (none for a C6v data set)
  double top_layer_deviation = 0.0;     // bohr, the largest of layer 1
  double second_layer_deviation = 0.0;  // bohr, the largest of layer 2
  int worst_atom = -1;                  // index in ScmSlab::atoms of the largest deviation; −1: none checked

  // The largest deviation of both layers (NaN if one is NaN).
  [[nodiscard]] double max_deviation() const;
  // True if the largest deviation is above scm_lattice_warning_distance (or NaN). Out of line, so
  // that the NaN test is compiled with strict floating point in every build.
  [[nodiscard]] bool warning() const;
};

// SCM D9: the ideal layer 1 on the sites of the CRP lattice (a₀, γ of `crp`, unrescaled) at z = 0;
// with a C3v data set layer 2 at z = −d₁ under the data set's second-layer site (not checked for
// C6v, a one-layer PES).
[[nodiscard]] ScmLatticeCheck check_slab_lattice(const ScmSlab& slab, const SystemParameters& crp,
                                                 const CrpDataSetSummary& data_set);

// Throws InputError("PosFilepath", ...), naming the worst atom, if an ideal atom lies farther than
// scm_lattice_error_distance from its site (or the deviation is NaN).
void require_slab_lattice(const ScmLatticeCheck& check, const ScmSlab& slab, const SystemParameters& crp);

}  // namespace spodvr

#pragma once
// The write2DPES output: pes_cut_2d.dat (V(r, Z) at one fixed x, y, theta, phi) and its set-up lines.
// Fortran: Potential.f90 Write2DPESToFile (2dpes.dat), Main.f90 (the lines it prints)

#include <filesystem>
#include <ostream>
#include <string_view>

#include "potential/pes_cut.hpp"
#include "report/format.hpp"

namespace spodvr {

// The file name of the cut in the output directory (the Fortran writes 2dpes.dat).
inline constexpr std::string_view pes_cut_file_name = "pes_cut_2d.dat";

// Writes the cut: '#' lines with x, y (bohr), theta, phi (degrees and radians) and the columns, then one line
// "r Z V(hartree) V(eV)" (%.16e) per point, Z fastest, with a blank line after each r block (the Fortran's
// order; V(hartree) = V(eV) / 27.2114, the legacy factor).
void write_pes_cut_file(const std::filesystem::path& path, const PesCut& cut);

// The set-up section of a write2DPES run: the cut's position, grids and file.
void print_pes_cut(std::ostream& out, const Style& style, const PesCut& cut, const std::filesystem::path& file);

}  // namespace spodvr

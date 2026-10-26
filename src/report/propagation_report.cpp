// Printed sections and data files of the wave-packet propagation.
#include "report/propagation_report.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace spodvr {

void print_propagation_setup(std::ostream& out, const Style& style, const SpatialGrids& grids,
                             const WavePacketPropagation& propagation) {
  const RotationalBasis& basis = propagation.rotational_basis();
  const InitialWavePacketInfo& wp = propagation.initial_wave_packet();
  const AnalysisSurface& surface = propagation.analysis_surface();
  const int ichn = basis.initial_index();
  const auto& z = grids.z_specular;

  print_section(out, style, "ROTATIONAL BASIS AND INITIAL WAVE PACKET");
  KeyValueList list;
  list.add("Rotational basis functions", std::to_string(basis.size()));
  list.add("Initial channel (index: J, mJ)",
           std::to_string(ichn) + ": " + std::to_string(basis.j(ichn)) + ", " + std::to_string(basis.mj(ichn)));
  list.add("Parallel kinetic energy", full_precision(propagation.tables().parallel_energy), "hartree");
  list.add("Norm of the initial wave packet", full_precision(wp.norm));
  list.add("|psi(k)| at k = 0", full_precision(wp.k_space_at_zero));
  list.add("|psi(k)| at the grid edge", full_precision(wp.k_space_at_edge));
  list.add("|psi(k)| maximum", full_precision(wp.k_space_maximum));
  list.add("Allowed region (indices)", std::to_string(wp.i_low) + " .. " + std::to_string(wp.i_high));
  list.add("Allowed region (Z)", fixed(z[wp.i_low], 6) + " .. " + fixed(z[wp.i_high], 6), "bohr");
  list.add("Norm inside the allowed region", full_precision(wp.norm_inside));
  if (surface.fixed()) {
    list.add("Analysis surface on the specular grid", "fixed at Z = " + fixed(z[surface.index()], 6), "bohr");
  } else {
    list.add("Analysis surface on the specular grid", "dynamic, starts at Z = " + fixed(z[surface.index()], 6), "bohr");
    list.add("  tolerance for fixing it", full_precision(surface.tolerance()));
  }
  list.print(out, style);
  out << "\n  Indices are 0-based.\n";
}

void AveragesTable::print_header() const {
  print_section(out_, style_, "PROPAGATION");
  out_ << std::setw(12) << "t (a.u.)" << std::setw(25) << "norm" << std::setw(25) << "<Z> (bohr)" << std::setw(25)
       << "<r> (bohr)" << std::setw(25) << "rest norm" << '\n';
}

void AveragesTable::print_row(const Averages& a) const {
  out_ << std::setw(12) << fixed(a.time, 2) << std::setw(25) << full_precision(a.norm) << std::setw(25)
       << full_precision(a.z_mean) << std::setw(25) << full_precision(a.r_mean) << std::setw(25)
       << full_precision(a.rest_norm) << std::endl;
}

void AveragesTable::print_note(const std::string& text) const { out_ << "  " << style_.emphasis(text) << std::endl; }

void write_averages_file(const std::filesystem::path& path, const std::vector<Averages>& averages) {
  std::ofstream file(path);
  if (!file) throw std::runtime_error("cannot open output file '" + path.string() + "'");
  file << "# Wave packet averages (Fortran: averages.out)\n# t(a.u.)  norm  <Z>(bohr)  <r>(bohr)  rest_norm\n";
  for (const Averages& a : averages)
    file << full_precision(a.time) << ' ' << full_precision(a.norm) << ' ' << full_precision(a.z_mean) << ' '
         << full_precision(a.r_mean) << ' ' << full_precision(a.rest_norm) << '\n';
}

}  // namespace spodvr

// Writes and prints the write2DPES cut of the PES (pes_cut_2d.dat).
// Fortran: Potential.f90 Write2DPESToFile (2dpes.dat), Main.f90 (the lines printed after it)
#include "report/pes_cut_file.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

#include "core/unit_conversion.hpp"

namespace spodvr {

void write_pes_cut_file(const std::filesystem::path& path, const PesCut& cut) {
  std::ofstream file(path);
  if (!file) throw std::runtime_error("cannot open output file '" + path.string() + "'");
  file << "# V(r, Z) of the PES at fixed x, y, theta, phi (write2DPES): the PES itself, without the switch to the "
          "gas phase and without the VMin/VMax limits\n";
  file << "# x = " << full_precision(cut.x) << " bohr, y = " << full_precision(cut.y)
       << " bohr (along the lattice vectors)\n";
  file << "# theta = " << full_precision(cut.theta_deg) << " deg (" << full_precision(cut.theta)
       << " rad), phi = " << full_precision(cut.phi_deg) << " deg (" << full_precision(cut.phi) << " rad)\n";
  file << "# r(bohr)  Z(bohr)  V(hartree)  V(eV); Z fastest, a blank line after each r\n";
  const int nr = static_cast<int>(cut.r.size());
  const int nz = static_cast<int>(cut.z.size());
  for (int ir = 0; ir < nr; ++ir) {
    for (int iz = 0; iz < nz; ++iz) {
      const double v = cut.v(ir, iz);
      file << full_precision(cut.r[ir]) << ' ' << full_precision(cut.z[iz]) << ' '
           << full_precision(units::legacy.ev_to_hartree(v)) << ' ' << full_precision(v) << '\n';
    }
    file << '\n';
  }
  if (!file) throw std::runtime_error("cannot write output file '" + path.string() + "'");
}

void print_pes_cut(std::ostream& out, const Style& style, const PesCut& cut, const std::filesystem::path& file) {
  print_section(out, style, "2D CUT OF THE PES (write2DPES)");
  KeyValueList list;
  list.add("x, y (along the lattice vectors)", fixed(cut.x, 6) + ", " + fixed(cut.y, 6), "bohr");
  list.add("theta, phi", fixed(cut.theta_deg, 6) + ", " + fixed(cut.phi_deg, 6), "deg");
  list.add("Points (r x Z)", std::to_string(cut.r.size()) + " x " + std::to_string(cut.z.size()));
  list.add("V", "the PES itself: no switching to the gas phase, no VMin/VMax");
  list.add("Written to", file.string());
  list.print(out, style);
}

}  // namespace spodvr

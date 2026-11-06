// Printed sections and data files of the wave-packet propagation.
#include "report/propagation_report.hpp"

#include <algorithm>
#include <cmath>
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

namespace {

constexpr int time_width = 12;
constexpr int value_width = 25;
constexpr int bar_width = 30;

// text right-aligned in a field of width characters (before any colour codes are added).
std::string right_aligned(const std::string& text, int width) {
  const auto size = static_cast<int>(text.size());
  return size >= width ? text : std::string(width - size, ' ') + text;
}

std::string repeated(std::string_view piece, int count) {
  std::string text;
  for (int i = 0; i < count; ++i) text += piece;
  return text;
}

}  // namespace

AveragesTable::AveragesTable(std::ostream& out, const Style& style, int n_steps, double time_step, double t_start)
    : out_(out),
      style_(style),
      n_steps_(n_steps),
      time_step_(time_step),
      t_start_(t_start),
      start_(std::chrono::steady_clock::now()) {}

void AveragesTable::print_header() const {
  print_section(out_, style_, "PROPAGATION");
  const std::string header = right_aligned("t (a.u.)", time_width) + right_aligned("norm", value_width) +
                             right_aligned("<Z> (bohr)", value_width) + right_aligned("<r> (bohr)", value_width) +
                             right_aligned("rest norm", value_width);
  out_ << style_.emphasis(header) << '\n';
  if (style_.color) out_ << style_.dim(repeated("\u2500", time_width + 4 * value_width)) << '\n';
}

void AveragesTable::print_row(const Averages& a) {
  clear_status();
  out_ << style_.accent(right_aligned(fixed(a.time, 2), time_width)) << std::setw(value_width) << full_precision(a.norm)
       << std::setw(value_width) << full_precision(a.z_mean) << std::setw(value_width) << full_precision(a.r_mean)
       << std::setw(value_width) << full_precision(a.rest_norm) << std::endl;
}

void AveragesTable::print_note(const std::string& text, double time) {
  clear_status();
  out_ << "  " << style_.warning("\u25b8 " + text) << style_.dim(" (t = " + fixed(time, 2) + " a.u.)") << std::endl;
}

void AveragesTable::print_status(double time) {
  if (!style_.color) return;
  const long long step = std::llround((time - t_start_) / time_step_);
  const double fraction = n_steps_ > 0 ? std::clamp(static_cast<double>(step) / n_steps_, 0.0, 1.0) : 1.0;
  const double seconds_per_step = step > 0 ? elapsed_seconds() / static_cast<double>(step) : 0.0;
  const double remaining = seconds_per_step * static_cast<double>(n_steps_ - step);
  const auto filled = static_cast<int>(fraction * bar_width);

  out_ << "\r\033[K  \u2595" << style_.accent(repeated("\u2588", filled))
       << style_.dim(repeated("\u2591", bar_width - filled)) << "\u258f "
       << right_aligned(fixed(100.0 * fraction, 0), 3) << "%  step " << step << '/' << n_steps_
       << "  t = " << fixed(time, 1) << " a.u.  " << fixed(seconds_per_step, 3) << " s/step  ETA "
       << format_duration(remaining) << std::flush;
  status_shown_ = true;
}

void AveragesTable::finish() {
  clear_status();
  const double elapsed = elapsed_seconds();
  out_ << "\n  Propagation time: " << fixed(elapsed, 1) << " s (" << format_duration(elapsed) << ", " << n_steps_
       << " steps";
  if (n_steps_ > 0) out_ << ", " << fixed(elapsed / n_steps_, 3) << " s/step";
  out_ << ")\n";
}

void AveragesTable::clear_status() {
  if (!status_shown_) return;
  out_ << "\r\033[K";
  status_shown_ = false;
}

double AveragesTable::elapsed_seconds() const {
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_;
  return elapsed.count();
}

void write_averages_file(const std::filesystem::path& path, const std::vector<Averages>& averages) {
  std::ofstream file(path);
  if (!file) throw std::runtime_error("cannot open output file '" + path.string() + "'");
  file << "# Wave packet averages (Fortran: averages.out)\n# t(a.u.)  norm  <Z>(bohr)  <r>(bohr)  rest_norm\n";
  for (const Averages& a : averages)
    file << full_precision(a.time) << ' ' << full_precision(a.norm) << ' ' << full_precision(a.z_mean) << ' '
         << full_precision(a.r_mean) << ' ' << full_precision(a.rest_norm) << '\n';
}

}  // namespace spodvr

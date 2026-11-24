// Printed sections of the asymptotic analysis: set-up and compact result tables.
#include "report/analysis_report.hpp"

#include <cmath>
#include <complex>
#include <string>
#include <vector>

namespace spodvr {

namespace {

constexpr int max_printed_energies = 50;  // longer energy lists are thinned on screen
constexpr int max_rovib_columns = 12;

std::string state_label(int v, int j, int mj, int n, int m) {
  return "v=" + std::to_string(v) + " j=" + std::to_string(j) + " mJ=" + std::to_string(mj) + " (n,m)=(" +
         std::to_string(n) + "," + std::to_string(m) + ")";
}

// Every k-th energy, so that at most max_printed_energies rows are printed.
int energy_stride(int n_energies) { return (n_energies + max_printed_energies - 1) / max_printed_energies; }

// Rovibrational state (v, j), summed over mJ and diffraction.
struct RovibColumn {
  int v;
  int j;
};

// The (v, j) states that are open at some analysis energy, in order of v, then j.
std::vector<RovibColumn> open_rovib_states(const AnalysisBasis& basis, const ScatteringResults& results) {
  std::vector<RovibColumn> columns;
  for (int v = 0; v < basis.n_v(); ++v)
    for (int j = 0; j <= basis.j_max(); ++j) {
      bool open = false;
      for (const EnergyResult& e : results.energies)
        for (int is = 0; is < basis.n_states(); ++is)
          if (basis.j(is) == j && e.open[is + basis.n_states() * v]) open = true;
      if (open) columns.push_back({v, j});
    }
  return columns;
}

// P(v, j) summed over mJ (and diffraction), and whether any of its channels is open.
double rovib_probability(const AnalysisBasis& basis, const EnergyResult& e, RovibColumn c, bool& open) {
  double p = 0.0;
  open = false;
  for (int irotana = 0; irotana < basis.n_rotation(); ++irotana) {
    const int is = irotana * basis.n_diffraction();
    if (basis.j(is) != c.j) continue;
    p = p + e.p_rotation[irotana + basis.n_rotation() * c.v];
    for (int idif = 0; idif < basis.n_diffraction(); ++idif)
      if (e.open[is + idif + basis.n_states() * c.v]) open = true;
  }
  return p;
}

void print_summary_table(std::ostream& out, const Style& style, const AnalysisBasis& basis,
                         const ScatteringResults& results) {
  std::vector<std::string> headers = {"E_coll (eV)", "P_scattered", "P_reaction"};
  for (int v = 0; v < basis.n_v(); ++v) headers.push_back("P(v=" + std::to_string(v) + ")");
  for (const std::string h : {"<E_vib> (eV)", "<E_rot> (eV)", "<E_dif> (eV)", ""}) headers.push_back(h);
  Table table(headers);
  const int stride = energy_stride(static_cast<int>(results.energies.size()));
  for (std::size_t i = 0; i < results.energies.size(); i += stride) {
    const EnergyResult& e = results.energies[i];
    std::vector<std::string> row = {fixed(e.collision_energy, 4), fixed(e.p_total, 7), fixed(e.p_reaction, 7)};
    for (const double p : e.p_vibration) row.push_back(fixed(p, 7));
    for (const double x : {e.e_vibration, e.e_rotation, e.e_diffraction}) row.push_back(fixed(x, 6));
    row.push_back(e.p_total > 1.01 ? "sum > 1.01" : "");
    table.add_row(row);
  }
  table.print(out, style);
}

void print_rovib_table(std::ostream& out, const Style& style, const AnalysisBasis& basis,
                       const ScatteringResults& results) {
  std::vector<RovibColumn> columns = open_rovib_states(basis, results);
  const bool truncated = static_cast<int>(columns.size()) > max_rovib_columns;
  if (truncated) columns.resize(max_rovib_columns);
  std::vector<std::string> headers = {"E_coll (eV)"};
  for (const RovibColumn& c : columns) headers.push_back("v=" + std::to_string(c.v) + " j=" + std::to_string(c.j));
  Table table(headers);
  const int stride = energy_stride(static_cast<int>(results.energies.size()));
  for (std::size_t i = 0; i < results.energies.size(); i += stride) {
    const EnergyResult& e = results.energies[i];
    std::vector<std::string> row = {fixed(e.collision_energy, 4)};
    for (const RovibColumn& c : columns) {
      bool open = false;
      const double p = rovib_probability(basis, e, c, open);
      row.push_back(open ? fixed(p, 7) : "-");
    }
    table.add_row(row);
  }
  table.print(out, style);
  out << "  Summed over mJ and diffraction; '-': closed.";
  if (truncated) out << " Further states: probabilities_rovib.tsv.";
  out << '\n';
}

void print_diffraction_table(std::ostream& out, const Style& style, const Config& config, const AnalysisBasis& basis,
                             const ScatteringResults& results) {
  const int v = config.input.v_in;
  const int j = config.input.j_in;
  const int n_orders = basis.max_order() + 1;
  std::vector<std::string> headers = {"E_coll (eV)"};
  for (int order = 0; order < n_orders; ++order) headers.push_back("order " + std::to_string(order));
  Table table(headers);
  const int stride = energy_stride(static_cast<int>(results.energies.size()));
  for (std::size_t i = 0; i < results.energies.size(); i += stride) {
    const EnergyResult& e = results.energies[i];
    std::vector<std::string> row = {fixed(e.collision_energy, 4)};
    for (int order = 0; order < n_orders; ++order) {
      double p = 0.0;
      for (int irotana = 0; irotana < basis.n_rotation(); ++irotana)
        if (basis.j(irotana * basis.n_diffraction()) == j)
          p = p + e.p_order[order + n_orders * (irotana + basis.n_rotation() * v)];
      row.push_back(fixed(p, 7));
    }
    table.add_row(row);
  }
  table.print(out, style);
}

// A number in full precision; a NaN as "nan" whatever its sign bit.
std::string number(double x) { return std::isnan(x) ? "nan" : full_precision(x); }

// "(v, j, mJ, n, m) = (0, 2, -1, 1, 0)" for analysed state (is, v).
std::string state_labels(const AnalysisBasis& basis, int is, int v) {
  const DiffractionChannel& c = basis.channel(is);
  return "(v, j, mJ, n, m) = (" + std::to_string(v) + ", " + std::to_string(basis.j(is)) + ", " +
         std::to_string(basis.mj(is)) + ", " + std::to_string(c.n) + ", " + std::to_string(c.m) + ")";
}

// "(v, j, mJ) = (0, 2, -1)" for analysed rotational state irotana in level v.
std::string rotation_labels(const AnalysisBasis& basis, int irotana, int v) {
  const int is = irotana * basis.n_diffraction();  // any state of that rotation: its first diffraction channel
  return "(v, j, mJ) = (" + std::to_string(v) + ", " + std::to_string(basis.j(is)) + ", " +
         std::to_string(basis.mj(is)) + ")";
}

// The quantity at value.index of one energy's results, with its labels, and the number itself, as
// "<quantity> is <number>". Layouts: EnergyResult.
std::string non_finite_quantity(const NonFiniteValue& value, const AnalysisBasis& basis, const EnergyResult& e) {
  const int i = value.index;
  const int n_states = basis.n_states();
  const int n_rot = basis.n_rotation();
  const int n_orders = basis.max_order() + 1;
  switch (value.quantity) {
    case ResultQuantity::total_energy: return "E_total is " + number(e.total_energy);
    case ResultQuantity::collision_energy: return "E_coll is " + number(e.collision_energy);
    case ResultQuantity::z_energy: return "E_z is " + number(e.z_energy);
    case ResultQuantity::total_momentum: return "the total momentum is " + number(e.total_momentum);
    case ResultQuantity::s_matrix: {
      const std::complex<double> s = e.s_matrix[i];
      return "S of " + state_labels(basis, i % n_states, i / n_states) + " is (" + number(s.real()) + ", " +
             number(s.imag()) + ")";
    }
    case ResultQuantity::probability:
      return "P of " + state_labels(basis, i % n_states, i / n_states) + " is " + number(e.probability[i]);
    case ResultQuantity::phase:
      return "the phase of " + state_labels(basis, i % n_states, i / n_states) + " is " + number(e.phase[i]);
    case ResultQuantity::p_vibration: return "P of v = " + std::to_string(i) + " is " + number(e.p_vibration[i]);
    case ResultQuantity::p_rotation:
      return "P of " + rotation_labels(basis, i % n_rot, i / n_rot) + " is " + number(e.p_rotation[i]);
    case ResultQuantity::p_order: {
      const int order = i % n_orders;
      const int rotation = i / n_orders;
      return "P of " + rotation_labels(basis, rotation % n_rot, rotation / n_rot) + ", diffraction order " +
             std::to_string(order) + " is " + number(e.p_order[i]);
    }
    case ResultQuantity::e_rotation_v:
      return "<E_rot> of v = " + std::to_string(i) + " is " + number(e.e_rotation_v[i]);
    case ResultQuantity::e_diffraction_v:
      return "<E_dif> of v = " + std::to_string(i) + " is " + number(e.e_diffraction_v[i]);
    case ResultQuantity::p_total: return "P_scattered is " + number(e.p_total);
    case ResultQuantity::p_reaction: return "P_reaction is " + number(e.p_reaction);
    case ResultQuantity::e_vibration: return "<E_vib> is " + number(e.e_vibration);
    case ResultQuantity::e_rotation: return "<E_rot> is " + number(e.e_rotation);
    case ResultQuantity::e_diffraction: return "<E_dif> is " + number(e.e_diffraction);
  }
  return "an unknown quantity is not finite";
}

}  // namespace

void print_analysis_setup(std::ostream& out, const Style& style, const Config& config, const AnalysisBasis& basis,
                          const AnalysisEnergies& energies) {
  const InputParameters& p = config.input;
  print_section(out, style, "ASYMPTOTIC ANALYSIS");
  KeyValueList list;
  list.add("Vibrational levels analysed (nviba)", std::to_string(basis.n_v()));
  list.add("Highest j analysed (jana)", std::to_string(basis.j_max()));
  list.add("Rotational states analysed (nana)", std::to_string(basis.n_rotation()));
  list.add("Highest diffraction order analysed (mdifa)", std::to_string(basis.max_order()));
  list.add("Diffraction channels analysed (ndifa)",
           std::to_string(basis.n_diffraction()) + " of " + std::to_string(basis.diffraction_channels().size()));
  list.add("Rotation-diffraction states (nstata)", std::to_string(basis.n_states()));
  const int is0 = basis.initial_state();
  list.add("Initial channel (index)",
           std::to_string(is0) + ": " + state_label(p.v_in, basis.j(is0), basis.mj(is0), 0, 0));
  list.add("Analysis surface (scattering grid)", fixed(p.z_analysis, 6), "bohr");
  if (energies.size() > 0) {
    list.add("Analysis energies (collision)",
             fixed(energies.collision.front(), 6) + " .. " + fixed(energies.collision.back(), 6), "eV");
    const double step = energies.size() > 1 ? energies.collision[1] - energies.collision[0] : 0.0;
    list.add("  number, step", std::to_string(energies.size()) + ", " + fixed(step, 6), "eV");
  }
  list.print(out, style);
  if (basis.diffraction_basis_too_small())
    out << "\n  " << style.warning("Warning")
        << ": the diffraction basis is too small to describe all open channels.\n";
}

std::string energy_thinning_note(int n_energies) {
  const int stride = energy_stride(n_energies);
  if (stride <= 1) return "";
  return "One energy in " + std::to_string(stride) + " is shown; all are in the data files.";
}

void print_scattering_results(std::ostream& out, const Style& style, const Config& config, const AnalysisBasis& basis,
                              const ScatteringResults& results) {
  print_section(out, style, "SCATTERING PROBABILITIES");
  const std::string thinning_note = energy_thinning_note(static_cast<int>(results.energies.size()));
  if (!thinning_note.empty()) out << "  " << thinning_note << "\n\n";
  out << style.emphasis("  Summary per energy") << "\n\n";
  print_summary_table(out, style, basis, results);
  out << '\n' << style.emphasis("  Rovibrational state-resolved probabilities") << "\n\n";
  print_rovib_table(out, style, basis, results);
  out << '\n'
      << style.emphasis("  Diffraction orders of the initial state (v=" + std::to_string(config.input.v_in) +
                        ", j=" + std::to_string(config.input.j_in) + ", summed over mJ)")
      << "\n\n";
  print_diffraction_table(out, style, config, basis, results);
}

std::string describe_non_finite(std::string_view set, const NonFiniteValue& value, const AnalysisBasis& basis,
                                const ScatteringResults& results) {
  const EnergyResult& e = results.energies[value.energy_index];
  return std::string(set) + ", energy " + std::to_string(value.energy_index + 1) + " of " +
         std::to_string(results.energies.size()) + " (E_coll = " + fixed(e.collision_energy, 6) +
         " eV): " + non_finite_quantity(value, basis, e);
}

}  // namespace spodvr

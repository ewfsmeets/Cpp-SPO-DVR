// Printed sections of the flux analysis: its set-up, its memory rule and the reaction probability.
// Fortran: FluxAnalyzor.f90 setupFluxAnalyzor (the set-up printout); the results replace no Fortran
// printout (the Fortran writes them to files only).
#include "report/flux_report.hpp"

#include <string>

#include "core/unit_conversion.hpp"

namespace spodvr {

namespace {

// The flux analysis's energies are converted as in FluxAnalyzor.f90 (flux Q5).
constexpr const UnitConversion& flux_units = units::precise;

// "N (with 1 slot)", or why there is no such number.
std::string max_energies_text(const FluxMemoryCheck& memory) {
  if (memory.max_energies) return std::to_string(*memory.max_energies) + " (with 1 slot)";
  if (!memory.physical_bytes) return "unknown (physical memory unknown)";
  return "not limited here: the run exceeds the physical memory without the flux analysis too";
}

}  // namespace

void print_flux_analysis_setup(std::ostream& out, const Style& style, const Config& config, const SpatialGrids& grids,
                               const FluxEnergies& energies, double initial_state_ratio,
                               const FluxMemoryCheck& memory) {
  const int ir_flux = *grids.ir_flux_analysis;
  print_section(out, style, "FLUX ANALYSIS");
  KeyValueList list;
  list.add("Flux surface r_c (index)", fixed(grids.r[ir_flux], 6) + " (" + std::to_string(ir_flux) + ")", "bohr");
  list.add("Z points summed (Z < Z asymptotic)", std::to_string(grids.nz_potential));
  if (energies.size() > 0) {
    list.add("Flux energies (E_z)",
             fixed(flux_units.hartree_to_ev(energies.z_energy.front()), 6) + " .. " +
                 fixed(flux_units.hartree_to_ev(energies.z_energy.back()), 6),
             "eV");
    list.add("  number, step", std::to_string(energies.size()) + ", " + fixed(config.input.delta_e_flux_analysis, 6),
             "eV");
  }
  list.add("chi_vin(r_c) / max |chi_vin|", scientific(initial_state_ratio, 3));
  list.add("Slots (analysis steps per flush)", std::to_string(memory.plan.slot_count));
  list.add("Memory per flux energy", format_bytes(memory.items.bytes_per_energy));
  list.add("Most flux energies that fit", max_energies_text(memory));
  list.print(out, style);
  out << "\n  Indices are 0-based. P_flux: the reaction probability from the flux through r = r_c.\n";
  for (const std::string& w : energies.warnings) out << "  " << style.warning("Warning") << ": " << w << '\n';
  if (!(initial_state_ratio <= flux_initial_state_ratio_limit))
    out << "  " << style.warning("Warning")
        << ": the flux analysis leaves out the initial channel (J_in, mJ_in, specular) at r_c, as the "
           "Fortran does (flux Q7); chi_vin(r_c) / max |chi_vin| = "
        << scientific(initial_state_ratio, 3) << " bounds only the incoming packet's part of it.\n";
}

void print_flux_results(std::ostream& out, const Style& style, const FluxResults& results) {
  print_section(out, style, "FLUX ANALYSIS: REACTION PROBABILITY");
  const int n_energies = static_cast<int>(results.energy.size());
  const int stride = energy_stride(n_energies);
  if (stride > 1) out << "  Every " << stride << "th energy is shown; all are in flux_reaction.tsv.\n\n";
  Table table({"E_z (eV)", "P_flux"});
  for (int n = 0; n < n_energies; n += stride)
    table.add_row({fixed(results.energy[n], 4), fixed(results.probability[n], 7)});
  table.print(out, style);
}

void print_flux_memory(std::ostream& out, const Style& style, const FluxMemoryCheck& memory) {
  const FluxMemoryItems& f = memory.items;
  out << '\n' << style.emphasis("  Flux analysis (host memory)") << "\n\n";
  KeyValueList list;
  const std::string base_label = memory.device == Device::Cuda ? "Run without it (GPU host): set-up / time loop"
                                                               : "Run without it: set-up / time loop";
  list.add(base_label, format_bytes(memory.base.setup) + " / " + format_bytes(memory.base.time_loop));
  list.add("Coefficients", std::to_string(f.n_energies) + " energies x " + format_bytes(f.bytes_per_energy) + " = " +
                               format_bytes(f.coefficient_bytes()));
  list.add("Slots", std::to_string(memory.plan.slot_count) + " x " + format_bytes(f.slot_bytes));
  list.add("Physical memory", memory.physical_bytes ? format_bytes(*memory.physical_bytes) : "unknown");
  list.add("Most flux energies that fit", max_energies_text(memory));
  list.print(out, style);
  if (!memory.plan.fits)
    out << "\n  " << style.warning("The flux analysis does not fit in the physical memory") << ": the run stops.\n";
}

}  // namespace spodvr

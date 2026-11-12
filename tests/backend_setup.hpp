#pragma once
// Shared test set-up: everything a propagation backend is built from, for a reference case of
// tests/reference/ with changed input keys (e.g. the flux analysis keys). Used by the CPU tests
// and by the GPU tests (tests/cuda/, which build the GPU backend from it).

#include <fstream>
#include <memory>
#include <string>

#include <nlohmann/json.hpp>

#include "fft/fftw_fft.hpp"
#include "gasphase/gas_phase.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "grids/spatial_grids.hpp"
#include "input/config.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "linalg/symmetric_eigensolver.hpp"
#include "potential/potential_factory.hpp"
#include "potential/potential_grid.hpp"
#include "potential/potential_symmetry.hpp"
#include "propagation/cpu_backend.hpp"
#include "propagation/initial_wave_packet.hpp"
#include "propagation/propagation_backend.hpp"
#include "propagation/propagation_tables.hpp"
#include "propagation/wave_function.hpp"

namespace spodvr::testing {

// The input.json of the reference case tests/reference/<case_name>.
inline nlohmann::json reference_case_input(const std::string& case_name) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + case_name + "/input.json");
  return nlohmann::json::parse(file);
}

// The input keys that switch on the flux analysis of a reference case: the flux surface at
// rFluxAnalysis = r_flux (bohr, a point of the case's r grid, <= ROptStart) and the energies
// EMinWP .. EMaxWP of its wave packet in steps of delta_e (eV).
inline nlohmann::json flux_analysis_keys(const std::string& case_name, double r_flux, double delta_e) {
  const nlohmann::json input = reference_case_input(case_name);
  return {{"doFluxAnalysis", true},
          {"rFluxAnalysis", r_flux},
          {"EMinFluxAna", input.at("EMinWP")},
          {"EMaxFluxAna", input.at("EMaxWP")},
          {"DeltaEFluxAna", delta_e}};
}

// Everything a backend is built from, for a reference case with `changes` applied to its input
// keys: grids, potential (with the case's surface symmetry), gas phase, tables and the initial
// wave packet on the specular grid.
struct BackendSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  FftwFft fft;
  RotationalBasis basis;
  PotentialGrid potential;
  GasPhase gas;
  PropagationTables tables;
  SpecularWave initial;

  static Config read(const std::string& case_name, const nlohmann::json& changes) {
    nlohmann::json json = reference_case_input(case_name);
    json.update(changes);
    return build_config(parse_input(json));
  }

  explicit BackendSetup(const std::string& case_name, const nlohmann::json& changes = nlohmann::json::object())
      : config(read(case_name, changes)),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)),
        basis(angular, config.input.j_in, config.input.mj_in),
        potential(compute_potential_grid(*make_potential(config.input), config, grids, angular,
                                         make_grid_symmetry(config, angular))),
        gas(compute_gas_phase(config, grids, angular, potential.asymptotic, fft, LapackEigensolver())),
        tables(make_propagation_tables(config, grids, basis, potential.asymptotic)),
        initial(make_initial_wave_packet(config, grids, gas, fft).specular) {}

  // The references in the result point into this set-up, which must outlive every backend built from it.
  [[nodiscard]] BackendInput input() const { return {config, grids, angular, basis, tables, fft}; }
  [[nodiscard]] std::unique_ptr<PropagationBackend> cpu() const {
    return make_cpu_backend(input(), potential, initial);
  }
};

}  // namespace spodvr::testing

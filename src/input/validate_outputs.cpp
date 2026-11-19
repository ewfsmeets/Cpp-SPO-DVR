// Input checks of the plot-data keys (C++ only): writeDensities. Each check names its key (exit 2).
#include "input/validate_outputs.hpp"

#include "input/validation_helpers.hpp"

namespace spodvr {

using validation::check;

void validate_output_keys(const Config& config) {
  const InputParameters& p = config.input;
  const bool propagates = config.derived.run_mode == 0 || config.derived.run_mode == 1;
  check(p.write_densities && !propagates, "writeDensities",
        "only a propagating run (a new run or a restart) writes densities");
}

}  // namespace spodvr

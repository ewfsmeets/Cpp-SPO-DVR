#pragma once

#include <string>

namespace spodvr {

// Configures the external numerical libraries once, before their first use (called by the FFT
// and eigensolver backends; repeated calls do nothing).
// With oneMKL: conditional numerical reproducibility (CNR). Unless the environment variable
// MKL_CBWR is set (it is then honoured by oneMKL itself), the branch AUTO is selected: the best
// code path for this CPU, with results that are reproducible from run to run on the same
// machine. They may differ in the last bits between machines with different vector extensions;
// set MKL_CBWR=AVX2 (or COMPATIBLE) for identical results across machines. See
// docs/performance.md.
void configure_numerical_libraries();

// One line naming the BLAS/LAPACK and FFT libraries in use and the oneMKL CNR branch.
[[nodiscard]] std::string numerical_libraries_summary();

}  // namespace spodvr

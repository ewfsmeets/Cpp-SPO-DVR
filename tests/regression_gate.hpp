#pragma once

// Gates 3 and 4 (docs/validation.md) for one reference case, with a given propagation backend: the C++ run
// of <case>/input.json against the Fortran output of the same case (skipped when missing).
// Used by test_regression_propagation.cpp (CPU) and cuda/test_cuda_regression.cpp (GPU).

#include <filesystem>
#include <iostream>
#include <string>

#include "analysis_reference.hpp"
#include "propagation_reference.hpp"

namespace spodvr::testing {

// Gate 4 tolerance of the BK coefficients, as a fraction of the reference scale (compare_bk_with_fortran).
inline constexpr double bk_gate_tolerance = 1e-10;

inline void run_gate(const std::string& case_name, const BackendFactory& backend = make_cpu_backend) {
  namespace fs = std::filesystem;
  const fs::path dir = fs::path(SPODVR_REFERENCE_DIR) / case_name;
  const auto ref = read_fortran_reference(dir / "fortran_stdout.txt", dir / "fortran_averages.out");
  if (!ref)
    GTEST_SKIP() << "no Fortran reference in " << dir << "; run "
                 << "tools/fortran_reference/run_propagation_reference.sh " << case_name;
  ASSERT_FALSE(ref->averages.empty()) << "no averages found in fortran_averages.out";
  const PropagationRun cpp = run_cpp_propagation(dir / "input.json", backend);
  compare_with_fortran(cpp, *ref);  // gate 3

  // Gate 4: the printed probabilities, and the BK coefficients of every analysis step.
  compare_probabilities_with_fortran(cpp, read_fortran_probabilities(dir / "fortran_stdout.txt"));
  if (fs::exists(dir / "fortran_fort21.bin"))
    compare_bk_with_fortran(cpp, read_fortran_bk_records(dir / "fortran_fort21.bin"), bk_gate_tolerance);
  else
    std::cout << "  no fortran_fort21.bin: BK coefficients not compared\n";

  // Restart gate: a Fortran run split by a restart gives the probabilities of the unsplit run
  // (tools/fortran_reference/run_restart_reference.sh). The C++ restart is bitwise exact (test_restart).
  if (fs::exists(dir / "fortran_restart_stdout.txt"))
    compare_probabilities_with_fortran(cpp, read_fortran_probabilities(dir / "fortran_restart_stdout.txt"));
}

}  // namespace spodvr::testing

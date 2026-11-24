#pragma once

// Gates 3 and 4 (docs/validation.md) for one reference case, with a given propagation backend: the C++ run
// of <case>/input.json against the Fortran output of the same case (skipped or failed when missing).
// Used by test_regression_propagation.cpp (CPU) and cuda/test_cuda_regression.cpp (GPU).

#include <filesystem>
#include <iostream>
#include <map>
#include <string>

#include "analysis_reference.hpp"
#include "fort25_reference.hpp"
#include "propagation_reference.hpp"

namespace spodvr::testing {

// Gate 4 tolerance of the BK coefficients, as a fraction of the reference scale (compare_bk_with_fortran).
inline constexpr double bk_gate_tolerance = 1e-10;

// Whether a case's Fortran files may be missing: without them the gates of an optional case skip
// (without fort21 they leave out the BK comparison), those of a required case fail.
enum class References { optional, required };

// Gate 4, the S-matrix elements: the code-difference term of the tolerance (Fort25Tolerance::s_scale_fraction), the
// smallest power of 10 at least 10 times the largest difference measured over every case with the CPU, the GPU and
// their fast-math builds (2026-11-24, as a fraction of max(|S|, s_floor = 1e-3)): 5.99e-12 in hardwall_short, whose
// elements are all rounding noise (|S| <= 2.4e-15; 6.0e-15 absolute), with the fast-math CPU build; 1.1e-12 in the
// Dai-Light cases with fast-math; at most 2.1e-13 with the default builds.
inline constexpr double fort25_code_difference = 1e-10;

// The counts gate 4 must reach per case with that tolerance (measured; the same with every build): the rows whose
// phase is compared, and the rows that fail with the initial state's phase referred to Zana, the convention of the
// Fortran up to r1719 (a realistic porting mistake the gate must catch): exactly the rows r1720 changed, 0 in
// hardwall_short, whose specular line never fixes.
inline const std::map<std::string, Fort25Pins>& fort25_pins() {
  static const std::map<std::string, Fort25Pins> pins = {
      {"hardwall_short", {0, 0}},
      {"hardwall_small", {51, 51}},
      {"cu111_dailight_short", {31624, 50}},
      {"cu111_dailight_short_sym", {31584, 50}},
      {"cu111_dailight_short_lowz", {31621, 50}},
      {"cu111_dailight_small", {36534, 50}},
      {"scm_cu111_srp48_small", {9486, 21}},
  };
  return pins;
}

// Gate 4 (S-matrix elements): every printed row of the Fortran's on-the-fly fort.25 against the C++ S, P and phase
// (compare_fort25), with the counts and pins of expect_fort25_gate.
inline void compare_fort25_with_fortran(const PropagationRun& cpp, const Fort25File& file,
                                        const std::string& case_name) {
  Fort25Tolerance tolerance;
  tolerance.s_scale_fraction = fort25_code_difference;
  const Fort25Comparison c = compare_fort25(cpp, file, tolerance);
  std::cout << "  fort.25: " << c.rows << " rows (print rule " << c.expected_rows << ", file " << c.row_like_lines
            << "), " << c.compared << " compared, " << c.phase_compared << " phases; largest code difference "
            << c.worst_ratio << " of max(|S|, " << tolerance.s_floor
            << ") = " << nan_aware_max(c.s_scale, tolerance.s_floor) << "\n";
  const auto pin = fort25_pins().find(case_name);
  ASSERT_NE(pin, fort25_pins().end()) << "no pinned gate-4 counts for case " << case_name;
  PropagationRun zana = cpp;
  zana.scattering = cpp.scattering_zana;
  expect_fort25_gate(c, compare_fort25(zana, file, tolerance).failed, pin->second);
}

inline void run_gate(const std::string& case_name, const BackendFactory& backend = make_cpu_backend,
                     References references = References::optional) {
  namespace fs = std::filesystem;
  const fs::path dir = fs::path(SPODVR_REFERENCE_DIR) / case_name;
  const fs::path stdout_file = dir / "fortran_stdout.txt";
  const fs::path averages_file = dir / "fortran_averages.out";
  const fs::path bk_file = dir / "fortran_fort21.bin";
  const fs::path restart_file = dir / "fortran_restart_stdout.txt";  // optional in every case
  const fs::path fort25_file = dir / "fortran_fort25.txt";
  if (references == References::required) {
    bool missing = false;  // one failure per missing file, then stop
    for (const fs::path& file : {stdout_file, averages_file, bk_file, fort25_file})
      if (!fs::exists(file)) {
        ADD_FAILURE() << "reference file " << file << " missing";
        missing = true;
      }
    if (missing) return;
  }
  const auto ref = read_fortran_reference(stdout_file, averages_file);
  if (!ref)
    GTEST_SKIP() << "no Fortran reference in " << dir << "; run "
                 << "tools/fortran_reference/run_propagation_reference.sh " << case_name;
  ASSERT_FALSE(ref->averages.empty()) << "no averages found in fortran_averages.out";
  const PropagationRun cpp = run_cpp_propagation(dir / "input.json", backend);
  compare_with_fortran(cpp, *ref);  // gate 3

  // Gate 4: the printed probabilities, and the BK coefficients of every analysis step.
  compare_probabilities_with_fortran(cpp, read_fortran_probabilities(stdout_file));
  if (fs::exists(bk_file))
    compare_bk_with_fortran(cpp, read_fortran_bk_records(bk_file), bk_gate_tolerance);
  else
    std::cout << "  no fortran_fort21.bin: BK coefficients not compared\n";
  // Gate 4: the S-matrix elements of every printed state at every energy (a case with a stdout has its fort.25).
  ASSERT_TRUE(fs::exists(fort25_file)) << "reference file " << fort25_file << " missing";
  compare_fort25_with_fortran(cpp, read_fortran_fort25_rows(fort25_file), case_name);

  // Restart gate: a Fortran run split by a restart gives the probabilities of the unsplit run
  // (tools/fortran_reference/run_restart_reference.sh). The C++ restart is bitwise exact (test_restart).
  if (fs::exists(restart_file)) compare_probabilities_with_fortran(cpp, read_fortran_probabilities(restart_file));
}

}  // namespace spodvr::testing

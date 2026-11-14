// Tests: the SCM correction ΔV (src/potential/scm/scm_correction.*) on data built in memory: exact zeros, the
// switch and r clamp of f, one displaced atom, translation and exchange symmetry, each tolerance test with a
// negative control (a porting mistake that fails it). Against the Fortran: tests/scm/test_scm_correction_reference.cpp
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "comparison.hpp"
#include "potential/scm/scm_correction.hpp"
#include "potential/scm/scm_files.hpp"

using namespace spodvr;
using spodvr::testing::nan_aware_max;

namespace {

constexpr double pi = std::numbers::pi;
// Hartree, for the analytic and symmetry tests: the rounding of the atom sums measured here is ≤ 6e-16, the
// smallest deviation of a negative control 1.1e-10 (the layer mix-up of the analytic test).
constexpr double tolerance = 1e-14;

// The 18 values of a layer (P1–P8, P9–P16, P17, P18), as in params.dat.
ScmLayerParameters layer_parameters(const std::vector<double>& values) {
  ScmLayerParameters layer;
  for (int k = 0; k < 8; ++k) {
    layer.intercept[k] = values[k];
    layer.slope[k] = values[8 + k];
  }
  layer.r_clamp_min = values[16];
  layer.r_clamp_max = values[17];
  return layer;
}

// Layer 1: the example params.dat of the SRP48 data set (both of its layers are this one).
// Layer 2: that of tests/scm/reference/synthetic/params.in, unlike layer 1, so that a mix-up of
// the layers shows.
ScmParameters parameters() {
  ScmParameters result;
  result.version = 2;
  result.layers[0] = layer_parameters({2.28970941935, 1.29096748614, -0.0703557026734, 0.0234889163443,
                                       -0.0632740237351, 0.0272448940344, 7.44024869632, 10.4600719899, 0.0236010392861,
                                       -0.00311521325544, 0.0165908842195, 0.0286500845806, -0.00722325023675,
                                       -0.00644192925251, 0.000769448712731, 4.63492299557e-05, 0.725, 1.847});
  result.layers[1] = layer_parameters(
      {2.63970941935, 1.097322363219, -0.04221342160403999, 0.014093349806579999, -0.03796441424105999,
       0.016346936420639998, 6.840248696320001, 9.9600719899, 0.01180051964305, -0.00155760662772, 0.00829544210975,
       0.0143250422903, -0.003611625118375, -0.003220964626255, 0.0003847243563655, 2.317461497785e-05, 0.8, 1.7});
  return result;
}

// The same parameters with layers 1 and 2 swapped: a correction that indexes the wrong layer.
ScmParameters swapped_layers() {
  ScmParameters result = parameters();
  std::swap(result.layers[0], result.layers[1]);
  return result;
}

constexpr double lattice_constant = 4.9;   // bohr, the nearest-neighbour distance of the slab
constexpr double skewing_angle = pi / 3.0;  // the hexagonal cell, as the C3v CRP data set gives it
constexpr double rescaling = 1.0125;       // an a_T/a₀ for the rescaled runs

constexpr int half_width = 3;                     // the slab's sites i, j run from −3 to 3
constexpr int sites_per_row = 2 * half_width + 1;  // 7
constexpr int atoms_per_layer = sites_per_row * sites_per_row;

// The index in ideal_slab() of site (i, j) of a layer.
int atom_index(int layer, int i, int j) {
  return (layer - 1) * atoms_per_layer + (j + half_width) * sites_per_row + (i + half_width);
}

// An undistorted FCC(111) slab: layer 1 at z = 0 on the sites i a₁ + j a₂ (|i|, |j| ≤ 3), layer 2
// a√(2/3) below it under the hcp sites. Distorted and ideal positions are equal.
ScmSlab ideal_slab() {
  ScmSlab slab;
  for (int layer = 1; layer <= 2; ++layer) {
    const double shift = (layer - 1) / 3.0;
    for (int j = -half_width; j <= half_width; ++j) {
      for (int i = -half_width; i <= half_width; ++i) {
        const double u = (i + shift) * lattice_constant;
        const double v = (j + shift) * lattice_constant;
        const Vec3 site{u + 0.5 * v, std::sqrt(3.0) / 2.0 * v, -(layer - 1) * std::sqrt(2.0 / 3.0) * lattice_constant};
        slab.atoms.push_back(ScmAtom{site, site, layer});
      }
    }
  }
  return slab;
}

// The ideal slab with every atom displaced by a smooth, deterministic pattern (≤ 0.25 bohr).
ScmSlab distorted_slab() {
  ScmSlab slab = ideal_slab();
  for (int k = 0; k < std::ssize(slab.atoms); ++k) {
    Vec3& p = slab.atoms[k].distorted;
    p.x += 0.25 * std::sin(1.3 * k + 0.1);
    p.y += 0.2 * std::cos(0.7 * k);
    p.z += 0.1 * std::sin(2.1 * k + 0.5);
  }
  return slab;
}

// A molecule: its centre in skewed lattice coordinates (u, v), height, bond length and angles.
struct Molecule {
  double u, v, z, r, theta, phi;
};

// Molecules near and between the atoms, with r below, inside and above the clamp ranges
// ([0.725, 1.847] and [0.8, 1.7]), upright and tilted, below and far above the top layer.
const std::vector<Molecule> molecules = {
    {0.0, 0.0, 1.5, 1.4, 0.0, 0.0},  {2.4, 0.0, 2.0, 0.6, 1.2, 0.4},   {1.6, 1.6, 1.0, 2.5, 2.0, 4.0},
    {3.3, 1.1, 3.5, 1.0, 0.7, 5.5},  {-1.2, 2.7, 0.8, 1.75, 2.9, 1.9}, {0.9, -2.2, 6.0, 1.4, 1.5, 3.1},
    {4.1, 3.7, -0.5, 1.2, 0.3, 2.2}, {2.0, 1.0, 8.5, 0.8, 1.0, 1.0},   {-3.0, -1.0, 2.6, 0.4, 2.4, 0.2},
    {1.2, 2.2, 1.9, 3.5, 1.1, 2.7},
};

// ΔV of a molecule, its (u, v) converted with the correction's own lattice → Cartesian map.
double delta_v(const ScmCorrection& correction, const Molecule& m) {
  const CartesianVector xy = correction.to_cartesian(m.u, m.v);
  return correction.value_hartree(xy.x, xy.y, m.z, m.r, m.theta, m.phi);
}

// The bit pattern of a double, for bitwise comparisons.
std::uint64_t bits(double x) { return std::bit_cast<std::uint64_t>(x); }

// The number of deviations that are not within the tolerance (NaN counts as one).
int count_failures(const std::vector<double>& deviations) {
  int failures = 0;
  for (double d : deviations)
    if (!(d <= tolerance)) ++failures;
  return failures;
}

// The largest deviation, NaN if any is NaN.
double worst(const std::vector<double>& deviations) {
  double result = 0.0;
  for (double d : deviations) result = nan_aware_max(result, d);
  return result;
}

// The H atoms of a molecule at Cartesian (x, y), computed here, independently of the library.
std::pair<Vec3, Vec3> hydrogen_atoms(double x, double y, const Molecule& m) {
  const double hx = 0.5 * m.r * std::sin(m.theta) * std::cos(m.phi);
  const double hy = 0.5 * m.r * std::sin(m.theta) * std::sin(m.phi);
  const double hz = 0.5 * m.r * std::cos(m.theta);
  return {Vec3{x - hx, y - hy, m.z - hz}, Vec3{x + hx, y + hy, m.z + hz}};
}

// f without the switch, written here for the mutated pair functions below.
double unswitched(const ScmPairTerms& t, double distance) {
  const double q = distance - t.reference_distance;
  const auto& c = t.polynomial;
  return -(c[0] + q * (c[1] + q * (c[2] + q * c[3]))) * std::exp(-t.decay * q);
}

// Porting mistake: the switch phase measured from p₇ instead of p₈ (ρ falls from 1 to 0).
double switch_from_p7(const ScmPairTerms& t, double distance) {
  if (distance > t.switch_end) return t.switched_value;
  if (!(distance > t.switch_start)) return unswitched(t, distance);
  const double rho = 0.5 * std::cos(((distance - t.switch_start) * pi) / (t.switch_end - t.switch_start)) + 0.5;
  return (1.0 - rho) * unswitched(t, distance) + rho * t.switched_value;
}

// Porting mistake: zero beyond p₈ (the cached value left out, as the Fortran comment "f is left
// out in this version" might suggest).
double zero_beyond_p8(const ScmPairTerms& t, double distance) {
  return distance > t.switch_end ? 0.0 : scm_pair_energy(t, distance);
}

// A pair function: the library's or a mutated one.
using PairEnergy = std::function<double(const ScmPairTerms&, double)>;

// The pair terms of both layers at bond lengths below, at, inside and above the clamp ranges.
std::vector<ScmPairTerms> pair_terms_samples() {
  const ScmParameters p = parameters();
  std::vector<ScmPairTerms> samples;
  for (const ScmLayerParameters& layer : p.layers)
    for (double r : {0.5, layer.r_clamp_min, 1.2, 1.4, layer.r_clamp_max, 3.0})
      samples.push_back(scm_pair_terms(layer, r));
  return samples;
}

}  // namespace

// Distorted = ideal, no rescaling: each atom adds (ΔV − f) + f, so ΔV is exactly +0, as for a slab without
// atoms (the start, the Fortran's potE = 0.0). Rescaling with s = 1 is not "none": its round trip rounds h.
TEST(ScmCorrection, UndistortedSlabWithoutRescalingGivesExactlyPlusZero) {
  const ScmCorrection correction(parameters(), ideal_slab(), std::nullopt, skewing_angle);
  for (const Molecule& m : molecules) EXPECT_EQ(bits(delta_v(correction, m)), bits(+0.0)) << m.u << " " << m.v;
  const ScmCorrection no_atoms(parameters(), ScmSlab{}, std::nullopt, skewing_angle);
  EXPECT_EQ(bits(delta_v(no_atoms, molecules[0])), bits(+0.0));

  const ScmCorrection round_trip(parameters(), ideal_slab(), 1.0, skewing_angle);
  int nonzero = 0;
  for (const Molecule& m : molecules)
    if (delta_v(round_trip, m) != 0.0) ++nonzero;
  EXPECT_GT(nonzero, 0);
}

// f(p₈) = C bitwise (ρ = 1 there), and f is continuous at p₇ and p₈: the jumps over one ulp are
// within the tolerance. Negative controls: a switch measured from p₇ jumps at p₇, and a pair
// function that is zero beyond p₈ jumps at p₈, for every sample.
TEST(ScmPairEnergy, ContinuousAtTheSwitchAndEqualToTheCachedValueAtP8) {
  const auto jumps = [](const PairEnergy& f, const ScmPairTerms& t) {
    const double inf = std::numeric_limits<double>::infinity();
    const double at_start = std::abs(f(t, std::nextafter(t.switch_start, inf)) - f(t, t.switch_start));
    const double at_end = std::abs(f(t, std::nextafter(t.switch_end, inf)) - f(t, std::nextafter(t.switch_end, -inf)));
    return std::pair{at_start, at_end};
  };
  std::vector<double> start, end, mutant_start, mutant_end;
  for (const ScmPairTerms& t : pair_terms_samples()) {
    EXPECT_EQ(bits(scm_pair_energy(t, t.switch_end)), bits(t.switched_value));
    EXPECT_NE(t.switched_value, 0.0);
    const auto [s, e] = jumps(scm_pair_energy, t);
    start.push_back(s);
    end.push_back(e);
    mutant_start.push_back(jumps(switch_from_p7, t).first);
    mutant_end.push_back(jumps(zero_beyond_p8, t).second);
  }
  std::cout << "worst jump at p7 " << worst(start) << ", at p8 " << worst(end) << " Ha; mutants " << worst(mutant_start)
            << ", " << worst(mutant_end) << "\n";
  EXPECT_EQ(count_failures(start), 0);
  EXPECT_EQ(count_failures(end), 0);
  EXPECT_EQ(count_failures(mutant_start), std::ssize(start));
  EXPECT_EQ(count_failures(mutant_end), std::ssize(end));
}

// The parameters use r clamped to [P17, P18]: below P17 they are those at P17, above P18 those at
// P18, bit for bit; inside they change with r (a step of one ulp is too small to change them).
TEST(ScmPairTerms, ClampRAtP17AndP18) {
  // The number of the 9 fields (p₁…p₈, C) in which two pair terms differ bitwise.
  const auto differing_fields = [](const ScmPairTerms& a, const ScmPairTerms& b) {
    const std::array<double, 9> x = {a.reference_distance, a.decay,         a.polynomial[0],
                                     a.polynomial[1],      a.polynomial[2], a.polynomial[3],
                                     a.switch_start,       a.switch_end,    a.switched_value};
    const std::array<double, 9> y = {b.reference_distance, b.decay,         b.polynomial[0],
                                     b.polynomial[1],      b.polynomial[2], b.polynomial[3],
                                     b.switch_start,       b.switch_end,    b.switched_value};
    int count = 0;
    for (int k = 0; k < 9; ++k)
      if (bits(x[k]) != bits(y[k])) ++count;
    return count;
  };
  const ScmParameters p = parameters();
  for (const ScmLayerParameters& layer : p.layers) {
    const ScmPairTerms at_min = scm_pair_terms(layer, layer.r_clamp_min);
    const ScmPairTerms at_max = scm_pair_terms(layer, layer.r_clamp_max);
    for (double r : {0.0, 0.3, layer.r_clamp_min - 1e-6})
      EXPECT_EQ(differing_fields(scm_pair_terms(layer, r), at_min), 0) << r;
    for (double r : {layer.r_clamp_max + 1e-6, 2.5, 10.0})
      EXPECT_EQ(differing_fields(scm_pair_terms(layer, r), at_max), 0) << r;
    EXPECT_GT(differing_fields(scm_pair_terms(layer, layer.r_clamp_min + 1e-6), at_min), 0);
    EXPECT_GT(differing_fields(scm_pair_terms(layer, layer.r_clamp_max - 1e-6), at_max), 0);
    EXPECT_GT(differing_fields(at_min, at_max), 0);
  }
}

// One displaced layer-2 atom, no rescaling: ΔV = Σ_a f(|h_a − R_distorted|) − f(|h_a − R_ideal|), h at the
// unclamped r. Controls: the layers swapped (passes only at Z = 8.5 bohr: beyond p₈ of both layers), and the
// r of h clamped too, as the CRP does (fails at the 5 molecules with r outside the clamp range).
TEST(ScmCorrection, OneDisplacedAtomGivesTheAnalyticDifference) {
  ScmSlab slab = ideal_slab();
  const int displaced = atom_index(2, 0, 0);  // under the hcp site next to the origin
  ASSERT_EQ(slab.atoms[displaced].layer, 2);
  slab.atoms[displaced].distorted.x += 0.3;
  slab.atoms[displaced].distorted.y -= 0.2;
  slab.atoms[displaced].distorted.z += 0.25;
  const ScmAtom& atom = slab.atoms[displaced];
  const ScmLayerParameters layer2 = parameters().layers[1];
  const ScmCorrection correction(parameters(), slab, std::nullopt, skewing_angle);
  const ScmCorrection layers_mixed_up(swapped_layers(), slab, std::nullopt, skewing_angle);

  std::vector<double> deviations, mixed_up, clamped_positions;
  int outside_clamp = 0;
  for (const Molecule& m : molecules) {
    const CartesianVector xy = correction.to_cartesian(m.u, m.v);
    const auto [h1, h2] = hydrogen_atoms(xy.x, xy.y, m);
    const ScmPairTerms terms = scm_pair_terms(layer2, m.r);
    const double expected =
        (scm_pair_energy(terms, distance(h1, atom.distorted)) - scm_pair_energy(terms, distance(h1, atom.ideal))) +
        (scm_pair_energy(terms, distance(h2, atom.distorted)) - scm_pair_energy(terms, distance(h2, atom.ideal)));
    deviations.push_back(std::abs(correction.value_hartree(xy.x, xy.y, m.z, m.r, m.theta, m.phi) - expected));
    mixed_up.push_back(std::abs(layers_mixed_up.value_hartree(xy.x, xy.y, m.z, m.r, m.theta, m.phi) - expected));
    const double clamped_r = std::clamp(m.r, layer2.r_clamp_min, layer2.r_clamp_max);
    if (clamped_r != m.r) ++outside_clamp;
    clamped_positions.push_back(
        std::abs(correction.value_hartree(xy.x, xy.y, m.z, clamped_r, m.theta, m.phi) - expected));
  }
  std::cout << "worst deviation " << worst(deviations) << " Ha; mutants: layers mixed up " << worst(mixed_up)
            << ", clamped positions " << worst(clamped_positions) << "\n";
  EXPECT_EQ(count_failures(deviations), 0);
  EXPECT_EQ(count_failures(mixed_up), std::ssize(molecules) - 1);
  EXPECT_EQ(outside_clamp, 5);
  EXPECT_EQ(count_failures(clamped_positions), outside_clamp);
}

// Moving the molecule by a₂ (v + a), the distorted atoms by its image T and the ideal ones by T/s (by T without
// rescaling) leaves ΔV unchanged to rounding. Controls: lattice coordinates passed as Cartesian ones, and, with
// rescaling, a correction that ignores it (compared with itself before the move).
TEST(ScmCorrection, InvariantUnderALatticeTranslationOfMoleculeAndBothAtomSets) {
  const ScmParameters p = parameters();
  const Vec3 shift{lattice_constant * std::cos(skewing_angle), lattice_constant * std::sin(skewing_angle), 0.0};
  for (const std::optional<double> s : {std::optional<double>{}, std::optional<double>{rescaling}}) {
    SCOPED_TRACE(s ? "rescaling" : "no rescaling");
    const double ideal_factor = s ? 1.0 / *s : 1.0;
    ScmSlab moved = distorted_slab();
    for (ScmAtom& atom : moved.atoms) {
      atom.distorted = Vec3{atom.distorted.x + shift.x, atom.distorted.y + shift.y, atom.distorted.z};
      atom.ideal = Vec3{atom.ideal.x + ideal_factor * shift.x, atom.ideal.y + ideal_factor * shift.y, atom.ideal.z};
    }
    const ScmCorrection original(p, distorted_slab(), s, skewing_angle);
    const ScmCorrection translated(p, moved, s, skewing_angle);
    const ScmCorrection ignored_before(p, distorted_slab(), std::nullopt, skewing_angle);
    const ScmCorrection ignored_after(p, moved, std::nullopt, skewing_angle);
    std::vector<double> deviations, not_cartesian, ignored;
    for (const Molecule& m : molecules) {
      Molecule m_moved = m;
      m_moved.v += lattice_constant;
      const double before = delta_v(original, m);
      deviations.push_back(std::abs(delta_v(translated, m_moved) - before));
      const double not_cartesian_before = original.value_hartree(m.u, m.v, m.z, m.r, m.theta, m.phi);
      not_cartesian.push_back(
          std::abs(translated.value_hartree(m_moved.u, m_moved.v, m.z, m.r, m.theta, m.phi) - not_cartesian_before));
      ignored.push_back(std::abs(delta_v(ignored_after, m_moved) - delta_v(ignored_before, m)));
    }
    std::cout << (s ? "rescaling" : "no rescaling") << ": worst deviation " << worst(deviations)
              << " Ha; mutants: not Cartesian " << worst(not_cartesian) << ", rescaling ignored " << worst(ignored)
              << "\n";
    EXPECT_EQ(count_failures(deviations), 0);
    EXPECT_EQ(count_failures(not_cartesian), std::ssize(molecules));
    EXPECT_EQ(count_failures(ignored), s ? std::ssize(molecules) : 0);
  }
}

// Exchanging the H atoms, (θ, φ) → (π − θ, φ + π), leaves ΔV unchanged to rounding. Control: θ and φ swapped on
// one side; it passes only at the upright molecule with φ = 0, whose swapped angles (0, 0) and (π, π) are again
// an exchanged pair.
TEST(ScmCorrection, SymmetricUnderExchangeOfTheTwoHAtoms) {
  for (const std::optional<double> s : {std::optional<double>{}, std::optional<double>{rescaling}}) {
    SCOPED_TRACE(s ? "rescaling" : "no rescaling");
    const ScmCorrection correction(parameters(), distorted_slab(), s, skewing_angle);
    std::vector<double> deviations, swapped;
    for (const Molecule& m : molecules) {
      const CartesianVector xy = correction.to_cartesian(m.u, m.v);
      const double theta = pi - m.theta;
      const double phi = m.phi + pi;
      deviations.push_back(std::abs(correction.value_hartree(xy.x, xy.y, m.z, m.r, m.theta, m.phi) -
                                    correction.value_hartree(xy.x, xy.y, m.z, m.r, theta, phi)));
      swapped.push_back(std::abs(correction.value_hartree(xy.x, xy.y, m.z, m.r, m.phi, m.theta) -
                                 correction.value_hartree(xy.x, xy.y, m.z, m.r, phi, theta)));
    }
    std::cout << (s ? "rescaling" : "no rescaling") << ": worst deviation " << worst(deviations)
              << " Ha; mutant (one-sided swap) " << worst(swapped) << "\n";
    EXPECT_EQ(count_failures(deviations), 0);
    EXPECT_EQ(count_failures(swapped), std::ssize(molecules) - 1);
  }
}

// The constructor refuses atoms outside layers 1–2 (their parameters do not exist) and a
// rescaling factor that is not positive.
TEST(ScmCorrection, RejectsAtomsOutsideLayers1And2AndANonPositiveFactor) {
  for (int layer : {0, 3}) {
    ScmSlab slab = ideal_slab();
    slab.atoms[5].layer = layer;
    EXPECT_THROW(ScmCorrection(parameters(), slab, std::nullopt, skewing_angle), std::invalid_argument) << layer;
  }
  for (double s : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN()})
    EXPECT_THROW(ScmCorrection(parameters(), ideal_slab(), s, skewing_angle), std::invalid_argument) << s;
}

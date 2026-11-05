// Surface point groups from generators. Fortran: CxvSurfaceSymmetry.f90 (generalised)
#include "symmetry/point_group.hpp"

#include <algorithm>
#include <numeric>

#include "core/error.hpp"
#include "core/real_compare.hpp"

namespace spodvr {

namespace {

constexpr SurfaceOperation rotation(int twelfths) { return {false, twelfths}; }
constexpr SurfaceOperation mirror(int twelfths) { return {true, twelfths}; }

bool canonical_order(const SurfaceOperation& a, const SurfaceOperation& b) {
  if (a.mirror != b.mirror) return !a.mirror;
  return a.angle_twelfths < b.angle_twelfths;
}

bool contains(const std::vector<SurfaceOperation>& list, const SurfaceOperation& g) {
  return std::find(list.begin(), list.end(), g) != list.end();
}

}  // namespace

std::vector<SurfaceOperation> group_closure(const std::vector<SurfaceOperation>& generators) {
  std::vector<SurfaceOperation> group = {identity_operation()};
  // Multiply every element by every generator until nothing new appears (at most 24 elements).
  for (std::size_t i = 0; i < group.size(); ++i)
    for (const SurfaceOperation& g : generators) {
      const SurfaceOperation product = compose(g, group[i]);
      if (!contains(group, product)) group.push_back(product);
    }
  std::sort(group.begin(), group.end(), canonical_order);
  return group;
}

bool is_group(const std::vector<SurfaceOperation>& operations) {
  if (!contains(operations, identity_operation())) return false;
  for (const SurfaceOperation& a : operations)
    for (const SurfaceOperation& b : operations)
      if (!contains(operations, compose(a, b))) return false;
  return true;
}

int phi_points_multiple(const std::vector<SurfaceOperation>& operations) {
  int common = 12;
  for (const SurfaceOperation& g : operations) common = std::gcd(common, g.angle_twelfths);
  return 12 / common;
}

PointGroup::PointGroup(std::string_view name) {
  using fortran_compare::seq;
  if (seq(name, "C2v")) {
    name_ = "C2v";
    generators_ = {rotation(6), mirror(0)};
  } else if (seq(name, "C4v")) {
    name_ = "C4v";
    generators_ = {rotation(3), mirror(0)};
  } else if (seq(name, "C6v")) {
    name_ = "C6v";
    generators_ = {rotation(2), mirror(0)};
  } else if (seq(name, "C3v")) {
    name_ = "C3v";
    generators_ = {rotation(4), mirror(6)};
  } else if (seq(name, "C3")) {
    name_ = "C3";
    generators_ = {rotation(4)};
  } else if (seq(name, "none")) {
    name_ = "none";
  } else {
    throw InputError("SystemID",
                     "unknown surface point group '" + std::string(name) + "' (known: C2v, C4v, C6v, C3v, C3, none)");
  }
  operations_ = group_closure(generators_);
}

}  // namespace spodvr

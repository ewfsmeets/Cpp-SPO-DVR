// The CRP PES of a data set (SystemID CRPlibrary): the adapter between the CRP library port
// (src/crplib/) and the Potential interface. Fortran: H2SurfaceSystem.f90 (CRPlibrary branch)
#include "potential/crp_potential.hpp"

#include <format>

#include "core/error.hpp"
#include "crplib/atom_surface.hpp"
#include "crplib/dataset.hpp"
#include "crplib/errors.hpp"
#include "crplib/qd_wrapper.hpp"
#include "crplib/units.hpp"

namespace spodvr {

namespace {

constexpr const char* data_set_key = "CRPlibDataSet";

std::string unit_name(crplib::LengthUnit unit) { return unit == crplib::LengthUnit::bohr ? "bohr" : "angstrom"; }

std::string unit_name(crplib::EnergyUnit unit) {
  switch (unit) {
    case crplib::EnergyUnit::hartree: return "hartree";
    case crplib::EnergyUnit::ev: return "eV";
    case crplib::EnergyUnit::kcal_per_mol: return "kcal/mol";
    case crplib::EnergyUnit::kj_per_mol: return "kJ/mol";
    case crplib::EnergyUnit::joule: return "J";
  }
  return "?";
}

crplib::DataSet read_checked(const std::filesystem::path& directory) {
  try {
    return crplib::read_data_set(directory);
  } catch (const crplib::DataError& e) {
    throw InputError(data_set_key, e.what());
  }
}

CrpDataSetSummary make_summary(const std::filesystem::path& directory, const crplib::DataSet& data,
                               const crplib::QdWrapper& wrapper) {
  const crplib::MoleculeSurface& surface = wrapper.surface();
  CrpDataSetSummary s;
  s.directory = directory.string();
  if (data.info) s.description = data.info->short_info;
  s.symmetry = crplib::symmetry_name(data.symmetry);
  s.file_units = unit_name(data.units.length) + ", " + unit_name(data.units.energy);
  s.surface_layers = data.num_layers();
  s.atom_cuts = static_cast<int>(data.atom_cuts.size());
  s.molecule_cuts = static_cast<int>(data.molecule_cuts.size());
  s.molecule_cuts_used = crplib::molecule_cut_count(data.symmetry);
  s.cuts_hold_interpolation = data.cuts_interpolation;
  s.r_min = surface.r_min();
  s.r_max = surface.r_max();
  s.z_min = surface.z_min();
  s.zso_min = surface.zso_min();
  s.zso_max = surface.zso_max();
  s.z_asymptotic = surface.z_asymptotic();
  s.energy_shift = -crplib::convert_energy(surface.h2_asymptote(), crplib::EnergyUnit::hartree, crplib::EnergyUnit::ev);
  const crplib::ClusterSize cluster = surface.atom_surface().cluster_size();
  s.cluster_atoms = cluster.atoms;
  s.cluster_atoms_within_reach = cluster.within_reach;
  s.warnings = surface.warnings();
  return s;
}

}  // namespace

CrpPotential::CrpPotential(const std::filesystem::path& data_set) {
  const crplib::DataSet data = read_checked(data_set);
  try {
    wrapper_ = std::make_unique<const crplib::QdWrapper>(data);
  } catch (const crplib::DataError& e) {
    throw InputError(data_set_key, e.what());
  } catch (const crplib::EvaluationError& e) {  // the setup evaluates the atom PES (reduced cuts)
    throw InputError(data_set_key, e.what());
  }
  const crplib::QdParameters p = wrapper_->parameters();
  parameters_ = {p.lattice_constant, p.lattice_constant, p.z_asymptotic, p.skewing_angle, p.point_group};
  summary_ = make_summary(data_set, data, *wrapper_);
}

CrpPotential::~CrpPotential() = default;

double CrpPotential::value(double x, double y, double z, double r, double theta, double phi, int /*iz*/,
                           int /*ir*/) const {
  try {
    return wrapper_->potential_ev(x, y, z, r, theta, phi);
  } catch (const crplib::EvaluationError& e) {
    throw InputError(data_set_key, std::format("the CRP PES cannot be evaluated at x = {}, y = {}, Z = {}, r = {} "
                                               "bohr, theta = {}, phi = {} rad (the grid is outside the data): {}",
                                               x, y, z, r, theta, phi, e.what()));
  }
}

CrpRangeReport CrpPotential::range_report(std::span<const double> z, std::span<const double> r) const {
  const crplib::RangeReport report = wrapper_->range_report(z, r);
  return {report.r_points, report.r_clamped, report.z_points, report.z_clamped, report.zr_points, report.zr_uncovered};
}

}  // namespace spodvr

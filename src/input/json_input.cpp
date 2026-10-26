#include "input/json_input.hpp"

#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <variant>

#include "core/error.hpp"
#include "core/real_compare.hpp"

namespace spodvr {

namespace {

using Json = nlohmann::json;

double get_real(const std::string& key, const Json& value) {
  if (!value.is_number()) throw InputError(key, "expected a real number");
  return value.get<double>();
}

int get_integer(const std::string& key, const Json& value) {
  if (!value.is_number_integer()) throw InputError(key, "expected an integer");
  const auto v = value.get<long long>();
  if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max())
    throw InputError(key, "integer out of range");
  return static_cast<int>(v);
}

bool get_logical(const std::string& key, const Json& value) {
  if (!value.is_boolean()) throw InputError(key, "expected true or false");
  return value.get<bool>();
}

std::string get_string(const std::string& key, const Json& value) {
  if (!value.is_string()) throw InputError(key, "expected a string");
  return value.get<std::string>();
}

// Member pointers of InputParameters, one alternative per value type.
using Field = std::variant<double InputParameters::*, int InputParameters::*,
                           bool InputParameters::*, std::string InputParameters::*,
                           std::optional<double> InputParameters::*,
                           std::optional<int> InputParameters::*>;

struct FieldSpec {
  Field field;
  bool required;  // Fortran: pre-registered (compulsory) field, plus TolSp
};

const std::map<std::string, FieldSpec>& field_table() {
  using P = InputParameters;
  static const std::map<std::string, FieldSpec> table = {
      // Compulsory fields (Fortran: InputNew.f90 ReadInput Set...Field calls).
      {"Title", {&P::title, true}},
      {"SystemID", {&P::system_id, true}},
      {"NewCalc", {&P::new_calculation, true}},
      {"ZDoF", {&P::z_dof, true}},
      {"RDoF", {&P::r_dof, true}},
      {"XDoF", {&P::x_dof, true}},
      {"YDoF", {&P::y_dof, true}},
      {"ThetaDoF", {&P::theta_dof, true}},
      {"PhiDoF", {&P::phi_dof, true}},
      {"Mass1", {&P::mass1, true}},
      {"Mass2", {&P::mass2, true}},
      {"HomoNuc", {&P::homonuclear, true}},
      {"ZStart", {&P::z_start, true}},
      {"dZ", {&P::dz, true}},
      {"nZ", {&P::nz, true}},
      {"nZSpecular", {&P::nz_specular, true}},
      {"RStart", {&P::r_start, true}},
      {"dR", {&P::dr, true}},
      {"nR", {&P::nr, true}},
      {"nX", {&P::nx, true}},
      {"nY", {&P::ny, true}},
      {"JMax", {&P::j_max, true}},
      {"ProjOp", {&P::projection_operator, true}},
      {"tStart", {&P::t_start, true}},
      {"tEnd", {&P::t_end, true}},
      {"dtStep", {&P::time_step, true}},
      {"nStep", {&P::n_steps, true}},
      {"dtAna", {&P::dt_analysis, true}},
      {"DynamZSpAna", {&P::dynamical_zsp_analysis, true}},
      {"dtAve", {&P::dt_averages, true}},
      {"Zana", {&P::z_analysis, true}},
      {"NrOfEs", {&P::n_energies, true}},
      {"EMinWP", {&P::e_min_wave_packet, true}},
      {"EMaxWP", {&P::e_max_wave_packet, true}},
      {"Z0", {&P::z0, true}},
      {"ZOptStart", {&P::z_opt_start, true}},
      {"ZOptEnd", {&P::z_opt_end, true}},
      {"ROptStart", {&P::r_opt_start, true}},
      {"ROptEnd", {&P::r_opt_end, true}},
      {"ZSpOptStart", {&P::zsp_opt_start, true}},
      {"ZSpOptEnd", {&P::zsp_opt_end, true}},
      {"EParallel", {&P::e_parallel, true}},
      {"PhiParallel", {&P::phi_parallel, true}},
      {"writeSMatrix", {&P::write_s_matrix, true}},
      {"writeProbs", {&P::write_probabilities, true}},
      {"writeRDProbs", {&P::write_rd_probabilities, true}},
      // TolSp has no default in the Fortran code (uninitialised if absent): required here.
      {"TolSp", {&P::tol_specular, true}},

      // Optional fields with Fortran defaults.
      {"MJMax", {&P::mj_max, false}},
      {"onlyAnalysis", {&P::only_analysis, false}},
      {"alsoAnalysis", {&P::also_analysis, false}},
      {"JAna", {&P::j_analysis, false}},
      {"VAna", {&P::v_analysis, false}},
      {"DifOrdAna", {&P::diffraction_order_analysis, false}},
      {"ZAnaSp", {&P::z_analysis_specular, false}},
      {"doFluxAnalysis", {&P::do_flux_analysis, false}},
      {"rFluxAnalysis", {&P::r_flux_analysis, false}},
      {"EMinFluxAna", {&P::e_min_flux_analysis, false}},
      {"EMaxFluxAna", {&P::e_max_flux_analysis, false}},
      {"DeltaEFluxAna", {&P::delta_e_flux_analysis, false}},
      {"VIn", {&P::v_in, false}},
      {"JIn", {&P::j_in, false}},
      {"MJIn", {&P::mj_in, false}},
      {"Potential", {&P::only_potential, false}},
      {"ImposeVMin", {&P::impose_v_min, false}},
      {"VMin", {&P::v_min, false}},
      {"ImposeVMax", {&P::impose_v_max, false}},
      {"VMax", {&P::v_max, false}},
      {"UseSwitching", {&P::use_switching, false}},
      {"ZminSwitch", {&P::z_min_switch, false}},
      {"ZmaxSwitch", {&P::z_max_switch, false}},
      {"ZOptEnergy", {&P::z_opt_energy, false}},
      {"ZOptFactor", {&P::z_opt_factor, false}},
      {"ROptEnergy", {&P::r_opt_energy, false}},
      {"ROptFactor", {&P::r_opt_factor, false}},
      {"ZSpOptEnergy", {&P::zsp_opt_energy, false}},
      {"ZSpOptFactor", {&P::zsp_opt_factor, false}},
      {"DoDump", {&P::do_dump, false}},
      {"FullSMatrix", {&P::full_s_matrix, false}},
  };
  return table;
}

// Fortran input fields of features that are not available in this program. Logical keys are
// accepted when false (the Fortran default); all others are rejected whenever present.
const std::set<std::string>& unsupported_logical_keys() {
  static const std::set<std::string> keys = {"SystemMix", "UseSymmetry", "nPhisOptimal",
                                             "UseRCutOff", "write2DPES"};
  return keys;
}
const std::set<std::string>& unsupported_keys() {
  static const std::set<std::string> keys = {"SystemID2", "SystemMixCoeff", "RCutOff", "x2DPES",
                                             "y2DPES",    "theta2DPES",     "phi2DPES", "PESMixingAlfa",
                                             "CRPlibDataSet"};
  return keys;
}

void set_field(InputParameters& p, const std::string& key, const Field& field, const Json& value) {
  std::visit(
      [&](auto member) {
        using Member = std::remove_reference_t<decltype(p.*member)>;
        if constexpr (std::is_same_v<Member, double>) {
          p.*member = get_real(key, value);
        } else if constexpr (std::is_same_v<Member, int>) {
          p.*member = get_integer(key, value);
        } else if constexpr (std::is_same_v<Member, bool>) {
          p.*member = get_logical(key, value);
        } else if constexpr (std::is_same_v<Member, std::string>) {
          p.*member = get_string(key, value);
        } else if constexpr (std::is_same_v<Member, std::optional<int>>) {
          p.*member = get_integer(key, value);
        } else {
          p.*member = get_real(key, value);
        }
      },
      field);
}

void require_present(const Json& json, const std::string& key, const std::string& message) {
  if (!json.contains(key)) throw InputError(key, message);
}

}  // namespace

InputParameters parse_input(const Json& json) {
  if (!json.is_object()) throw InputError("", "the input must be a JSON object");

  InputParameters p;
  const auto& table = field_table();
  for (const auto& [key, value] : json.items()) {
    if (unsupported_logical_keys().contains(key)) {
      if (get_logical(key, value))
        throw NotSupportedError(key + " = true is not supported by this program");
      continue;
    }
    if (unsupported_keys().contains(key))
      throw NotSupportedError(key + " is not supported by this program");
    const auto it = table.find(key);
    if (it == table.end()) throw InputError(key, "unknown input key");
    set_field(p, key, it->second.field, value);
  }

  for (const auto& [key, spec] : table) {
    if (spec.required && !json.contains(key))
      throw InputError(key, "required input key is missing");
  }

  // Fortran: InputNew.f90 ReadInput ERROR calls.
  if (p.projection_operator == 0)
    throw InputError("ProjOp", "ProjOp = 0 is not supported anymore.");
  if (p.do_flux_analysis) {
    for (const char* key : {"rFluxAnalysis", "EMinFluxAna", "EMaxFluxAna", "DeltaEFluxAna"})
      require_present(json, key, std::string(key) + " must be present");
  }
  // One of energy / factor is required for every CAP. Fortran reads the factor only when the
  // energy is absent; specifying both is ambiguous and rejected here.
  const std::pair<const char*, const char*> cap_keys[] = {
      {"ZOptEnergy", "ZOptFactor"}, {"ROptEnergy", "ROptFactor"}, {"ZSpOptEnergy", "ZSpOptFactor"}};
  for (const auto& [energy, factor] : cap_keys) {
    const bool has_energy = json.contains(energy);
    const bool has_factor = json.contains(factor);
    if (!has_energy && !has_factor)
      throw InputError(energy, std::string(factor) + " or " + energy + " should be present.");
    if (has_energy && has_factor)
      throw InputError(energy, std::string("give either ") + energy + " or " + factor + ", not both");
  }
  using fortran_compare::rlt;
  if (p.z_opt_energy && rlt(*p.z_opt_energy, 0.0)) throw InputError("ZOptEnergy", "ZOptEnergy < 0.");
  if (p.r_opt_energy && rlt(*p.r_opt_energy, 0.0)) throw InputError("ROptEnergy", "ROptEnergy < 0.");
  if (p.zsp_opt_energy && rlt(*p.zsp_opt_energy, 0.0))
    throw InputError("ZSpOptEnergy", "ZSpOptEnergy < 0.");
  if (p.z_analysis_specular) {
    const double z_ana_sp = *p.z_analysis_specular;
    if (p.dynamical_zsp_analysis)
      throw InputError("ZAnaSp", "ZAnaSp can only be used with DynamZSpAna set to FALSE");
    if (z_ana_sp > p.zsp_opt_start) throw InputError("ZAnaSp", "ZAnaSp was set > ZSpOptStart");
    if (z_ana_sp < p.z0) throw InputError("ZAnaSp", "ZAnaSp was set < Z0");
  }
  // A negative MJMax means: equal to JMax. Fortran: InputNew.f90 ReadInput
  if (p.mj_max < 0) p.mj_max = p.j_max;
  return p;
}

InputParameters read_input_file(const std::filesystem::path& path) {
  std::ifstream file(path);
  if (!file) throw InputError("", "cannot open input file '" + path.string() + "'");
  Json json;
  try {
    json = Json::parse(file, nullptr, true, /*ignore_comments=*/true);
  } catch (const Json::parse_error& e) {
    throw InputError("", "invalid JSON in '" + path.string() + "': " + e.what());
  }
  return parse_input(json);
}

}  // namespace spodvr

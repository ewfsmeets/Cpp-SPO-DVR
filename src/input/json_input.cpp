// Reads the JSON input file into InputParameters and performs the checks that the Fortran
// reader does while reading. Fortran: InputNew.f90 ReadInput
#include "input/json_input.hpp"

#include <fstream>
#include <limits>
#include <set>
#include <string>

#include "core/error.hpp"
#include "core/real_compare.hpp"

namespace spodvr {

namespace {

using Json = nlohmann::json;

enum class Presence { Required, Optional };

// Reads typed values from the JSON object and remembers which keys were used, so that unknown
// keys can be reported afterwards. An optional key that is absent leaves the target unchanged
// (i.e. at its Fortran default).
class JsonReader {
 public:
  explicit JsonReader(const Json& json) : json_(json) {}

  void read(const std::string& key, double& target, Presence presence) {
    if (present(key, presence)) target = as_real(key, json_.at(key));
  }
  void read(const std::string& key, int& target, Presence presence) {
    if (present(key, presence)) target = as_integer(key, json_.at(key));
  }
  void read(const std::string& key, bool& target, Presence presence) {
    if (present(key, presence)) target = as_logical(key, json_.at(key));
  }
  void read(const std::string& key, std::string& target, Presence presence) {
    if (present(key, presence)) target = as_string(key, json_.at(key));
  }
  void read(const std::string& key, std::optional<double>& target) {
    if (present(key, Presence::Optional)) target = as_real(key, json_.at(key));
  }
  void read(const std::string& key, std::optional<int>& target) {
    if (present(key, Presence::Optional)) target = as_integer(key, json_.at(key));
  }

  // A key of a Fortran feature this program does not support. Logical keys are accepted when
  // false (the Fortran default, i.e. the feature is off); other keys are rejected if present.
  void reject_logical_if_true(const std::string& key) {
    bool enabled = false;
    read(key, enabled, Presence::Optional);
    if (enabled) throw NotSupportedError(key + " = true is not supported by this program");
  }
  void reject_if_present(const std::string& key) {
    if (json_.contains(key)) throw NotSupportedError(key + " is not supported by this program");
  }

  [[nodiscard]] bool contains(const std::string& key) const { return json_.contains(key); }

  void check_no_unknown_keys() const {
    for (const auto& item : json_.items())
      if (!used_keys_.contains(item.key())) throw InputError(item.key(), "unknown input key");
  }

 private:
  // True if the key is present; a missing required key is an error.
  bool present(const std::string& key, Presence presence) {
    used_keys_.insert(key);
    if (json_.contains(key)) return true;
    if (presence == Presence::Required) throw InputError(key, "required input key is missing");
    return false;
  }

  static double as_real(const std::string& key, const Json& value) {
    if (!value.is_number()) throw InputError(key, "expected a real number");
    return value.get<double>();
  }
  static int as_integer(const std::string& key, const Json& value) {
    if (!value.is_number_integer()) throw InputError(key, "expected an integer");
    const auto v = value.get<long long>();
    if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max())
      throw InputError(key, "integer out of range");
    return static_cast<int>(v);
  }
  static bool as_logical(const std::string& key, const Json& value) {
    if (!value.is_boolean()) throw InputError(key, "expected true or false");
    return value.get<bool>();
  }
  static std::string as_string(const std::string& key, const Json& value) {
    if (!value.is_string()) throw InputError(key, "expected a string");
    return value.get<std::string>();
  }

  const Json& json_;
  std::set<std::string> used_keys_;
};

// Compulsory keys: the fields the Fortran code pre-registers (Set...Field in ReadInput), plus TolSp,
// which has no Fortran default.
void read_required_keys(JsonReader& in, InputParameters& p) {
  constexpr auto required = Presence::Required;
  in.read("Title", p.title, required);
  in.read("SystemID", p.system_id, required);
  in.read("NewCalc", p.new_calculation, required);

  in.read("ZDoF", p.z_dof, required);
  in.read("RDoF", p.r_dof, required);
  in.read("XDoF", p.x_dof, required);
  in.read("YDoF", p.y_dof, required);
  in.read("ThetaDoF", p.theta_dof, required);
  in.read("PhiDoF", p.phi_dof, required);

  in.read("Mass1", p.mass1, required);
  in.read("Mass2", p.mass2, required);
  in.read("HomoNuc", p.homonuclear, required);

  in.read("ZStart", p.z_start, required);
  in.read("dZ", p.dz, required);
  in.read("nZ", p.nz, required);
  in.read("nZSpecular", p.nz_specular, required);
  in.read("RStart", p.r_start, required);
  in.read("dR", p.dr, required);
  in.read("nR", p.nr, required);
  in.read("nX", p.nx, required);
  in.read("nY", p.ny, required);
  in.read("JMax", p.j_max, required);
  in.read("ProjOp", p.projection_operator, required);

  in.read("tStart", p.t_start, required);
  in.read("tEnd", p.t_end, required);
  in.read("dtStep", p.time_step, required);
  in.read("nStep", p.n_steps, required);

  in.read("dtAna", p.dt_analysis, required);
  in.read("DynamZSpAna", p.dynamical_zsp_analysis, required);
  in.read("dtAve", p.dt_averages, required);
  in.read("Zana", p.z_analysis, required);
  in.read("NrOfEs", p.n_energies, required);

  in.read("EMinWP", p.e_min_wave_packet, required);
  in.read("EMaxWP", p.e_max_wave_packet, required);
  in.read("Z0", p.z0, required);

  in.read("ZOptStart", p.z_opt_start, required);
  in.read("ZOptEnd", p.z_opt_end, required);
  in.read("ROptStart", p.r_opt_start, required);
  in.read("ROptEnd", p.r_opt_end, required);
  in.read("ZSpOptStart", p.zsp_opt_start, required);
  in.read("ZSpOptEnd", p.zsp_opt_end, required);

  in.read("EParallel", p.e_parallel, required);
  in.read("PhiParallel", p.phi_parallel, required);

  in.read("writeSMatrix", p.write_s_matrix, required);
  in.read("writeProbs", p.write_probabilities, required);
  in.read("writeRDProbs", p.write_rd_probabilities, required);

  in.read("TolSp", p.tol_specular, required);
}

// Optional keys; absent keys keep the Fortran defaults of InputParameters.
void read_optional_keys(JsonReader& in, InputParameters& p) {
  constexpr auto optional = Presence::Optional;
  in.read("MJMax", p.mj_max, optional);

  in.read("onlyAnalysis", p.only_analysis, optional);
  in.read("alsoAnalysis", p.also_analysis, optional);
  in.read("JAna", p.j_analysis);
  in.read("VAna", p.v_analysis);
  in.read("DifOrdAna", p.diffraction_order_analysis);
  in.read("ZAnaSp", p.z_analysis_specular);

  in.read("doFluxAnalysis", p.do_flux_analysis, optional);
  in.read("rFluxAnalysis", p.r_flux_analysis, optional);
  in.read("EMinFluxAna", p.e_min_flux_analysis, optional);
  in.read("EMaxFluxAna", p.e_max_flux_analysis, optional);
  in.read("DeltaEFluxAna", p.delta_e_flux_analysis, optional);

  in.read("VIn", p.v_in, optional);
  in.read("JIn", p.j_in, optional);
  in.read("MJIn", p.mj_in, optional);

  in.read("Potential", p.only_potential, optional);
  in.read("ImposeVMin", p.impose_v_min, optional);
  in.read("VMin", p.v_min, optional);
  in.read("ImposeVMax", p.impose_v_max, optional);
  in.read("VMax", p.v_max, optional);
  in.read("UseSwitching", p.use_switching, optional);
  in.read("ZminSwitch", p.z_min_switch, optional);
  in.read("ZmaxSwitch", p.z_max_switch, optional);

  in.read("ZOptEnergy", p.z_opt_energy);
  in.read("ZOptFactor", p.z_opt_factor);
  in.read("ROptEnergy", p.r_opt_energy);
  in.read("ROptFactor", p.r_opt_factor);
  in.read("ZSpOptEnergy", p.zsp_opt_energy);
  in.read("ZSpOptFactor", p.zsp_opt_factor);

  in.read("DoDump", p.do_dump, optional);
  in.read("FullSMatrix", p.full_s_matrix, optional);
}

// Fortran input fields of features that are not available in this program.
void reject_unsupported_keys(JsonReader& in) {
  for (const char* key : {"SystemMix", "UseSymmetry", "nPhisOptimal", "UseRCutOff", "write2DPES"})
    in.reject_logical_if_true(key);
  for (const char* key : {"SystemID2", "SystemMixCoeff", "RCutOff", "x2DPES", "y2DPES", "theta2DPES",
                          "phi2DPES", "PESMixingAlfa", "CRPlibDataSet"})
    in.reject_if_present(key);
}

// The ERROR calls of the Fortran reader. Fortran: InputNew.f90 ReadInput
void check_while_reading(const JsonReader& in, const InputParameters& p) {
  using fortran_compare::rlt;

  if (p.projection_operator == 0) throw InputError("ProjOp", "ProjOp = 0 is not supported anymore.");

  if (p.do_flux_analysis) {
    for (const char* key : {"rFluxAnalysis", "EMinFluxAna", "EMaxFluxAna", "DeltaEFluxAna"})
      if (!in.contains(key)) throw InputError(key, std::string(key) + " must be present");
  }

  // Every optical potential needs its energy or its factor. Fortran reads the factor only when
  // the energy is absent; giving both is ambiguous and rejected here.
  const std::pair<std::string, std::string> cap_keys[] = {
      {"ZOptEnergy", "ZOptFactor"}, {"ROptEnergy", "ROptFactor"}, {"ZSpOptEnergy", "ZSpOptFactor"}};
  for (const auto& [energy, factor] : cap_keys) {
    const bool has_energy = in.contains(energy);
    const bool has_factor = in.contains(factor);
    if (!has_energy && !has_factor) throw InputError(energy, factor + " or " + energy + " should be present.");
    if (has_energy && has_factor) throw InputError(energy, "give either " + energy + " or " + factor + ", not both");
  }
  if (p.z_opt_energy && rlt(*p.z_opt_energy, 0.0)) throw InputError("ZOptEnergy", "ZOptEnergy < 0.");
  if (p.r_opt_energy && rlt(*p.r_opt_energy, 0.0)) throw InputError("ROptEnergy", "ROptEnergy < 0.");
  if (p.zsp_opt_energy && rlt(*p.zsp_opt_energy, 0.0)) throw InputError("ZSpOptEnergy", "ZSpOptEnergy < 0.");

  if (p.z_analysis_specular) {
    const double z_ana_sp = *p.z_analysis_specular;
    if (p.dynamical_zsp_analysis)
      throw InputError("ZAnaSp", "ZAnaSp can only be used with DynamZSpAna set to FALSE");
    if (z_ana_sp > p.zsp_opt_start) throw InputError("ZAnaSp", "ZAnaSp was set > ZSpOptStart");
    if (z_ana_sp < p.z0) throw InputError("ZAnaSp", "ZAnaSp was set < Z0");
  }
}

}  // namespace

InputParameters parse_input(const Json& json) {
  if (!json.is_object()) throw InputError("", "the input must be a JSON object");

  JsonReader in(json);
  InputParameters p;
  reject_unsupported_keys(in);
  read_required_keys(in, p);
  read_optional_keys(in, p);
  in.check_no_unknown_keys();
  check_while_reading(in, p);

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

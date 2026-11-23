// Tests: JSON reading, every input consistency check (valid and violating case), derived parameters.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "input/validate_outputs.hpp"
#include "input/wave_packet_params.hpp"
#include "potential/potential_factory.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using Json = nlohmann::json;

namespace {

Json reference_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return Json::parse(file);
}

// The input of the reference case `name`; a CRP data-set path (relative to the repository root) is made
// absolute.
Json case_json(const std::string& name) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + name + "/input.json");
  Json json = Json::parse(file, nullptr, true, /*ignore_comments=*/true);
  if (json.contains("CRPlibDataSet")) {
    const std::filesystem::path root = std::filesystem::path(SPODVR_REFERENCE_DIR).parent_path().parent_path();
    json["CRPlibDataSet"] = (root / json["CRPlibDataSet"].get<std::string>()).string();
  }
  return json;
}

Config build(const Json& json) { return build_config(parse_input(json)); }

// A PlotSites list of (name, x, y) entries.
struct SiteSpec {
  std::string name;
  double x;
  double y;
};
Json site_list(const std::vector<SiteSpec>& sites) {
  Json list = Json::array();
  for (const SiteSpec& site : sites) list.push_back({{"name", site.name}, {"x", site.x}, {"y", site.y}});
  return list;
}

// Valid energy-resolved keys: two sites a, b of PlotSites and one energy inside the packet's range.
void energy_resolved(Json& j) {
  j["PlotSites"] = site_list({{"a", 0.0, 0.0}, {"b", 0.5, 0.0}});
  j["EnergyResolvedSites"] = {"a", "b"};
  j["EnergyResolvedEnergies"] = {0.3};
}

// Expects an InputError naming `parameter` whose message contains `message`, the fixed text of the
// check that must fire: several checks share a parameter, so the message shows that the input reached
// its own check and not an earlier one. Returns the message ("" if none).
std::string expect_input_error(const Json& json, const std::string& parameter, const std::string& message) {
  if (message.empty()) {  // an empty fragment would match any message
    ADD_FAILURE() << "no expected message for " << parameter;
    return "";
  }
  try {
    const Config config = build(json);
    (void)make_spatial_grids(config);
    ADD_FAILURE() << "expected an InputError for " << parameter << " with input " << json.dump();
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), parameter) << e.what();
    EXPECT_NE(std::string(e.what()).find(message), std::string::npos)
        << "expected the check '" << message << "', got: " << e.what();
    return e.what();
  }
  return "";
}

}  // namespace

TEST(Input, ReferenceInputIsValid) {
  const Config config = build(reference_json());
  EXPECT_EQ(config.input.nz, 176);
  EXPECT_EQ(config.input.mj_max, 4);
  EXPECT_EQ(config.derived.mj_max_basis, 4);
  EXPECT_EQ(config.derived.run_mode, 0);
  EXPECT_EQ(config.derived.steps_per_analysis, 5);
  EXPECT_EQ(config.derived.analyses_per_average, 10);
  EXPECT_EQ(config.derived.z_asymptotic, 7.0);
}

TEST(Input, UnknownKeyIsAnError) {
  Json json = reference_json();
  json["NotAKey"] = 1;
  expect_input_error(json, "NotAKey", "unknown input key");
}

// A key given twice in one JSON object is an error (nlohmann alone keeps the last value; the
// Fortran reader stops on a repeated field too). Inline texts: a tracked .json file with a repeated
// key would be refused by the check-json hook.
namespace {

Json parse_text(const std::string& text) {
  std::istringstream in(text);
  return parse_input_text(in, "test.json");
}

// Expects parse_input_text to refuse `text` with an InputError naming `parameter`; returns its message.
std::string expect_repeated_key_error(const std::string& text, const std::string& parameter) {
  try {
    (void)parse_text(text);
    ADD_FAILURE() << "expected an InputError for " << text;
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), parameter) << e.what();
    return e.what();
  }
  return "";
}

}  // namespace

TEST(Input, DuplicateKeyIsAnError) {
  const std::string message = expect_repeated_key_error(R"({"nZ": 176, "dZ": 0.08, "nZ": 192})", "nZ");
  EXPECT_NE(message.find("given twice"), std::string::npos) << message;
  EXPECT_NE(message.find("test.json"), std::string::npos) << message;
  // The keys of the top-level object are kept while a nested object is read.
  expect_repeated_key_error(R"({"a": 1, "b": {}, "a": 2})", "a");
}

TEST(Input, DuplicateKeyInsideAnEntryNamesTheTopLevelKey) {
  const std::string message =
      expect_repeated_key_error(R"({"PlotSites": [{"name": "a", "x": 0.0, "x": 0.5, "y": 0.0}]})", "PlotSites");
  EXPECT_NE(message.find("'x'"), std::string::npos) << message;
  expect_repeated_key_error(R"({"a": 1, "PlotSites": [{"x": 0.0, "x": 0.5}]})", "PlotSites");
  expect_repeated_key_error(R"({"a": {"b": 1, "b": 2}})", "a");
  expect_repeated_key_error(R"({"a": [[{"b": 1, "b": 2}]]})", "a");
}

TEST(Input, SameKeyInDifferentObjectsIsValid) {
  EXPECT_EQ(parse_text(R"({"a": {"c": 1}, "c": 2})")["c"], 2);
  EXPECT_EQ(parse_text(R"({"a": [{"b": 1}, {"b": 2}]})")["a"][1]["b"], 2);
}

TEST(Input, DuplicateKeyCheckSkipsComments) {
  // A commented-out line is not a key; a real repeat in a file with comments still is.
  EXPECT_EQ(parse_text("{\n  // \"nZ\": 128,\n  \"nZ\": 176\n}")["nZ"], 176);
  expect_repeated_key_error("{\n  // the grid\n  \"nZ\": 176,\n  /* again */ \"nZ\": 192\n}", "nZ");
}

TEST(Input, ReadInputFileRefusesADuplicateKey) {
  const test_support::TemporaryDirectory directory("repeated_key");
  std::string text = reference_json().dump(2);
  text.insert(text.find('{') + 1, "\n  \"nZ\": 192,");
  const std::filesystem::path file = directory.write("input.json", text);
  try {
    (void)read_input_file(file);
    ADD_FAILURE() << "a repeated nZ was accepted";
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), "nZ") << e.what();
    EXPECT_NE(std::string(e.what()).find(file.string()), std::string::npos) << e.what();
  }
}

TEST(Input, MissingRequiredKeyIsAnError) {
  for (const std::string key : {"Title", "nZ", "TolSp", "writeRDProbs", "ZSpOptEnd"}) {
    Json json = reference_json();
    json.erase(key);
    expect_input_error(json, key, "required input key is missing");
  }
}

TEST(Input, WrongTypesAreErrors) {
  Json json = reference_json();
  json["nZ"] = 176.0;
  expect_input_error(json, "nZ", "expected an integer");
  json = reference_json();
  json["NewCalc"] = 1;
  expect_input_error(json, "NewCalc", "expected true or false");
  json = reference_json();
  json["dZ"] = "0.08";
  expect_input_error(json, "dZ", "expected a real number");
}

TEST(Input, IntegerAcceptedForRealKey) {
  Json json = reference_json();
  json["ZStart"] = -2;
  EXPECT_EQ(build(json).input.z_start, -2.0);
}

TEST(Input, UnsupportedFeatures) {
  for (const std::string key : {"SystemMix", "UseRCutOff"}) {
    Json json = reference_json();
    json[key] = false;
    EXPECT_NO_THROW((void)build(json)) << key;
    json[key] = true;
    EXPECT_THROW((void)build(json), NotSupportedError) << key;
  }
  Json json = reference_json();
  json["RCutOff"] = 3.0;
  EXPECT_THROW((void)build(json), NotSupportedError);
}

// dtFrames needs dtAna / dtStep to be whole as the program counts it: it truncates (Fortran INT), and 2.4 / 0.8 is
// 2.9999999999999996, so the analysis step would be 2 time steps (1.6 a.u.), not 3. Without dtFrames the input is
// accepted as the Fortran accepts it; with dtFrames it is refused, naming dtAna.
TEST(Input, DtFramesNeedsAWholeNumberOfStepsPerAnalysis) {
  Json json = reference_json();
  json["dtStep"] = 0.8;
  json["nStep"] = 18000;
  json["tEnd"] = 14400.0;
  json["dtAna"] = 2.4;
  json["dtAve"] = 24.0;
  EXPECT_EQ(build(json).derived.steps_per_analysis, 2);
  json["dtFrames"] = 24.0;
  expect_input_error(json, "dtAna", "dtAna / dtStep is not a whole number with dtFrames (the program truncates it");
  json.erase("dtFrames");
  json["NormBudget"] = true;  // its records follow the analysis steps too
  expect_input_error(json, "dtAna", "dtAna / dtStep is not a whole number with NormBudget");
}

// write2DPES and its coordinates (Fortran defaults false, 0, 0, 90, 0). The coordinates are checked only with
// write2DPES, with the Fortran's tolerant .RLT./.RGT. (1e-12): a value 5e-13 beyond a bound passes, 2e-12 does
// not. The bounds of x and y are the PES's lattice constants (hard wall: 6.239660, 5.239660 bohr).
TEST(Input, Write2DPESKeys) {
  const InputParameters defaults = build(reference_json()).input;
  EXPECT_FALSE(defaults.write_2d_pes);
  EXPECT_EQ(defaults.x_2d_pes, 0.0);
  EXPECT_EQ(defaults.y_2d_pes, 0.0);
  EXPECT_EQ(defaults.theta_2d_pes, 90.0);
  EXPECT_EQ(defaults.phi_2d_pes, 0.0);
  Json json = reference_json();
  json["x2DPES"] = -1.0;  // not checked without write2DPES (the Fortran checks only inside "if (write2DPES)")
  json["phi2DPES"] = 400.0;
  EXPECT_NO_THROW((void)build(json));
  json = reference_json();
  json["write2DPES"] = true;
  json["x2DPES"] = 6.239660 + 5e-13;
  json["y2DPES"] = -5e-13;
  json["theta2DPES"] = 360.0;
  json["phi2DPES"] = 0.0;
  const InputParameters p = build(json).input;
  EXPECT_TRUE(p.write_2d_pes);
  EXPECT_EQ(p.x_2d_pes, 6.239660 + 5e-13);
  EXPECT_EQ(p.theta_2d_pes, 360.0);
  json["x2DPES"] = 6.239660 + 2e-12;
  expect_input_error(json, "x2DPES", "x2DPES > alat.");
  json["x2DPES"] = "1.0";
  EXPECT_THROW((void)build(json), InputError);
}

// CRPlibDataSet (directory of the CRP data set, relative to the working directory): required with
// SystemID CRPlibrary, an error with any other PES. build_config reads the data set (the PES gives
// the system parameters), so a missing one is reported there; the data-set checks themselves are
// in tests/crplib and test_crp_potential.
TEST(Input, CrpDataSetKey) {
  Json json = reference_json();  // hard wall
  json["CRPlibDataSet"] = "some/data/set";
  expect_input_error(json, "CRPlibDataSet", "is only used with SystemID CRPlibrary or CRPDistrtSurface");

  json = reference_json();
  json["SystemID"] = "CRPlibrary";
  expect_input_error(json, "CRPlibDataSet", "SystemID CRPlibrary needs CRPlibDataSet");  // missing
  json["CRPlibDataSet"] = "";
  expect_input_error(json, "CRPlibDataSet", "SystemID CRPlibrary needs CRPlibDataSet");
  json["CRPlibDataSet"] = "no/such/data/set";
  expect_input_error(json, "CRPlibDataSet", "does not exist (no file CRPData)");

  const std::string root = std::filesystem::path(SPODVR_REFERENCE_DIR).parent_path().parent_path().string();
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/crp_c6v_synthetic_small/input.json");
  json = Json::parse(file);
  json["SystemID"] = "crplibrary";  // case-insensitive, like the Fortran
  json["CRPlibDataSet"] = root + "/" + json["CRPlibDataSet"].get<std::string>();
  const Config config = build(json);
  EXPECT_EQ(config.input.crp_dataset, json["CRPlibDataSet"].get<std::string>());
  EXPECT_EQ(config.derived.system.point_group, "C6v");
  EXPECT_EQ(config.derived.z_asymptotic_unswitched, 10.0);
}

// ---------------------------------------------------------------------------------------------
// The SCM PES (SystemID CRPDistrtSurface): its keys and checks, each with a valid and a violating
// case. The potential is built in build_config, so the file, slab number, cut-off and lattice
// checks of the potential appear here as input errors too.
// ---------------------------------------------------------------------------------------------

namespace {

const std::filesystem::path example_directory =
    std::filesystem::path(SPODVR_REFERENCE_DIR).parent_path().parent_path() / "potentials/H2Cu111SRP_C3v.dir";
const std::string example_slab = (example_directory / "examplePosfile.dat").string();

const std::vector<std::string> scm_keys = {"DistortedPESMode", "PosFileSlabNum", "EnableSurfExpan",
                                           "SCMCutoff",        "SCMParFilepath", "PosFilepath"};

// The CRP gate case on the SRP48 data set (crp_cu111_srp48_small) as an SCM input with the example
// files of the data set (all keys given, paths absolute).
Json scm_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/crp_cu111_srp48_small/input.json");
  Json json = Json::parse(file);
  json["SystemID"] = "CRPDistrtSurface";
  json["CRPlibDataSet"] = (example_directory / "H2Cu111SRP_C3v").string();
  json["DistortedPESMode"] = 1;
  json["SCMParFilepath"] = (example_directory / "params.dat").string();
  json["PosFilepath"] = example_slab;
  json["PosFileSlabNum"] = 1;
  json["EnableSurfExpan"] = true;
  json["SCMCutoff"] = 16.0;
  return json;
}

// The example slab file with `edit` applied to each line (1-based numbers), written to `directory`.
std::string edited_example_slab(const test_support::TemporaryDirectory& directory,
                                const std::function<std::string(int, const std::string&)>& edit) {
  std::ifstream in(example_slab);
  std::string text, line;
  for (int number = 1; std::getline(in, line); ++number) text += edit(number, line) + "\n";
  return directory.write("slab.dat", text).string();
}

}  // namespace

// The SCM input with every key, and with only the required ones (the Fortran's defaults: mode 1,
// configuration 1, rescaling on, cut-off 16 bohr). With rescaling the lattice constants are the
// slab's a_T, without them the CRP's.
TEST(Input, ScmKeysValid) {
  Config c = build(scm_json());
  EXPECT_EQ(c.input.scm.parameter_file, (example_directory / "params.dat").string());
  EXPECT_EQ(c.input.scm.position_file, example_slab);
  EXPECT_EQ(c.derived.system.point_group, "C3v");
  EXPECT_EQ(c.derived.system.lattice_constant_x, 4.8896663475);  // a_T, line 3 of the slab file

  Json json = scm_json();
  for (const std::string key : {"DistortedPESMode", "PosFileSlabNum", "EnableSurfExpan", "SCMCutoff"}) json.erase(key);
  c = build(json);
  EXPECT_EQ(c.input.scm.mode, 1);
  EXPECT_EQ(c.input.scm.configuration, 1);
  EXPECT_TRUE(c.input.scm.rescale);
  EXPECT_EQ(c.input.scm.cutoff, 16.0);

  json["EnableSurfExpan"] = false;
  EXPECT_NEAR(build(json).derived.system.lattice_constant_x, 2.601516559 / 0.529177249, 1e-14);  // a₀ of SRP48
}

// SCM D1: SystemID is case-insensitive, as for every PES; the SCM keys are an error with any other
// SystemID (the Fortran ignores them there).
TEST(Input, ScmD1KeysOnlyWithCrpDistrtSurface) {
  Json json = scm_json();
  json["SystemID"] = "crpdistrtsurface";
  EXPECT_EQ(build(json).derived.system.point_group, "C3v");

  for (const std::string& key : scm_keys) {
    SCOPED_TRACE(key);
    json = scm_json();
    json["SystemID"] = "CRPlibrary";
    for (const std::string& other : scm_keys)
      if (other != key) json.erase(other);
    expect_input_error(json, key, "is only used with SystemID CRPDistrtSurface");
    json.erase(key);
    EXPECT_EQ(build(json).derived.system.point_group, "C3v");  // the plain CRP input
  }
  json = reference_json();  // hard wall
  json["SCMCutoff"] = 16.0;
  expect_input_error(json, "SCMCutoff", "is only used with SystemID CRPDistrtSurface");
}

// CRPlibDataSet is required with CRPDistrtSurface: the CRP potential the correction is added to.
TEST(Input, ScmCrpDataSetRequired) {
  Json json = scm_json();
  json.erase("CRPlibDataSet");
  expect_input_error(json, "CRPlibDataSet", "SystemID CRPDistrtSurface needs CRPlibDataSet");
  json["CRPlibDataSet"] = "";
  expect_input_error(json, "CRPlibDataSet", "SystemID CRPDistrtSurface needs CRPlibDataSet");
  json["CRPlibDataSet"] = "no/such/data/set";
  expect_input_error(json, "CRPlibDataSet", "does not exist (no file CRPData)");
}

// SCM D7: only DistortedPESMode 1 (SCM) exists, as in the Fortran.
TEST(Input, ScmD7OnlyMode1) {
  for (const int mode : {0, 2, -1}) {
    Json json = scm_json();
    json["DistortedPESMode"] = mode;
    expect_input_error(json, "DistortedPESMode", "only mode 1 (SCM) exists");
  }
  Json json = scm_json();
  json["DistortedPESMode"] = 1;
  EXPECT_EQ(build(json).input.scm.mode, 1);
}

// SCM D12: both paths are required (the Fortran defaults to ./params.dat and ./SlabPos.dat); a path
// that does not name a readable file is an error naming its key.
TEST(Input, ScmD12PathsRequired) {
  for (const std::string key : {"SCMParFilepath", "PosFilepath"}) {
    SCOPED_TRACE(key);
    Json json = scm_json();
    json.erase(key);
    expect_input_error(json, key, "SystemID CRPDistrtSurface needs " + key);
    json[key] = "";
    expect_input_error(json, key, "SystemID CRPDistrtSurface needs " + key);
    json[key] = "no/such/file.dat";
    expect_input_error(json, key, "no such file");
  }
  EXPECT_NO_THROW((void)build(scm_json()));
}

// SCM D2: UseSymmetry is an error with CRPDistrtSurface (the distorted slab has no symmetry; the
// Fortran switches it off for C3v data sets and would wrongly apply it for C6v ones).
TEST(Input, ScmD2UseSymmetryRejected) {
  Json json = scm_json();
  json["UseSymmetry"] = true;
  expect_input_error(json, "UseSymmetry", "UseSymmetry = true is not supported with SystemID CRPDistrtSurface");
  json["UseSymmetry"] = false;
  EXPECT_FALSE(build(json).input.use_symmetry);
}

// SCM D13: HomoNuc false (HD) is an error with CRPDistrtSurface; the Fortran places the atoms
// symmetrically about the grid point.
TEST(Input, ScmD13HomoNucFalseRejected) {
  Json json = scm_json();
  json["HomoNuc"] = false;
  json["Mass2"] = 2.0141;
  expect_input_error(json, "HomoNuc", "HomoNuc = false is not supported with SystemID CRPDistrtSurface");
  json["HomoNuc"] = true;
  json["Mass2"] = json["Mass1"];
  EXPECT_TRUE(build(json).input.homonuclear);
}

// SCM D14: the slab is FCC111, so the CRP data set must be hexagonal (C3v or C6v); a C4v one is an
// error naming CRPlibDataSet (the message tells it from a data set that fails to load). Valid: SRP48.
TEST(Input, ScmD14HexagonalDataSet) {
  Json json = scm_json();
  json["CRPlibDataSet"] =
      (std::filesystem::path(SPODVR_REFERENCE_DIR).parent_path() / "crplib/reference/c4v_pairwise/dataset").string();
  const std::string message =
      expect_input_error(json, "CRPlibDataSet", "the SCM slab (PosFilepath) is an FCC111 surface");
  EXPECT_NE(message.find("C3v or C6v"), std::string::npos) << message;
}

// The slab number must be 1 … NumSlabConfigs (1 in the example; the Fortran reads configuration 1
// for 0, SCM D6) and the cut-off positive (as in the Fortran).
TEST(Input, ScmSlabNumberAndCutoff) {
  for (const int configuration : {0, -1, 2}) {
    Json json = scm_json();
    json["PosFileSlabNum"] = configuration;
    expect_input_error(json, "PosFileSlabNum",
                       configuration < 1 ? "the configuration number must be at least 1"
                                         : "asked for, but the file has NumSlabConfigs");
  }
  for (const double cutoff : {0.0, -16.0}) {
    Json json = scm_json();
    json["SCMCutoff"] = cutoff;
    expect_input_error(json, "SCMCutoff", "the cut-off must be positive");
  }
  Json json = scm_json();
  json["SCMCutoff"] = 8.0;
  EXPECT_EQ(build(json).input.scm.cutoff, 8.0);
}

// SCM D6: with EnableSurfExpan the slab's lattice constant a_T (line 3) must be positive: it
// becomes the grid's lattice constant and s = a_T/a₀. Without rescaling it is not used.
TEST(Input, ScmD6SlabLatticeConstantPositiveWhenRescaling) {
  const test_support::TemporaryDirectory directory("scm_input");
  for (const std::string lattice_constant : {"0.0", "-4.8896663475"}) {
    SCOPED_TRACE(lattice_constant);
    Json json = scm_json();
    json["PosFilepath"] = edited_example_slab(
        directory, [&](int number, const std::string& line) { return number == 3 ? lattice_constant : line; });
    expect_input_error(json, "PosFilepath", "the lattice constant a_T on line 3 is");
    json["EnableSurfExpan"] = false;
    EXPECT_FALSE(build(json).input.scm.rescale);
  }
}

// SCM D9: the slab's ideal atoms must lie on the CRP lattice. Moving every ideal position by
// 0.01 bohr in x is an error naming PosFilepath; the example slab (within 1.2e-6 bohr) is valid.
TEST(Input, ScmD9IdealAtomsOnTheCrpLattice) {
  const test_support::TemporaryDirectory directory("scm_input");
  Json json = scm_json();
  json["PosFilepath"] = edited_example_slab(directory, [](int number, const std::string& line) {
    if (number <= 3) return line;
    std::istringstream fields(line);
    std::string xd, yd, zd, layer;
    double xi = 0.0, yi = 0.0, zi = 0.0;
    fields >> xd >> yd >> zd >> xi >> yi >> zi >> layer;
    return std::format("{} {} {} {:.17g} {:.17g} {:.17g} {}", xd, yd, zd, xi + 0.01, yi, zi, layer);
  });
  expect_input_error(json, "PosFilepath", "from its site of the CRP lattice");
  EXPECT_NO_THROW((void)build(scm_json()));
}

// UseSymmetry, nPhisOptimal and SymmetryTolerance (eV, default 1e-5, must be > 0). With both
// UseSymmetry and nPhisOptimal nphi becomes a multiple of the point group's multiple.
TEST(Input, SurfaceSymmetryKeys) {
  Json json = reference_json();  // hard wall: C2v
  Config c = build(json);
  EXPECT_FALSE(c.input.use_symmetry);
  EXPECT_FALSE(c.input.phis_optimal);
  EXPECT_EQ(c.input.symmetry_tolerance, 1e-5);
  EXPECT_EQ(c.derived.phi_points_multiple, 1);

  json["UseSymmetry"] = true;
  EXPECT_EQ(build(json).derived.phi_points_multiple, 1);  // only with nPhisOptimal
  json["nPhisOptimal"] = true;
  c = build(json);
  EXPECT_TRUE(c.input.use_symmetry);
  EXPECT_TRUE(c.input.phis_optimal);
  EXPECT_EQ(c.derived.phi_points_multiple, 2);
  json["UseSymmetry"] = false;
  EXPECT_EQ(build(json).derived.phi_points_multiple, 1);  // nPhisOptimal alone has no effect

  json["SymmetryTolerance"] = 1e-3;
  EXPECT_EQ(build(json).input.symmetry_tolerance, 1e-3);
  json["SymmetryTolerance"] = 0.0;
  expect_input_error(json, "SymmetryTolerance", "SymmetryTolerance <= 0.");
  json["SymmetryTolerance"] = -1e-5;
  expect_input_error(json, "SymmetryTolerance", "SymmetryTolerance <= 0.");
}

TEST(Input, PhiMultipleOfTheHexagonalCase) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/cu111_dailight_short_sym/input.json");
  const Config c = build(Json::parse(file));
  EXPECT_EQ(c.derived.phi_points_multiple, 6);  // C6v
}

// Run mode precedence as in Fortran: Potential (6) > onlyAnalysis (2) > NewCalc (0 / 1).
TEST(Input, AnalysisOnlyRunMode) {
  Json json = reference_json();
  json["onlyAnalysis"] = true;
  json["RestartDir"] = "earlier_run";
  EXPECT_EQ(build(json).derived.run_mode, 2);
  json["Potential"] = true;
  EXPECT_EQ(build(json).derived.run_mode, 6);
}

TEST(Input, MjMaxRules) {
  Json json = reference_json();
  json.erase("MJMax");  // default -1 -> JMax
  EXPECT_EQ(build(json).input.mj_max, 4);
  json["MJMax"] = 2;
  EXPECT_EQ(build(json).derived.mj_max_basis, 2);
  json["MJMax"] = 0;  // Fortran quirk: 0 means JMax for the basis
  EXPECT_EQ(build(json).input.mj_max, 0);
  EXPECT_EQ(build(json).derived.mj_max_basis, 4);
}

TEST(Input, CapEnergyOrFactor) {
  Json json = reference_json();
  json.erase("ZOptEnergy");
  expect_input_error(json, "ZOptEnergy", "ZOptFactor or ZOptEnergy should be present.");
  json["ZOptFactor"] = 0.01;
  EXPECT_EQ(build(json).derived.z_opt_factor, 0.01);
  json["ZOptEnergy"] = 0.3;
  expect_input_error(json, "ZOptEnergy", "give either ZOptEnergy or ZOptFactor, not both");  // both given
  json.erase("ZOptFactor");
  json["ZOptEnergy"] = -0.1;
  expect_input_error(json, "ZOptEnergy", "ZOptEnergy <= 0.");
}

TEST(Input, SwitchingReplacesZAsymptotic) {
  Json json = reference_json();
  json["UseSwitching"] = true;
  json["ZminSwitch"] = 5.0;
  json["ZmaxSwitch"] = 6.0;
  const Config config = build(json);
  EXPECT_EQ(config.derived.z_asymptotic, 6.0);
  EXPECT_EQ(config.derived.z_asymptotic_unswitched, 7.0);
  json["ZmaxSwitch"] = 5.0;
  expect_input_error(json, "ZminSwitch", "ZminSwitch >= ZmaxSwitch.");
}

// flux Q1: a flux energy range that ends at or below EMaxWP is valid (the Fortran check rejected
// every such range up to r1713; fixed in r1714); the wave packet spans 0.2..0.6 eV.
TEST(Input, FluxQ1EnergyRangeInsideTheWavePacketIsValid) {
  Json json = reference_json();
  json["doFluxAnalysis"] = true;
  json["rFluxAnalysis"] = 3.05;
  json["EMinFluxAna"] = 0.2;
  json["DeltaEFluxAna"] = 0.01;
  for (const double e_max : {0.6, 0.5}) {  // EMaxFluxAna = EMaxWP, < EMaxWP
    json["EMaxFluxAna"] = e_max;
    const Config config = build(json);
    EXPECT_NO_THROW((void)make_spatial_grids(config)) << "EMaxFluxAna = " << e_max;
    EXPECT_EQ(config.input.e_max_flux_analysis, e_max);
  }
}

// NrOfEs < 2 is an error (as the Fortran since r1719; analysis D1 is resolved); two energies are the
// smallest valid grid.
TEST(Input, NrOfEsTwoIsValid) {
  Json json = reference_json();
  json["NrOfEs"] = 2;
  const Config config = build(json);
  EXPECT_EQ(config.derived.analysis_energy_start, 0.2);
  EXPECT_EQ(config.derived.analysis_energy_step, 0.6 - 0.2);

  json["NrOfEs"] = 1;
  expect_input_error(json, "NrOfEs", "NrOfEs = 1: the analysis needs at least 2 energies");
}

// The valid side of the exact absorber-energy checks (r1718): just above 0 is accepted, and a low-Z
// absorber on the scattering grid may have the energy 0 (its tolerant .RLT. check).
TEST(Input, AbsorberEnergiesJustAboveZeroAreValid) {
  for (const std::string key : {"ZOptEnergy", "ROptEnergy", "ZSpOptEnergy"}) {
    Json json = reference_json();
    json[key] = 1.0e-13;
    EXPECT_NO_THROW((void)make_spatial_grids(build(json))) << key;
  }
  for (const auto& [prefix, energy] :
       std::vector<std::pair<std::string, double>>{{"ZSpOptLow", 1.0e-13}, {"ZOptLow", 0.0}, {"ZOptLow", -5.0e-13}}) {
    Json json = reference_json();
    json[prefix + "Start"] = -2.0;
    json[prefix + "End"] = 8.0;
    json[prefix + "Energy"] = energy;
    EXPECT_NO_THROW((void)make_spatial_grids(build(json))) << prefix << " " << energy;
  }
}

// The valid side of the exact flux spacing check (r1719): a spacing just above 0 is accepted.
TEST(Input, FluxDeltaEJustAboveZeroIsValid) {
  Json json = reference_json();
  json["doFluxAnalysis"] = true;
  json["rFluxAnalysis"] = 3.05;
  json["EMinFluxAna"] = 0.2;
  json["EMaxFluxAna"] = 0.6;
  json["DeltaEFluxAna"] = 5.0e-13;
  EXPECT_NO_THROW((void)make_spatial_grids(build(json)));
}

// The valid side of the analysis limits of r1726: JIn < JMax, VIn < VAna, JIn < JAna (JAna is the
// highest analysed j, of JIn's parity here).
TEST(Input, AnalysisLimitsJustInsideAreValid) {
  Json json = reference_json();
  json["JIn"] = 2;
  EXPECT_NO_THROW((void)make_spatial_grids(build(json)));  // JMax 4
  json = reference_json();
  json["VIn"] = 0;
  json["VAna"] = 1;
  EXPECT_NO_THROW((void)make_spatial_grids(build(json)));
  json = reference_json();
  json["JAna"] = 2;
  EXPECT_NO_THROW((void)make_spatial_grids(build(json)));  // JIn 0
}

// RStart <= 0 is rejected with the 1e-12 absolute tolerance of .RLE. (5.24 InputNew.f90
// ValidateInput), so a start 2e-12 above 0 is valid; 0 and 1e-13 are in the violation table.
TEST(Input, RStartJustAboveZeroIsValid) {
  Json json = reference_json();
  json["RStart"] = 2.0e-12;
  json["ROptEnd"] = 8.25;  // the last r point is now 2e-12 + 55 * 0.15
  const Config config = build(json);
  EXPECT_NO_THROW((void)make_spatial_grids(config));
  EXPECT_EQ(config.input.r_start, 2.0e-12);
}

// Low-Z absorbers on both grids from exactly ZStart to exactly Zana are valid: the Fortran checks
// Start < ZStart and End > Zana with plain comparisons (1e-13 beyond them is in the violation table).
TEST(Input, LowZAbsorbersAtTheirBoundsAreValid) {
  Json json = reference_json();
  json["ZOptLowStart"] = -2.0;
  json["ZOptLowEnd"] = 8.0;
  json["ZOptLowEnergy"] = 0.3;
  json["ZSpOptLowStart"] = -2.0;
  json["ZSpOptLowEnd"] = 8.0;
  json["ZSpOptLowFactor"] = 0.01;
  const Config config = build(json);
  EXPECT_NO_THROW((void)make_spatial_grids(config));
  EXPECT_TRUE(config.input.z_opt_low.enabled());
  EXPECT_TRUE(config.input.zsp_opt_low.enabled());
}

// Strength of a low-Z absorber: from the energy over End - Start (as the high-end absorbers),
// else the factor; an energy of 0 gives 0, and so does one within the 1e-12 tolerance below 0
// (accepted by the Fortran's .RLT. check). Without low-Z keys both absorbers are off.
TEST(Input, LowZAbsorberStrengths) {
  const Config off = build(reference_json());
  EXPECT_FALSE(off.input.z_opt_low.enabled());
  EXPECT_FALSE(off.input.zsp_opt_low.enabled());
  EXPECT_EQ(off.derived.z_opt_low_factor, 0.0);
  EXPECT_EQ(off.derived.zsp_opt_low_factor, 0.0);

  Json json = reference_json();
  json["ZOptLowStart"] = -1.0;
  json["ZOptLowEnd"] = 4.0;
  json["ZOptLowEnergy"] = 0.3;
  json["ZSpOptLowStart"] = -2.0;
  json["ZSpOptLowEnd"] = 6.0;
  json["ZSpOptLowFactor"] = 0.0125;
  Config config = build(json);
  const double ekin = units::legacy.ev_to_hartree(0.3);
  EXPECT_EQ(config.derived.z_opt_low_factor, optical_potential_strength(ekin, 4.0 - (-1.0), config.derived.total_mass));
  EXPECT_EQ(config.derived.zsp_opt_low_factor, 0.0125);

  for (const double energy : {0.0, -1.0e-13}) {
    json["ZOptLowEnergy"] = energy;
    config = build(json);
    EXPECT_EQ(config.derived.z_opt_low_factor, 0.0) << energy;
  }
}

// writePotential with sites on the grid is accepted in a new run and a potential-only run; a site off the grid is
// refused naming the nearest grid point.
TEST(Input, WritePotentialWithSitesOnTheGrid) {
  Json json = reference_json();
  json["writePotential"] = true;
  json["PlotSites"] = site_list({{"origin", 0.0, 0.0}, {"far", 0.75, 0.5}});
  const Config config = build(json);
  ASSERT_EQ(config.input.plot_sites.size(), 2u);
  EXPECT_EQ(plot_site_ixy(config.input.plot_sites[1], config.input.nx, config.input.ny), 2 + config.input.ny * 3);
  json["Potential"] = true;
  EXPECT_TRUE(build(json).input.write_potential);
  json["PlotSites"] = site_list({{"between", 0.3, 0.0}});
  try {
    build(json);
    FAIL() << "no InputError";
  } catch (const InputError& e) {
    EXPECT_NE(std::string(e.what()).find("the nearest point is x = 0.25 (index 1)"), std::string::npos) << e.what();
  }
}

// A site within 1e-9 of a grid point (in units of the grid step) is on it, one just outside is not; a site that
// rounds to 1 is refused (it is the next cell's 0), as one below 0 is.
TEST(Input, PlotSitesOnTheGridWithinTheTolerance) {
  Json json = reference_json();  // nX = 4: x * 4 within 1e-9 of a whole number
  json["writePotential"] = true;
  json["PlotSites"] = site_list({{"inside", 0.25 + 2e-10, 0.0}});
  EXPECT_NO_THROW(build(json));
  for (const double x : {0.25 + 5e-10, 1.0 - 1e-12, -1e-12}) {
    json["PlotSites"] = site_list({{"edge", x, 0.0}});
    EXPECT_THROW(build(json), InputError) << x;
  }
}

// writeDensities in a new run (and in a restart, the other propagating mode) is accepted.
TEST(Input, WriteDensitiesInAPropagatingRun) {
  Json json = reference_json();
  json["writeDensities"] = true;
  EXPECT_TRUE(build(json).input.write_densities);
  const double duration = json["tEnd"].get<double>() - json["tStart"].get<double>();
  json["NewCalc"] = false;
  json["RestartDir"] = "earlier";
  json["tStart"] = 1000.0;
  json["tEnd"] = 1000.0 + duration;
  EXPECT_TRUE(build(json).input.write_densities);
}

// FrameExtras with dtFrames: every known name is accepted, in any order; none by default.
TEST(Input, FrameExtrasWithFrames) {
  Json json = reference_json();
  EXPECT_TRUE(build(json).input.frame_extras.empty());
  json["dtFrames"] = 100.0;
  json["FrameExtras"] = {"momentum", "angular", "diffraction", "rotation_moments", "vibration"};
  EXPECT_EQ(build(json).input.frame_extras.size(), 5u);
}

// NormBudget in a new run and in a restart (the propagating modes) is accepted; it is off by default.
TEST(Input, NormBudgetInAPropagatingRun) {
  Json json = reference_json();
  EXPECT_FALSE(build(json).input.norm_budget);
  json["NormBudget"] = true;
  EXPECT_TRUE(build(json).input.norm_budget);
  const double duration = json["tEnd"].get<double>() - json["tStart"].get<double>();
  json["NewCalc"] = false;
  json["RestartDir"] = "earlier";
  json["tStart"] = 1000.0;
  json["tEnd"] = 1000.0 + duration;
  EXPECT_TRUE(build(json).input.norm_budget);
}

// EnergyResolvedSites with EnergyResolvedEnergies: sites of PlotSites (which need no other plot key then) and collision
// energies at both ends of the packet's range (EMinWP 0.2 .. EMaxWP 0.6 eV), in a new run and a restart; none by
// default.
TEST(Input, EnergyResolvedSumsInAPropagatingRun) {
  Json json = reference_json();
  EXPECT_TRUE(build(json).input.energy_resolved_sites.empty());
  json["PlotSites"] = site_list({{"a", 0.0, 0.0}, {"b", 0.5, 0.0}});
  json["EnergyResolvedSites"] = {"b", "a"};
  json["EnergyResolvedEnergies"] = {0.2, 0.45, 0.6};
  const Config config = build(json);
  EXPECT_EQ(config.input.energy_resolved_sites, (std::vector<std::string>{"b", "a"}));
  EXPECT_EQ(config.input.energy_resolved_energies, (std::vector<double>{0.2, 0.45, 0.6}));
  const double duration = json["tEnd"].get<double>() - json["tStart"].get<double>();
  json["NewCalc"] = false;
  json["RestartDir"] = "earlier";
  json["tStart"] = 1000.0;
  json["tEnd"] = 1000.0 + duration;
  EXPECT_EQ(build(json).input.energy_resolved_energies.size(), 3u);
}

// One violation per check of ValidateInput / ReadInput / Setup.f90: the input change, the parameter the
// error names, and the fixed text of the check that must fire.
struct Violation {
  std::string name;
  std::function<void(Json&)> modify;
  std::string parameter;
  std::string message;
};

void PrintTo(const Violation& v, std::ostream* os) { *os << v.name; }

class InputViolation : public ::testing::TestWithParam<Violation> {};

TEST_P(InputViolation, IsRejected) {
  Json json = reference_json();
  GetParam().modify(json);
  expect_input_error(json, GetParam().parameter, GetParam().message);
}

INSTANTIATE_TEST_SUITE_P(
    Checks, InputViolation,
    ::testing::Values(
        Violation{"SystemID", [](Json& j) { j["SystemID"] = "H2Pt111"; }, "SystemID", "unknown system ID"},
        Violation{"ZDoF", [](Json& j) { j["ZDoF"] = false; }, "ZDoF", "Currently only 6D calculations are supported."},
        Violation{"PhiDoF", [](Json& j) { j["PhiDoF"] = false; }, "PhiDoF", "Currently only 6D calculations are supported."},
        Violation{"Mass1", [](Json& j) { j["Mass1"] = 0.0; j["Mass2"] = 0.0; }, "Mass1", "Mass of atom 1 <= 0."},
        Violation{"HomoMasses", [](Json& j) { j["Mass2"] = 2.0; }, "HomoNuc", "masses are different, but molecule is homo-nuclear."},
        Violation{"HeteroMasses", [](Json& j) { j["HomoNuc"] = false; }, "HomoNuc", "masses are the same, but molecule is hetero-nuclear."},
        Violation{"nZFft", [](Json& j) { j["nZ"] = 177; }, "nZ", "nZ is not a valid FFT size."},
        Violation{"nZ", [](Json& j) { j["nZ"] = 0; }, "nZ", "nZ < 1."},
        Violation{"dZ", [](Json& j) { j["dZ"] = 0.0; }, "dZ", "dZ <= 0."},
        Violation{"nZSpecularFft", [](Json& j) { j["nZSpecular"] = 281; }, "nZSpecular", "nZSpecular is not a valid FFT size."},
        Violation{"nZSpecularSmall", [](Json& j) { j["nZSpecular"] = 160; j["ZSpOptStart"] = 8.0; j["ZSpOptEnd"] = 10.0; j["Z0"] = 7.5; }, "nZSpecular", "nZSpecular < nZ (the specular grid must contain the scattering grid)"},
        Violation{"RStart", [](Json& j) { j["RStart"] = -0.1; }, "RStart", "RStart <= 0."},
        // RStart <= 0 within 1e-12; ROptEnd moves onto the shifted r grid, so only RStart is violated.
        Violation{"RStartZero", [](Json& j) { j["RStart"] = 0.0; j["ROptEnd"] = 8.25; }, "RStart", "RStart <= 0."},
        Violation{"RStartWithinTolerance", [](Json& j) { j["RStart"] = 1.0e-13; j["ROptEnd"] = 8.25; }, "RStart", "RStart <= 0."},
        Violation{"nR", [](Json& j) { j["nR"] = 57; }, "nR", "nR is not a valid FFT size."},
        Violation{"dR", [](Json& j) { j["dR"] = -0.15; }, "dR", "dR <= 0."},
        Violation{"nX", [](Json& j) { j["nX"] = 5; }, "nX", "nX is not a valid FFT size."},
        Violation{"nY", [](Json& j) { j["nY"] = 0; }, "nY", "nY < 1."},
        Violation{"JIn", [](Json& j) { j["JIn"] = -2; }, "JIn", "JIn < 0."},
        Violation{"JMax", [](Json& j) { j["JMax"] = -1; j["MJMax"] = 0; }, "JMax", "JMax < 0."},
        Violation{"MJMax", [](Json& j) { j["MJMax"] = 6; }, "MJMax", "|MJMax| > JMax."},
        Violation{"MJIn", [](Json& j) { j["MJIn"] = 1; }, "MJIn", "|MJIn| > JIn."},
        Violation{"JInGtJMax", [](Json& j) { j["JIn"] = 6; }, "JIn", "JIn >= JMax."},
        Violation{"JInEqJMax", [](Json& j) { j["JIn"] = 4; }, "JIn", "JIn >= JMax."},
        Violation{"MJInGtMJMax", [](Json& j) { j["JIn"] = 2; j["MJIn"] = 2; j["MJMax"] = 0; }, "MJIn", "|MJIn| > |MJMax|."},
        Violation{"ParityJ", [](Json& j) { j["JMax"] = 5; j["MJMax"] = 5; }, "JMax", "the molecule is homo-nuclear but JMax is even and JIn is odd, or vice versa."},
        Violation{"ParityMJ", [](Json& j) { j["MJMax"] = 3; }, "MJMax", "the molecule is homo-nuclear, but JMax is even and |MJMax| is odd, or vice versa."},
        Violation{"ProjOp0", [](Json& j) { j["ProjOp"] = 0; }, "ProjOp", "ProjOp = 0 is not supported anymore."},
        Violation{"ProjOp1", [](Json& j) { j["ProjOp"] = 1; }, "ProjOp", "ProjOp = 1 is a legacy option and not supported; use ProjOp = 2."},
        Violation{"ProjOp3", [](Json& j) { j["ProjOp"] = 3; }, "ProjOp", "ProjOp must be 0, 1 or 2."},
        Violation{"tStartNew", [](Json& j) { j["tStart"] = 10.0; j["tEnd"] = 14010.0; }, "tStart", "This is a new calculation but tStart /= 0."},
        Violation{"tStartRestart", [](Json& j) { j["NewCalc"] = false; }, "tStart", "This is a restart but tStart = 0."},
        Violation{"OnlyAnalysisNeedsRestartDir", [](Json& j) { j["onlyAnalysis"] = true; }, "RestartDir", "An analysis-only run (onlyAnalysis) needs RestartDir: the directory with the BK coefficient file."},
        // Plot data (validate_outputs.cpp): only a propagating run writes them.
        Violation{"WriteDensitiesPotentialOnly", [](Json& j) { j["Potential"] = true; j["writeDensities"] = true; }, "writeDensities", "only a propagating run (a new run or a restart) writes densities"},
        Violation{"WriteDensitiesAnalysisOnly",
                  [](Json& j) { j["onlyAnalysis"] = true; j["RestartDir"] = "earlier"; j["writeDensities"] = true; }, "writeDensities", "only a propagating run (a new run or a restart) writes densities"},
        Violation{"WritePotentialAnalysisOnly",
                  [](Json& j) { j["onlyAnalysis"] = true; j["RestartDir"] = "earlier"; j["writePotential"] = true; }, "writePotential", "an analysis-only run (onlyAnalysis) computes no 6D potential to write"},
        // PlotSites: the list and its fields (json_input.cpp), then the values (validate_outputs.cpp).
        Violation{"PlotSitesNotAList", [](Json& j) { j["writePotential"] = true; j["PlotSites"] = 3; }, "PlotSites", "expected a list of sites"},
        Violation{"PlotSitesNotObjects", [](Json& j) { j["writePotential"] = true; j["PlotSites"] = Json::array({1}); }, "PlotSites", "each site must be an object with a name, x and y"},
        Violation{"PlotSitesUnknownField",
                  [](Json& j) { j["writePotential"] = true; j["PlotSites"] = site_list({{"a", 0.0, 0.0}}); j["PlotSites"][0]["z"] = 1.0; }, "PlotSites", "unknown site field 'z'"},
        Violation{"PlotSitesMissingField",
                  [](Json& j) { j["writePotential"] = true; j["PlotSites"] = site_list({{"a", 0.0, 0.0}}); j["PlotSites"][0].erase("y"); }, "PlotSites", "a site without 'y'"},
        Violation{"PlotSitesWrongType",
                  [](Json& j) { j["writePotential"] = true; j["PlotSites"] = site_list({{"a", 0.0, 0.0}}); j["PlotSites"][0]["name"] = 1; }, "PlotSites", "expected a string"},
        Violation{"PlotSitesWithoutUser", [](Json& j) { j["PlotSites"] = site_list({{"a", 0.0, 0.0}}); }, "PlotSites", "no plot key uses the sites: set writePotential, dtFrames or EnergyResolvedSites"},
        // dtFrames (frames/): propagating runs only, on whole multiples of dtAna (itself of dtStep), from a tStart
        // and to a tEnd on the frame grid.
        Violation{"DtFramesNotPositive", [](Json& j) { j["dtFrames"] = 0.0; }, "dtFrames", "dtFrames <= 0."},
        Violation{"DtFramesBelowDtAna", [](Json& j) { j["dtFrames"] = 5.0; }, "dtFrames", "dtFrames < dtAna."},
        Violation{"DtFramesNotAMultipleOfDtAna", [](Json& j) { j["dtFrames"] = 15.0; }, "dtFrames", "dtFrames must be a whole multiple of dtAna."},
        Violation{"DtFramesTEndOffTheGrid", [](Json& j) { j["dtFrames"] = 300.0; }, "dtFrames", "tEnd must be a whole multiple of dtFrames."},
        Violation{"DtFramesTStartOffTheGrid",
                  [](Json& j) {
                    j["NewCalc"] = false;
                    j["RestartDir"] = "earlier";
                    j["tStart"] = 150.0;
                    j["tEnd"] = 150.0 + j["nStep"].get<double>() * j["dtStep"].get<double>();
                    j["dtFrames"] = 100.0;  // tStart 150 is not on its grid (checked before tEnd)
                  }, "dtFrames", "tStart must be a whole multiple of dtFrames (the frames of a restart continue those of the run before)"},
        Violation{"DtFramesDtAnaNotAMultipleOfDtStep", [](Json& j) { j["dtAna"] = 5.0; j["dtFrames"] = 100.0; }, "dtAna", "dtAna / dtStep is not a whole number with dtFrames (the program truncates it, so its records would come at other times than dtAna says)"},
        Violation{"DtFramesAnalysisOnly",
                  [](Json& j) { j["onlyAnalysis"] = true; j["RestartDir"] = "earlier"; j["dtFrames"] = 100.0; }, "dtFrames", "only a propagating run (a new run or a restart) writes frames"},
        Violation{"DtFramesPotentialOnly", [](Json& j) { j["Potential"] = true; j["dtFrames"] = 100.0; }, "dtFrames", "only a propagating run (a new run or a restart) writes frames"},
        Violation{"DtFramesWrite2DPES", [](Json& j) { j["write2DPES"] = true; j["dtFrames"] = 100.0; }, "dtFrames", "write2DPES stops the run after the cut of the PES, before the propagation"},
        // FrameExtras (frames/): only with dtFrames, known names, each once, a list of names.
        Violation{"FrameExtrasWithoutFrames", [](Json& j) { j["FrameExtras"] = {"angular"}; }, "FrameExtras", "the extras are datasets of frames/: set dtFrames"},
        Violation{"FrameExtrasUnknownName",
                  [](Json& j) { j["dtFrames"] = 100.0; j["FrameExtras"] = {"angular", "spin"}; }, "FrameExtras", "unknown extra 'spin': known are angular, diffraction, vibration, rotation_moments, momentum"},
        Violation{"FrameExtrasTwice",
                  [](Json& j) { j["dtFrames"] = 100.0; j["FrameExtras"] = {"momentum", "momentum"}; }, "FrameExtras", "'momentum' is listed twice"},
        Violation{"FrameExtrasNotAList", [](Json& j) { j["dtFrames"] = 100.0; j["FrameExtras"] = "angular"; }, "FrameExtras", "expected a list of names"},
        Violation{"FrameExtrasNotNames", [](Json& j) { j["dtFrames"] = 100.0; j["FrameExtras"] = {"angular", 3}; }, "FrameExtras", "expected a string"},
        // NormBudget (norm_budget/): propagating runs only, one record per analysis step of whole time steps.
        Violation{"NormBudgetDtAnaNotAMultipleOfDtStep", [](Json& j) { j["dtAna"] = 5.0; j["NormBudget"] = true; }, "dtAna", "dtAna / dtStep is not a whole number with NormBudget (the program truncates it, so its records would come at other times than dtAna says)"},
        Violation{"NormBudgetAnalysisOnly",
                  [](Json& j) { j["onlyAnalysis"] = true; j["RestartDir"] = "earlier"; j["NormBudget"] = true; }, "NormBudget", "only a propagating run (a new run or a restart) has a norm budget"},
        Violation{"NormBudgetPotentialOnly", [](Json& j) { j["Potential"] = true; j["NormBudget"] = true; }, "NormBudget", "only a propagating run (a new run or a restart) has a norm budget"},
        Violation{"NormBudgetWrite2DPES", [](Json& j) { j["write2DPES"] = true; j["NormBudget"] = true; }, "NormBudget", "write2DPES stops the run after the cut of the PES, before the propagation"},
        // EnergyResolvedSites / EnergyResolvedEnergies (energy_resolved/): both or neither, propagating runs only,
        // sites of PlotSites once each, energies in the packet's range, at most 256 pairs, lists of their type.
        Violation{"EnergyResolvedSitesWithoutEnergies",
                  [](Json& j) { j["PlotSites"] = site_list({{"a", 0.0, 0.0}}); j["EnergyResolvedSites"] = {"a"}; }, "EnergyResolvedSites", "set EnergyResolvedEnergies too (the sums need sites and energies)"},
        Violation{"EnergyResolvedEnergiesWithoutSites", [](Json& j) { j["EnergyResolvedEnergies"] = {0.3}; }, "EnergyResolvedEnergies", "set EnergyResolvedSites too (the sums need sites and energies)"},
        Violation{"EnergyResolvedPotentialOnly",
                  [](Json& j) { energy_resolved(j); j["Potential"] = true; }, "EnergyResolvedSites", "only a propagating run (a new run or a restart) sums the wave packet"},
        Violation{"EnergyResolvedAnalysisOnly",
                  [](Json& j) { energy_resolved(j); j["onlyAnalysis"] = true; j["RestartDir"] = "earlier"; }, "EnergyResolvedSites", "only a propagating run (a new run or a restart) sums the wave packet"},
        Violation{"EnergyResolvedWrite2DPES", [](Json& j) { energy_resolved(j); j["write2DPES"] = true; }, "EnergyResolvedSites", "write2DPES stops the run after the cut of the PES, before the propagation"},
        Violation{"EnergyResolvedUnknownSite", [](Json& j) { energy_resolved(j); j["EnergyResolvedSites"] = {"a", "c"}; }, "EnergyResolvedSites", "'c' is not a site of PlotSites"},
        Violation{"EnergyResolvedSiteTwice", [](Json& j) { energy_resolved(j); j["EnergyResolvedSites"] = {"a", "a"}; }, "EnergyResolvedSites", "'a' is listed twice"},
        Violation{"EnergyResolvedEnergyBelowThePacket",
                  [](Json& j) { energy_resolved(j); j["EnergyResolvedEnergies"] = {0.3, 0.19}; }, "EnergyResolvedEnergies", "lies outside the packet's collision energies"},
        Violation{"EnergyResolvedEnergyAboveThePacket",
                  [](Json& j) { energy_resolved(j); j["EnergyResolvedEnergies"] = {0.61}; }, "EnergyResolvedEnergies", "lies outside the packet's collision energies"},
        Violation{"EnergyResolvedTooManyPairs",
                  [](Json& j) {
                    energy_resolved(j);
                    j["PlotSites"] = Json::array();
                    j["EnergyResolvedSites"] = Json::array();
                    for (int k = 0; k < max_plot_sites; ++k) {
                      j["PlotSites"].push_back(site_list({{"s" + std::to_string(k), 0.0, 0.0}})[0]);
                      j["EnergyResolvedSites"].push_back("s" + std::to_string(k));
                    }
                    j["EnergyResolvedEnergies"] = Json::array();
                    for (int k = 0; k <= 16; ++k) j["EnergyResolvedEnergies"].push_back(0.2 + 0.02 * k);  // 16 x 17
                  }, "EnergyResolvedSites", "(site, energy) pairs: at most 256"},
        Violation{"EnergyResolvedEnergiesNotAList", [](Json& j) { energy_resolved(j); j["EnergyResolvedEnergies"] = 0.3; }, "EnergyResolvedEnergies", "expected a list of numbers"},
        Violation{"EnergyResolvedEnergiesNotNumbers",
                  [](Json& j) { energy_resolved(j); j["EnergyResolvedEnergies"] = {0.3, "0.4"}; }, "EnergyResolvedEnergies", "expected a real number"},
        Violation{"PlotSitesOutsideTheCell",
                  [](Json& j) { j["writePotential"] = true; j["PlotSites"] = site_list({{"a", 1.0, 0.0}}); }, "PlotSites", "is not in [0, 1) (a fraction of the lattice vector)"},
        Violation{"PlotSitesOffTheGrid",
                  [](Json& j) { j["writePotential"] = true; j["PlotSites"] = site_list({{"a", 0.0, 0.3}}); }, "PlotSites", "is not on the grid"},
        Violation{"PlotSitesSameName",
                  [](Json& j) { j["writePotential"] = true; j["PlotSites"] = site_list({{"a", 0.0, 0.0}, {"a", 0.25, 0.0}}); }, "PlotSites", "two sites are named 'a'"},
        Violation{"PlotSitesWithoutName",
                  [](Json& j) { j["writePotential"] = true; j["PlotSites"] = site_list({{"", 0.0, 0.0}}); }, "PlotSites", "a site without a name"},
        Violation{"PlotSitesTooMany",
                  [](Json& j) {
                    j["writePotential"] = true;
                    j["PlotSites"] = Json::array();
                    for (int k = 0; k <= max_plot_sites; ++k) j["PlotSites"].push_back(site_list({{"s" + std::to_string(k), 0.0, 0.0}})[0]);
                  }, "PlotSites", "sites: at most 16"},
        Violation{"tEnd", [](Json& j) { j["tEnd"] = 14002.0; }, "tEnd", "tEnd /= tStart + nStep*dtStep."},
        Violation{"dtStep", [](Json& j) { j["dtStep"] = 0.0; j["nStep"] = 0; j["tEnd"] = 0.0; }, "dtStep", "dtStep <= 0."},
        // tEnd < 0 comes first; nStep < 0 itself is reached only inside the 1e-12 band (nStepWithinTolerance).
        Violation{"nStep", [](Json& j) { j["nStep"] = -1; j["tEnd"] = -2.0; }, "tEnd", "tEnd < 0."},
        Violation{"dtAna", [](Json& j) { j["dtAna"] = 1.0; }, "dtAna", "dtAna < dtStep."},
        Violation{"dtAnaInteger", [](Json& j) { j["dtAna"] = 30.0; j["dtAve"] = 90.0; j["dtStep"] = 2.0; j["nStep"] = 7001; j["tEnd"] = 14002.0; }, "dtAna", "Propagation time must be an integer times dtAna."},
        Violation{"dtAve", [](Json& j) { j["dtAve"] = 1.0; }, "dtAve", "dtAve < dtStep."},
        Violation{"dtAveLtAna", [](Json& j) { j["dtAve"] = 4.0; }, "dtAve", "dtAve < dtAna."},
        Violation{"dtAveInteger", [](Json& j) { j["dtAve"] = 25.0; }, "dtAve", "dtAve must be an integer times dtAna."},
        Violation{"JAnaNeg", [](Json& j) { j["JAna"] = -1; }, "JAna", "JAna < 0."},
        Violation{"VAnaNeg", [](Json& j) { j["VAna"] = -1; }, "VAna", "VAna < 0."},
        Violation{"DifOrdNeg", [](Json& j) { j["DifOrdAna"] = -1; }, "DifOrdAna", "DifOrdAna < 0."},
        Violation{"NrOfEs", [](Json& j) { j["NrOfEs"] = -1; }, "NrOfEs", "NrOfEs = -1: the analysis needs at least 2 energies"},
        Violation{"NrOfEsNen95", [](Json& j) { j["NrOfEs"] = -95; }, "NrOfEs", "NrOfEs = -95: the analysis needs at least 2 energies"},
        Violation{"AnalysisD1NrOfEsZero", [](Json& j) { j["NrOfEs"] = 0; }, "NrOfEs", "NrOfEs = 0: the analysis needs at least 2 energies, (EMaxWP-EMinWP)/(NrOfEs-1) apart; use NrOfEs >= 2."},
        Violation{"AnalysisD1NrOfEsOne", [](Json& j) { j["NrOfEs"] = 1; }, "NrOfEs", "NrOfEs = 1: the analysis needs at least 2 energies, (EMaxWP-EMinWP)/(NrOfEs-1) apart; use NrOfEs >= 2."},
        Violation{"ZanaLtZasymp", [](Json& j) { j["Zana"] = 6.0; }, "Zana", "ZAna < ZAsymptotic."},
        Violation{"ZanaOffGrid", [](Json& j) { j["Zana"] = 8.01; }, "Zana", "ZAna not on a grid point."},
        Violation{"ZanaGtZOptStart", [](Json& j) { j["Zana"] = 9.04; }, "Zana", "ZAnalysis > ZOptStart."},
        Violation{"JAnaGtJMax", [](Json& j) { j["JAna"] = 6; }, "JAna", "JAna > JMax."},
        Violation{"VAnaGtNr", [](Json& j) { j["VAna"] = 56; }, "VAna", "VAna > nR-1."},
        Violation{"DifOrdGt", [](Json& j) { j["DifOrdAna"] = 7; }, "DifOrdAna", "DifOrdAna > nX+nY-2 (the maximum diffraction order possible)."},
        Violation{"dtAnaFine", [](Json& j) { j["EMaxWP"] = 9.0; }, "dtAna", "dtAna is too big for fine analysis."},
        Violation{"FluxMissing", [](Json& j) { j["doFluxAnalysis"] = true; }, "rFluxAnalysis", "rFluxAnalysis must be present"},
        Violation{"FluxOffGrid", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.01; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.01; }, "rFluxAnalysis", "rFluxAnalysis not on a grid point."},
        Violation{"FluxEMin", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.1; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.01; }, "EMinFluxAna", "EMinFluxAna < EMinWP"},
        // The range checks of r1714, and their 1e-12 bands, which pass on to the next check.
        Violation{"FluxEMinAboveEMaxWP", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.7; j["EMaxFluxAna"] = 0.8; j["DeltaEFluxAna"] = 0.01; }, "EMinFluxAna", "EMinFluxAna > EMaxWP"},
        Violation{"FluxEMaxBelowEMinWP", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.1; j["DeltaEFluxAna"] = 0.01; }, "EMaxFluxAna", "EMaxFluxAna < EMinWP"},
        Violation{"FluxEMinWithinToleranceAboveEMaxWPReachesTheNextCheck", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.6 + 5.0e-13; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.01; }, "EMinFluxAna", "EMinFluxAna >= EMaxFluxAna"},
        Violation{"FluxEMaxWithinToleranceBelowEMinWPReachesTheNextCheck", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.2 - 5.0e-13; j["DeltaEFluxAna"] = 0.01; }, "EMinFluxAna", "EMinFluxAna >= EMaxFluxAna"},
        // flux Q1: EMaxFluxAna above EMaxWP (0.6 eV) is rejected.
        Violation{"FluxQ1EMax", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.7; j["DeltaEFluxAna"] = 0.01; }, "EMaxFluxAna", "EMaxFluxAna > EMaxWP"},
        // flux Q2: a zero or negative energy spacing is rejected (the Fortran accepted 0 up to r1713).
        Violation{"FluxQ2DeltaZero", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.0; }, "DeltaEFluxAna", "DeltaEFluxAna <= 0"},
        Violation{"FluxQ2DeltaNegative", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = -0.01; }, "DeltaEFluxAna", "DeltaEFluxAna <= 0"},
        Violation{"FluxDelta", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.5; }, "DeltaEFluxAna", "DeltaEFluxAna > EMaxFluxAna-EMinFluxAna"},
        Violation{"FluxROpt", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 5.0; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.01; }, "rFluxAnalysis", "rFluxAnalysis > ROptStart"},
        Violation{"EMinWP", [](Json& j) { j["EMinWP"] = 0.0; }, "EMinWP", "EMinWP <= 0."},
        Violation{"EMinGeEMax", [](Json& j) { j["EMinWP"] = 0.6; }, "EMinWP", "EMinWP >= EMaxWP."},
        Violation{"Z0LtZasymp", [](Json& j) { j["Z0"] = 6.0; }, "Z0", "Z0 < ZAsymptotic."},
        Violation{"Z0GtZSpOpt", [](Json& j) { j["Z0"] = 17.0; }, "Z0", "Z0 > ZSpOptStart."},
        Violation{"VIn", [](Json& j) { j["VIn"] = -1; }, "VIn", "VIn < 0."},
        Violation{"VInGtNr", [](Json& j) { j["VIn"] = 56; }, "VIn", "VIn > nR-1."},
        Violation{"VInGtVAna", [](Json& j) { j["VIn"] = 2; j["VAna"] = 1; }, "VIn", "VIn >= VAna."},
        Violation{"VInEqVAna", [](Json& j) { j["VIn"] = 1; j["VAna"] = 1; }, "VIn", "VIn >= VAna."},
        Violation{"VAnaZeroGiven", [](Json& j) { j["VAna"] = 0; }, "VIn", "VIn >= VAna."},
        Violation{"JInGtJAna", [](Json& j) { j["JIn"] = 2; j["JAna"] = 0; }, "JIn", "JIn >= JAna."},
        Violation{"JInEqJAna", [](Json& j) { j["JAna"] = 0; }, "JIn", "JIn >= JAna."},
        Violation{"VMinGtVMax", [](Json& j) { j["ImposeVMin"] = true; j["ImposeVMax"] = true; j["VMin"] = 2.0; j["VMax"] = 1.0; }, "VMin", "VMin > VMax."},
        Violation{"ZOptOrder", [](Json& j) { j["ZOptEnd"] = 8.5; }, "ZOptStart", "ZOptStart > ZOptEnd."},
        Violation{"ZOptEnd", [](Json& j) { j["ZOptEnd"] = 13.0; }, "ZOptEnd", "ZOptEnd is too big for the scattering grid in Z."},
        Violation{"ROptOrder", [](Json& j) { j["ROptEnd"] = 4.0; }, "ROptStart", "ROptStart > ROptEnd."},
        Violation{"ROptEnd", [](Json& j) { j["ROptEnd"] = 9.2; }, "ROptEnd", "ROptEnd is too big for the grid in r."},
        Violation{"ZSpOptOrder", [](Json& j) { j["ZSpOptEnd"] = 16.0; }, "ZSpOptStart", "ZSpOptStart > ZSpOptEnd."},
        Violation{"ZSpOptEnd", [](Json& j) { j["ZSpOptEnd"] = 21.0; }, "ZSpOptEnd", "ZSpOptEnd is too big for the specular grid in Z."},
        // ZSpOptStart < ZAsymptotic cannot be reached on its own: Z0 >= ZAsymptotic and
        // Z0 <= ZSpOptStart are checked first.
        // Z0 > ZSpOptStart comes first; ZSpOptStart < ZAsymptotic itself is reached only inside the 1e-12 band
        // (ZSpOptStartBelowZAsymptotic).
        Violation{"ZSpOptStart", [](Json& j) { j["ZSpOptStart"] = 6.0; }, "Z0", "Z0 > ZSpOptStart."},
        Violation{"EParallel", [](Json& j) { j["EParallel"] = -0.1; }, "EParallel", "EParallel < 0."},
        Violation{"PhiParallelNeg", [](Json& j) { j["PhiParallel"] = -1.0; }, "PhiParallel", "PhiParallel < 0."},
        Violation{"PhiParallelBig", [](Json& j) { j["PhiParallel"] = 361.0; }, "PhiParallel", "PhiParallel > 360 degrees."},
        Violation{"X2DPESNeg", [](Json& j) { j["write2DPES"] = true; j["x2DPES"] = -0.1; }, "x2DPES", "x2DPES < 0."},
        Violation{"X2DPESBig", [](Json& j) { j["write2DPES"] = true; j["x2DPES"] = 6.3; }, "x2DPES", "x2DPES > alat."},
        Violation{"Y2DPESNeg", [](Json& j) { j["write2DPES"] = true; j["y2DPES"] = -0.1; }, "y2DPES", "y2DPES < 0."},
        Violation{"Y2DPESBig", [](Json& j) { j["write2DPES"] = true; j["y2DPES"] = 5.3; }, "y2DPES", "y2DPES > blat."},
        Violation{"Theta2DPESNeg", [](Json& j) { j["write2DPES"] = true; j["theta2DPES"] = -1.0; }, "theta2DPES", "theta2DPES < 0."},
        Violation{"Theta2DPESBig", [](Json& j) { j["write2DPES"] = true; j["theta2DPES"] = 361.0; }, "theta2DPES", "theta2DPES > 360."},
        Violation{"Phi2DPESNeg", [](Json& j) { j["write2DPES"] = true; j["phi2DPES"] = -1.0; }, "phi2DPES", "phi2DPES < 0."},
        Violation{"Phi2DPESBig", [](Json& j) { j["write2DPES"] = true; j["phi2DPES"] = 361.0; }, "phi2DPES", "phi2DPES > 360."},
        Violation{"ZAnaSpDynamic", [](Json& j) { j["ZAnaSp"] = 12.0; }, "ZAnaSp", "ZAnaSp can only be used with DynamZSpAna set to FALSE"},
        Violation{"ZAnaSpGtOpt", [](Json& j) { j["DynamZSpAna"] = false; j["ZAnaSp"] = 17.0; }, "ZAnaSp", "ZAnaSp was set > ZSpOptStart"},
        Violation{"ZAnaSpLtZ0", [](Json& j) { j["DynamZSpAna"] = false; j["ZAnaSp"] = 11.0; }, "ZAnaSp", "ZAnaSp was set < Z0"},
        // Low-Z absorbers. low-Z D3: an incomplete set of keys (the Fortran ignores it).
        Violation{"LowZD3OnlyStart", [](Json& j) { j["ZOptLowStart"] = -2.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowEnd", "ZOptLowEnd is missing: a low-Z absorber needs ZOptLowStart and ZOptLowEnd; give Start, End and one of Energy/Factor, or none of them"},
        Violation{"LowZD3OnlyEnd", [](Json& j) { j["ZOptLowEnd"] = 8.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowStart", "ZOptLowStart is missing: a low-Z absorber needs ZOptLowStart and ZOptLowEnd; give Start, End and one of Energy/Factor, or none of them"},
        Violation{"LowZD3OnlyStrength", [](Json& j) { j["ZSpOptLowFactor"] = 0.01; }, "ZSpOptLowStart", "ZSpOptLowStart is missing: a low-Z absorber needs ZSpOptLowStart and ZSpOptLowEnd; give Start, End and one of Energy/Factor, or none of them"},
        // low-Z D1 (Start > End, ignored by the Fortran) and low-Z D2 (Start = End: NaN in the Fortran up to
        // r1713, ignored since r1717).
        Violation{"LowZD1StartAboveEnd", [](Json& j) { j["ZOptLowStart"] = 4.0; j["ZOptLowEnd"] = 3.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowStart", "the absorber would be empty"},
        Violation{"LowZD1SpecularStartAboveEnd", [](Json& j) { j["ZSpOptLowStart"] = 4.0; j["ZSpOptLowEnd"] = 3.0; j["ZSpOptLowEnergy"] = 0.3; }, "ZSpOptLowStart", "the absorber would be empty"},
        Violation{"LowZD2StartEqualsEnd", [](Json& j) { j["ZOptLowStart"] = 3.0; j["ZOptLowEnd"] = 3.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowStart", "the absorber needs a positive length"},
        Violation{"LowZD2SpecularStartEqualsEnd", [](Json& j) { j["ZSpOptLowStart"] = 3.0; j["ZSpOptLowEnd"] = 3.0; j["ZSpOptLowEnergy"] = 0.3; }, "ZSpOptLowStart", "the absorber needs a positive length"},
        // The Fortran's plain < and >: 1e-13 beyond ZStart or Zana is rejected (.RLT. would accept it).
        Violation{"LowZStartBelowZStart", [](Json& j) { j["ZOptLowStart"] = -2.0 - 1.0e-13; j["ZOptLowEnd"] = 8.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowStart", "ZOptLowStart was set smaller than ZStart"},
        Violation{"LowZEndAboveZana", [](Json& j) { j["ZOptLowStart"] = -2.0; j["ZOptLowEnd"] = 8.0 + 1.0e-13; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowEnd", "ZOptLowEnd was set higher than Zana"},
        Violation{"LowZSpecularStartBelowZStart", [](Json& j) { j["ZSpOptLowStart"] = -2.0 - 1.0e-13; j["ZSpOptLowEnd"] = 8.0; j["ZSpOptLowEnergy"] = 0.3; }, "ZSpOptLowStart", "ZSpOptLowStart was set smaller than ZStart"},
        Violation{"LowZSpecularEndAboveZana", [](Json& j) { j["ZSpOptLowStart"] = -2.0; j["ZSpOptLowEnd"] = 8.0 + 1.0e-13; j["ZSpOptLowEnergy"] = 0.3; }, "ZSpOptLowEnd", "ZSpOptLowEnd was set higher than Zana"},
        Violation{"LowZNoEnergyOrFactor", [](Json& j) { j["ZOptLowStart"] = -2.0; j["ZOptLowEnd"] = 8.0; }, "ZOptLowEnergy", "ZOptLowFactor or ZOptLowEnergy should be present."},
        Violation{"LowZEnergyAndFactor", [](Json& j) { j["ZSpOptLowStart"] = -2.0; j["ZSpOptLowEnd"] = 8.0; j["ZSpOptLowEnergy"] = 0.3; j["ZSpOptLowFactor"] = 0.01; }, "ZSpOptLowEnergy", "give either ZSpOptLowEnergy or ZSpOptLowFactor, not both"},
        Violation{"LowZEnergyNegative", [](Json& j) { j["ZOptLowStart"] = -2.0; j["ZOptLowEnd"] = 8.0; j["ZOptLowEnergy"] = -0.1; }, "ZOptLowEnergy", "ZOptLowEnergy < 0"},
        Violation{"LowZSpecularEnergyNegative", [](Json& j) { j["ZSpOptLowStart"] = -2.0; j["ZSpOptLowEnd"] = 8.0; j["ZSpOptLowEnergy"] = -0.1; }, "ZSpOptLowEnergy", "ZSpOptLowEnergy <= 0"},
        Violation{"LowZSpecularEnergyZero", [](Json& j) { j["ZSpOptLowStart"] = -2.0; j["ZSpOptLowEnd"] = 8.0; j["ZSpOptLowEnergy"] = 0.0; }, "ZSpOptLowEnergy", "ZSpOptLowEnergy <= 0"},
        // Checks whose only violation used to trip an earlier check of the same key (or none at all).
        Violation{"RDoF", [](Json& j) { j["RDoF"] = false; }, "RDoF", "Currently only 6D calculations are supported."},
        Violation{"XDoF", [](Json& j) { j["XDoF"] = false; }, "XDoF", "Currently only 6D calculations are supported."},
        Violation{"YDoF", [](Json& j) { j["YDoF"] = false; }, "YDoF", "Currently only 6D calculations are supported."},
        Violation{"ThetaDoF", [](Json& j) { j["ThetaDoF"] = false; }, "ThetaDoF", "Currently only 6D calculations are supported."},
        Violation{"Mass2", [](Json& j) { j["Mass2"] = 0.0; }, "Mass2", "Mass of atom 2 <= 0."},
        // rle: within 1e-12 of 0 is still "<= 0" (a plain <= would let 5e-13 through to the HomoNuc check).
        Violation{"Mass2WithinTolerance", [](Json& j) { j["Mass2"] = 5.0e-13; }, "Mass2", "Mass of atom 2 <= 0."},
        Violation{"nZSpecularZero", [](Json& j) { j["nZSpecular"] = 0; }, "nZSpecular", "nZSpecular < 1."},
        Violation{"nRZero", [](Json& j) { j["nR"] = 0; }, "nR", "nR < 1."},
        Violation{"nXZero", [](Json& j) { j["nX"] = 0; }, "nX", "nX < 1."},
        Violation{"nYFft", [](Json& j) { j["nY"] = 5; }, "nY", "nY is not a valid FFT size."},
        // nX = nY on the square and hexagonal cells (C4v, C6v, C3v); the hard wall (C2v) may differ.
        Violation{"nXnYC6v", [](Json& j) { j = case_json("cu111_dailight_small"); j["nY"] = 4; }, "nY", "nX /= nY."},
        Violation{"nXnYC3v", [](Json& j) { j = case_json("crp_c3v_synthetic_small"); j["nY"] = 8; }, "nY", "nX /= nY."},
        Violation{"nXnYC4v", [](Json& j) { j = case_json("crp_c4v_synthetic_small"); j["nY"] = 6; }, "nY", "nX /= nY."},
        Violation{"tStartNegative", [](Json& j) { j["tStart"] = -10.0; j["tEnd"] = 13990.0; }, "tStart", "tStart < 0."},
        Violation{"tStartNegativeBeyondTolerance", [](Json& j) { j["tStart"] = -2.0e-12; }, "tStart", "tStart < 0."},
        Violation{"tStartGtTEnd",
                  [](Json& j) { j["NewCalc"] = false; j["RestartDir"] = "earlier"; j["tStart"] = 100.0; j["tEnd"] = 50.0; },
                  "tStart", "tStart > tEnd."},
        Violation{"RestartFormat", [](Json& j) { j["RestartFormat"] = "binary"; }, "RestartFormat",
                  "RestartFormat must be \"native\" or \"fortran\"."},
        Violation{"CheckpointFormat", [](Json& j) { j["CheckpointFormat"] = "binary"; }, "CheckpointFormat",
                  "CheckpointFormat must be \"native\", \"fortran\" or \"both\"."},
        Violation{"dtAnaZero", [](Json& j) { j["dtAna"] = 0.0; }, "dtAna", "dtAna <= 0."},
        Violation{"dtAnaWithinTolerance", [](Json& j) { j["dtAna"] = 5.0e-13; }, "dtAna", "dtAna <= 0."},
        Violation{"dtAveZero", [](Json& j) { j["dtAve"] = 0.0; }, "dtAve", "dtAve <= 0."},
        Violation{"dtAveWithinTolerance", [](Json& j) { j["dtAve"] = 5.0e-13; }, "dtAve", "dtAve <= 0."},
        // ZStart above Zana (Zana stays above ZAsymptotic and on the shifted grid).
        Violation{"ZanaLtZStart", [](Json& j) { j["ZStart"] = 8.08; }, "Zana", "ZAna < ZStart."},
        Violation{"ZanaBeyondTheGrid", [](Json& j) { j["Zana"] = 12.08; }, "Zana",
                  "ZAnalysis is too big for the scattering grid in Z."},
        Violation{"EMaxWPZero", [](Json& j) { j["EMaxWP"] = 0.0; }, "EMaxWP", "EMaxWP <= 0."},
        Violation{"EMaxWPWithinTolerance", [](Json& j) { j["EMaxWP"] = 5.0e-13; }, "EMaxWP", "EMaxWP <= 0."},
        Violation{"ROptStrengthMissing", [](Json& j) { j.erase("ROptEnergy"); }, "ROptEnergy",
                  "ROptFactor or ROptEnergy should be present."},
        Violation{"ROptEnergyNegative", [](Json& j) { j["ROptEnergy"] = -0.1; }, "ROptEnergy", "ROptEnergy <= 0."},
        Violation{"ROptEnergyZero", [](Json& j) { j["ROptEnergy"] = 0.0; }, "ROptEnergy", "ROptEnergy <= 0."},
        Violation{"ROptEnergyJustBelowZero", [](Json& j) { j["ROptEnergy"] = -1.0e-13; }, "ROptEnergy", "ROptEnergy <= 0."},
        Violation{"ZSpOptStrengthMissing", [](Json& j) { j.erase("ZSpOptEnergy"); }, "ZSpOptEnergy",
                  "ZSpOptFactor or ZSpOptEnergy should be present."},
        Violation{"ZSpOptEnergyNegative", [](Json& j) { j["ZSpOptEnergy"] = -0.1; }, "ZSpOptEnergy", "ZSpOptEnergy <= 0."},
        Violation{"ZSpOptEnergyZero", [](Json& j) { j["ZSpOptEnergy"] = 0.0; }, "ZSpOptEnergy", "ZSpOptEnergy <= 0."},
        Violation{"ZSpOptEnergyJustBelowZero", [](Json& j) { j["ZSpOptEnergy"] = -1.0e-13; }, "ZSpOptEnergy",
                  "ZSpOptEnergy <= 0."},
        Violation{"LowZSpecularNoEnergyOrFactor", [](Json& j) { j["ZSpOptLowStart"] = -2.0; j["ZSpOptLowEnd"] = 8.0; },
                  "ZSpOptLowEnergy", "ZSpOptLowFactor or ZSpOptLowEnergy should be present."},
        // Checks reached only inside the 1e-12 band of an earlier tolerant check.
        Violation{"nStepWithinTolerance",
                  [](Json& j) { j["dtStep"] = 1.5e-12; j["nStep"] = -1; j["tEnd"] = -9.0e-13; }, "nStep", "nStep < 0."},
        Violation{"ZSpOptStartBelowZAsymptotic",
                  [](Json& j) { j["ZSpOptStart"] = 7.0 - 1.5e-12; j["Z0"] = 7.0 - 8.0e-13; }, "ZSpOptStart",
                  "ZSpOptStart < ZAsymptotic."},
        Violation{"RestartDirMissing", [](Json& j) { j["NewCalc"] = false; j["tStart"] = 10.0; j["tEnd"] = 14010.0; },
                  "RestartDir", "This is a restart (NewCalc = false) but RestartDir is not given."},
        Violation{"FluxRBelowTheGrid",
                  [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 0.5; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.01; },
                  "rFluxAnalysis", "rFluxAnalysis < rmin"},
        Violation{"FluxRBeyondTheGrid",
                  [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 9.2; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.01; },
                  "rFluxAnalysis", "rFluxAnalysis > rmax"},
        Violation{"FluxEMinGeEMax",
                  [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.4; j["EMaxFluxAna"] = 0.3; j["DeltaEFluxAna"] = 0.01; },
                  "EMinFluxAna", "EMinFluxAna >= EMaxFluxAna"},
        // The tolerant siblings of the boundary rows above, and the exact checks just below 0 (r1718).
        Violation{"Mass1WithinTolerance", [](Json& j) { j["Mass1"] = 5.0e-13; }, "Mass1", "Mass of atom 1 <= 0."},
        Violation{"ZOptEnergyZero", [](Json& j) { j["ZOptEnergy"] = 0.0; }, "ZOptEnergy", "ZOptEnergy <= 0."},
        Violation{"ZOptEnergyJustBelowZero", [](Json& j) { j["ZOptEnergy"] = -1.0e-13; }, "ZOptEnergy", "ZOptEnergy <= 0."},
        Violation{"LowZEnergyBeyondTolerance",
                  [](Json& j) { j["ZOptLowStart"] = -2.0; j["ZOptLowEnd"] = 8.0; j["ZOptLowEnergy"] = -2.0e-12; },
                  "ZOptLowEnergy", "ZOptLowEnergy < 0"},
        Violation{"LowZSpecularEnergyJustBelowZero",
                  [](Json& j) { j["ZSpOptLowStart"] = -2.0; j["ZSpOptLowEnd"] = 8.0; j["ZSpOptLowEnergy"] = -1.0e-13; },
                  "ZSpOptLowEnergy", "ZSpOptLowEnergy <= 0"}),
    [](const ::testing::TestParamInfo<Violation>& param_info) { return std::string(param_info.param.name); });

// The valid side of the tolerant checks above: within 1e-12 of the bound an rlt check accepts (a plain <
// would refuse), and the hard wall (C2v) may have nX /= nY.
TEST(Input, TolerantChecksAcceptWithinTheirTolerance) {
  Json json = reference_json();
  json["tStart"] = -5.0e-13;
  EXPECT_NO_THROW((void)make_spatial_grids(build(json)));
  json = reference_json();  // ZOptLowEnergy keeps the Fortran's tolerant .RLT. 0 (its siblings are exact)
  json["ZOptLowStart"] = -2.0;
  json["ZOptLowEnd"] = 8.0;
  json["ZOptLowEnergy"] = -5.0e-13;
  EXPECT_NO_THROW((void)make_spatial_grids(build(json)));
  json = reference_json();
  json["nY"] = 8;
  EXPECT_NO_THROW((void)make_spatial_grids(build(json)));
}

// The reader's own errors: invalid JSON (including a number that overflows a double), a JSON value that is
// not an object, a file that cannot be opened, an integer beyond the int range.
TEST(Input, ReaderErrors) {
  for (const std::string text : {"{\"nZ\": 176", R"({"dZ": 1e400})"}) {
    try {
      (void)parse_text(text);
      ADD_FAILURE() << "accepted " << text;
    } catch (const InputError& e) {
      EXPECT_NE(std::string(e.what()).find("invalid JSON in 'test.json'"), std::string::npos) << e.what();
    }
  }
  EXPECT_THROW((void)parse_input(Json::array({1})), InputError);
  try {
    (void)read_input_file("/nonexistent/input.json");
    ADD_FAILURE() << "read a file that does not exist";
  } catch (const InputError& e) {
    EXPECT_NE(std::string(e.what()).find("cannot open input file"), std::string::npos) << e.what();
  }
  Json json = reference_json();
  json["nZ"] = 10000000000LL;
  expect_input_error(json, "nZ", "integer out of range");
}

// run_program builds the PES once and passes it to build_config; the one-argument form builds it itself.
// Both give the same system parameters; derived.scm is set only for the SCM PES.
TEST(Input, BuildConfigWithAGivenPotentialEqualsTheOneArgumentForm) {
  for (const Json& json : {reference_json(), scm_json()}) {
    const InputParameters input = parse_input(json);
    const Config one = build_config(input);
    const Config two = build_config(input, *make_potential(input));
    EXPECT_EQ(two.derived.system.point_group, one.derived.system.point_group);
    EXPECT_EQ(two.derived.system.lattice_constant_x, one.derived.system.lattice_constant_x);
    EXPECT_EQ(two.derived.system.lattice_constant_y, one.derived.system.lattice_constant_y);
    EXPECT_EQ(two.derived.system.z_asymptotic, one.derived.system.z_asymptotic);
    EXPECT_EQ(two.derived.system.skewing_angle, one.derived.system.skewing_angle);
    EXPECT_EQ(two.derived.z_asymptotic_unswitched, one.derived.z_asymptotic_unswitched);
    EXPECT_EQ(two.derived.z_asymptotic, one.derived.z_asymptotic);
    ASSERT_EQ(two.derived.scm.has_value(), one.derived.scm.has_value());
    if (two.derived.scm) {
      EXPECT_EQ(two.derived.scm->slab_tag, one.derived.scm->slab_tag);
    }
  }
  EXPECT_FALSE(build(reference_json()).derived.scm);
  EXPECT_TRUE(build(scm_json()).derived.scm);
}

TEST(Input, DerivedQuantities) {
  const Config config = build(reference_json());
  const DerivedParameters& d = config.derived;
  const auto& u = units::legacy;
  EXPECT_EQ(d.total_mass, (1.0078 + 1.0078) * 1822.89);
  EXPECT_EQ(d.reduced_mass, 1.0078 * 1.0078 / (1.0078 + 1.0078) * 1822.89);
  EXPECT_EQ(d.dx, 6.239660 / 4.0);
  EXPECT_EQ(d.v_in_index, 1);
  // Wave packet in k: kmin..kmax maps to EMinWP..EMaxWP exactly at the 2 sqrt(alpha) edges.
  EXPECT_NEAR(d.confidence_e_min, 0.2, 1e-12);
  EXPECT_NEAR(d.confidence_e_max, 0.6, 1e-12);
  EXPECT_EQ(d.theta_incidence, 0.0);
  EXPECT_EQ(d.kz, d.k_total);
  EXPECT_DOUBLE_EQ(d.max_analysis_energy, u.hartree_to_ev(3.141592653589793 / 10.0));
}

TEST(Input, OpticalPotentialStrengthInterpolates) {
  // reduced length exactly between table points 2.0 and 2.5 -> average of 2.51 and 2.38.
  const double mass = 1000.0;
  const double ekin = 0.01;
  const double wavelength = 2.0 * std::acos(-1.0) / std::sqrt(2.0 * mass * ekin);
  EXPECT_NEAR(optical_potential_strength(ekin, 2.25 * wavelength, mass), 0.5 * (2.51 + 2.38) * ekin, 1e-15);
  EXPECT_EQ(optical_potential_strength(ekin, 0.5 * wavelength, mass), 2.75 * ekin);
  EXPECT_EQ(optical_potential_strength(ekin, 20.0 * wavelength, mass), 0.663 * ekin);
}

// absorber_reduced_length gives the L/lambda the strength is interpolated at, and says where it
// lies: below or above the table (1 .. 15; the Fortran's two warnings), inside, or exactly at 1.
TEST(Input, AbsorberReducedLengthNamesTheTableEnd) {
  const double mass = 1000.0;
  const double ekin = 0.01;
  const double wavelength = 2.0 * std::acos(-1.0) / std::sqrt(2.0 * mass * ekin);
  const AbsorberReducedLength below = absorber_reduced_length(ekin, 0.5 * wavelength, mass);
  EXPECT_EQ(below.value, 0.5 * wavelength / wavelength);
  EXPECT_EQ(below.range, AbsorberTableRange::below);
  EXPECT_EQ(absorber_reduced_length(ekin, 20.0 * wavelength, mass).range, AbsorberTableRange::above);
  EXPECT_EQ(absorber_reduced_length(ekin, wavelength, mass).range, AbsorberTableRange::at_first_point);
  // Not 15.0 itself: 15·λ/λ may round above 15 in the fast-math build.
  for (const double table_point : {2.25, 14.5})
    EXPECT_EQ(absorber_reduced_length(ekin, table_point * wavelength, mass).range, AbsorberTableRange::inside)
        << table_point;
}

// absorber D2: L/lambda exactly 1 matches no interval of the table: strength 0 (the Fortran returns
// its uninitialised local); the set-up warns that the absorber is off (Report test).
TEST(Input, AbsorberD2ReducedLengthOneGivesStrengthZero) {
  const double mass = 1000.0;
  const double ekin = 0.01;
  const double wavelength = 2.0 * std::acos(-1.0) / std::sqrt(2.0 * mass * ekin);
  const AbsorberReducedLength at_one = absorber_reduced_length(ekin, wavelength, mass);
  ASSERT_EQ(at_one.value, 1.0);
  EXPECT_EQ(at_one.range, AbsorberTableRange::at_first_point);
  EXPECT_EQ(optical_potential_strength(ekin, at_one), 0.0);
  EXPECT_EQ(optical_potential_strength(ekin, wavelength, mass), 0.0);
}

// Every absorber tuned by its energy (all five energy keys) whose L/lambda lies outside the table
// gets the end value of the table and is recorded for the set-up warning, in the order of the
// Fortran's ProcessInput (Fortran: getACoeffOptPotential warns). For ROpt and ZSpOpt the Fortran
// 5.24 printed reduced_length = 0.58847038336148938 and 15.945765873846009 for this input: the
// C++ values equal these doubles bitwise in the default build (2026-11-15); the tolerance 1e-13,
// about 50 ulp at 16, leaves room for the reassociations of the fast-math build.
TEST(Input, AbsorberOutsideTheStrengthTableIsRecorded) {
  Json json = reference_json();
  EXPECT_TRUE(build(json).derived.absorber_table_warnings.empty());

  json["ZOptEnergy"] = 5.0;      // 3.12 bohr: L/lambda = 18.2
  json["ZOptLowStart"] = -1.0;   // 0.5 bohr at 0.3 eV: L/lambda = 0.72
  json["ZOptLowEnd"] = -0.5;
  json["ZOptLowEnergy"] = 0.3;
  json["ROptEnergy"] = 0.01;     // 4.5 bohr with the reduced mass: L/lambda = 0.59
  json["ZSpOptEnergy"] = 3.0;    // 3.52 bohr: L/lambda = 15.9
  json["ZSpOptLowStart"] = -1.9; // 0.5 bohr at 0.3 eV: L/lambda = 0.72
  json["ZSpOptLowEnd"] = -1.4;
  json["ZSpOptLowEnergy"] = 0.3;
  const Config config = build(json);
  const DerivedParameters& d = config.derived;
  const auto& u = units::legacy;
  const std::vector<AbsorberTableWarning>& warnings = d.absorber_table_warnings;
  ASSERT_EQ(warnings.size(), 5);
  const std::vector<std::string> keys = {"ZOptEnergy", "ZOptLowEnergy", "ROptEnergy", "ZSpOptEnergy",
                                         "ZSpOptLowEnergy"};
  const std::vector<AbsorberTableRange> ranges = {AbsorberTableRange::above, AbsorberTableRange::below,
                                                  AbsorberTableRange::below, AbsorberTableRange::above,
                                                  AbsorberTableRange::below};
  for (int i = 0; i < 5; ++i) {
    EXPECT_EQ(warnings[i].energy_key, keys[i]) << i;
    EXPECT_EQ(warnings[i].range, ranges[i]) << keys[i];
  }
  EXPECT_NEAR(warnings[2].reduced_length, 0.58847038336148938, 1e-13);
  EXPECT_NEAR(warnings[3].reduced_length, 15.945765873846009, 1e-13);
  // The table's end values; within 4 ulp because the fast-math build rounds ekin differently
  // (1 ulp for ZOptEnergy 5 eV, 2026-11-15).
  EXPECT_DOUBLE_EQ(d.z_opt_factor, 0.663 * u.ev_to_hartree(5.0));
  EXPECT_DOUBLE_EQ(d.z_opt_low_factor, 2.75 * u.ev_to_hartree(0.3));
  EXPECT_DOUBLE_EQ(d.r_opt_factor, 2.75 * u.ev_to_hartree(0.01));
  EXPECT_DOUBLE_EQ(d.zsp_opt_factor, 0.663 * u.ev_to_hartree(3.0));
  EXPECT_DOUBLE_EQ(d.zsp_opt_low_factor, 2.75 * u.ev_to_hartree(0.3));

  json["ZOptEnergy"] = 0.3;
  json.erase("ROptEnergy");  // a factor: no strength computed from an energy, so nothing to warn about
  json["ROptFactor"] = 0.01;
  json["ZSpOptEnergy"] = 0.3;
  for (const std::string key :
       {"ZOptLowStart", "ZOptLowEnd", "ZOptLowEnergy", "ZSpOptLowStart", "ZSpOptLowEnd", "ZSpOptLowEnergy"})
    json.erase(key);
  EXPECT_TRUE(build(json).derived.absorber_table_warnings.empty());
}

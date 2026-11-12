// Tests: JSON reading, every input consistency check (valid and violating case), derived parameters.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>

#include <nlohmann/json.hpp>

#include "core/error.hpp"
#include "core/unit_conversion.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "input/wave_packet_params.hpp"

using namespace spodvr;
using Json = nlohmann::json;

namespace {

Json reference_json() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  return Json::parse(file);
}

Config build(const Json& json) { return build_config(parse_input(json)); }

// Expects an InputError naming `parameter` for the modified input.
void expect_input_error(const Json& json, const std::string& parameter) {
  try {
    const Config config = build(json);
    (void)make_spatial_grids(config);
    ADD_FAILURE() << "expected an InputError for " << parameter << " with input " << json.dump();
  } catch (const InputError& e) {
    EXPECT_EQ(e.parameter(), parameter) << e.what();
  }
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
  expect_input_error(json, "NotAKey");
}

TEST(Input, MissingRequiredKeyIsAnError) {
  for (const char* key : {"Title", "nZ", "TolSp", "writeRDProbs", "ZSpOptEnd"}) {
    Json json = reference_json();
    json.erase(key);
    expect_input_error(json, key);
  }
}

TEST(Input, WrongTypesAreErrors) {
  Json json = reference_json();
  json["nZ"] = 176.0;
  expect_input_error(json, "nZ");
  json = reference_json();
  json["NewCalc"] = 1;
  expect_input_error(json, "NewCalc");
  json = reference_json();
  json["dZ"] = "0.08";
  expect_input_error(json, "dZ");
}

TEST(Input, IntegerAcceptedForRealKey) {
  Json json = reference_json();
  json["ZStart"] = -2;
  EXPECT_EQ(build(json).input.z_start, -2.0);
}

TEST(Input, UnsupportedFeatures) {
  for (const char* key : {"SystemMix", "UseRCutOff", "write2DPES"}) {
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

// CRPlibDataSet (directory of the CRP data set, relative to the working directory): required with
// SystemID CRPlibrary, an error with any other PES. build_config reads the data set (the PES gives
// the system parameters), so a missing one is reported there; the data-set checks themselves are
// in tests/crplib and test_crp_potential.
TEST(Input, CrpDataSetKey) {
  Json json = reference_json();  // hard wall
  json["CRPlibDataSet"] = "some/data/set";
  expect_input_error(json, "CRPlibDataSet");

  json = reference_json();
  json["SystemID"] = "CRPlibrary";
  expect_input_error(json, "CRPlibDataSet");  // missing
  json["CRPlibDataSet"] = "";
  expect_input_error(json, "CRPlibDataSet");
  json["CRPlibDataSet"] = "no/such/data/set";
  expect_input_error(json, "CRPlibDataSet");

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
  expect_input_error(json, "SymmetryTolerance");
  json["SymmetryTolerance"] = -1e-5;
  expect_input_error(json, "SymmetryTolerance");
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
  expect_input_error(json, "ZOptEnergy");
  json["ZOptFactor"] = 0.01;
  EXPECT_EQ(build(json).derived.z_opt_factor, 0.01);
  json["ZOptEnergy"] = 0.3;
  expect_input_error(json, "ZOptEnergy");  // both given
  json.erase("ZOptFactor");
  json["ZOptEnergy"] = -0.1;
  expect_input_error(json, "ZOptEnergy");
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
  expect_input_error(json, "ZminSwitch");
}

// flux Q1: a flux energy range that ends at or below EMaxWP is valid (the Fortran check rejects
// every such range); the wave packet spans 0.2..0.6 eV.
TEST(Input, FluxEnergyRangeInsideTheWavePacketIsValid) {
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

// One violation per check of ValidateInput / ReadInput / Setup.f90.
struct Violation {
  const char* name;
  std::function<void(Json&)> modify;
  const char* parameter;
};

void PrintTo(const Violation& v, std::ostream* os) { *os << v.name; }

class InputViolation : public ::testing::TestWithParam<Violation> {};

TEST_P(InputViolation, IsRejected) {
  Json json = reference_json();
  GetParam().modify(json);
  expect_input_error(json, GetParam().parameter);
}

INSTANTIATE_TEST_SUITE_P(
    Checks, InputViolation,
    ::testing::Values(
        Violation{"SystemID", [](Json& j) { j["SystemID"] = "H2Pt111"; }, "SystemID"},
        Violation{"ZDoF", [](Json& j) { j["ZDoF"] = false; }, "ZDoF"},
        Violation{"PhiDoF", [](Json& j) { j["PhiDoF"] = false; }, "PhiDoF"},
        Violation{"Mass1", [](Json& j) { j["Mass1"] = 0.0; j["Mass2"] = 0.0; }, "Mass1"},
        Violation{"HomoMasses", [](Json& j) { j["Mass2"] = 2.0; }, "HomoNuc"},
        Violation{"HeteroMasses", [](Json& j) { j["HomoNuc"] = false; }, "HomoNuc"},
        Violation{"nZFft", [](Json& j) { j["nZ"] = 177; }, "nZ"},
        Violation{"nZ", [](Json& j) { j["nZ"] = 0; }, "nZ"},
        Violation{"dZ", [](Json& j) { j["dZ"] = 0.0; }, "dZ"},
        Violation{"nZSpecularFft", [](Json& j) { j["nZSpecular"] = 281; }, "nZSpecular"},
        Violation{"nZSpecularSmall", [](Json& j) { j["nZSpecular"] = 160; j["ZSpOptStart"] = 8.0; j["ZSpOptEnd"] = 10.0; j["Z0"] = 7.5; }, "nZSpecular"},
        Violation{"RStart", [](Json& j) { j["RStart"] = -0.1; }, "RStart"},
        // RStart <= 0 within 1e-12; ROptEnd moves onto the shifted r grid, so only RStart is violated.
        Violation{"RStartZero", [](Json& j) { j["RStart"] = 0.0; j["ROptEnd"] = 8.25; }, "RStart"},
        Violation{"RStartWithinTolerance", [](Json& j) { j["RStart"] = 1.0e-13; j["ROptEnd"] = 8.25; }, "RStart"},
        Violation{"nR", [](Json& j) { j["nR"] = 57; }, "nR"},
        Violation{"dR", [](Json& j) { j["dR"] = -0.15; }, "dR"},
        Violation{"nX", [](Json& j) { j["nX"] = 5; }, "nX"},
        Violation{"nY", [](Json& j) { j["nY"] = 0; }, "nY"},
        Violation{"JIn", [](Json& j) { j["JIn"] = -2; }, "JIn"},
        Violation{"JMax", [](Json& j) { j["JMax"] = -1; j["MJMax"] = 0; }, "JMax"},
        Violation{"MJMax", [](Json& j) { j["MJMax"] = 6; }, "MJMax"},
        Violation{"MJIn", [](Json& j) { j["MJIn"] = 1; }, "MJIn"},
        Violation{"JInGtJMax", [](Json& j) { j["JIn"] = 6; }, "JIn"},
        Violation{"MJInGtMJMax", [](Json& j) { j["JIn"] = 2; j["MJIn"] = 2; j["MJMax"] = 0; }, "MJIn"},
        Violation{"ParityJ", [](Json& j) { j["JMax"] = 5; j["MJMax"] = 5; }, "JMax"},
        Violation{"ParityMJ", [](Json& j) { j["MJMax"] = 3; }, "MJMax"},
        Violation{"ProjOp0", [](Json& j) { j["ProjOp"] = 0; }, "ProjOp"},
        Violation{"ProjOp1", [](Json& j) { j["ProjOp"] = 1; }, "ProjOp"},
        Violation{"ProjOp3", [](Json& j) { j["ProjOp"] = 3; }, "ProjOp"},
        Violation{"tStartNew", [](Json& j) { j["tStart"] = 10.0; j["tEnd"] = 14010.0; }, "tStart"},
        Violation{"tStartRestart", [](Json& j) { j["NewCalc"] = false; }, "tStart"},
        Violation{"OnlyAnalysisNeedsRestartDir", [](Json& j) { j["onlyAnalysis"] = true; }, "RestartDir"},
        Violation{"tEnd", [](Json& j) { j["tEnd"] = 14002.0; }, "tEnd"},
        Violation{"dtStep", [](Json& j) { j["dtStep"] = 0.0; j["nStep"] = 0; j["tEnd"] = 0.0; }, "dtStep"},
        Violation{"nStep", [](Json& j) { j["nStep"] = -1; j["tEnd"] = -2.0; }, "tEnd"},
        Violation{"dtAna", [](Json& j) { j["dtAna"] = 1.0; }, "dtAna"},
        Violation{"dtAnaInteger", [](Json& j) { j["dtAna"] = 30.0; j["dtAve"] = 90.0; j["dtStep"] = 2.0; j["nStep"] = 7001; j["tEnd"] = 14002.0; }, "dtAna"},
        Violation{"dtAve", [](Json& j) { j["dtAve"] = 1.0; }, "dtAve"},
        Violation{"dtAveLtAna", [](Json& j) { j["dtAve"] = 4.0; }, "dtAve"},
        Violation{"dtAveInteger", [](Json& j) { j["dtAve"] = 25.0; }, "dtAve"},
        Violation{"JAnaNeg", [](Json& j) { j["JAna"] = -1; }, "JAna"},
        Violation{"VAnaNeg", [](Json& j) { j["VAna"] = -1; }, "VAna"},
        Violation{"DifOrdNeg", [](Json& j) { j["DifOrdAna"] = -1; }, "DifOrdAna"},
        Violation{"NrOfEs", [](Json& j) { j["NrOfEs"] = -1; }, "NrOfEs"},
        Violation{"ZanaLtZasymp", [](Json& j) { j["Zana"] = 6.0; }, "Zana"},
        Violation{"ZanaOffGrid", [](Json& j) { j["Zana"] = 8.01; }, "Zana"},
        Violation{"ZanaGtZOptStart", [](Json& j) { j["Zana"] = 9.04; }, "Zana"},
        Violation{"JAnaGtJMax", [](Json& j) { j["JAna"] = 6; }, "JAna"},
        Violation{"VAnaGtNr", [](Json& j) { j["VAna"] = 56; }, "VAna"},
        Violation{"DifOrdGt", [](Json& j) { j["DifOrdAna"] = 7; }, "DifOrdAna"},
        Violation{"dtAnaFine", [](Json& j) { j["EMaxWP"] = 9.0; }, "dtAna"},
        Violation{"FluxMissing", [](Json& j) { j["doFluxAnalysis"] = true; }, "rFluxAnalysis"},
        Violation{"FluxOffGrid", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.01; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.01; }, "rFluxAnalysis"},
        Violation{"FluxEMin", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.1; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.01; }, "EMinFluxAna"},
        // flux Q1: EMaxFluxAna above EMaxWP (0.6 eV) is rejected.
        Violation{"FluxEMax", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.7; j["DeltaEFluxAna"] = 0.01; }, "EMaxFluxAna"},
        // flux Q2: a zero or negative energy spacing is rejected (the Fortran accepts 0).
        Violation{"FluxDeltaZero", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.0; }, "DeltaEFluxAna"},
        Violation{"FluxDeltaNegative", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = -0.01; }, "DeltaEFluxAna"},
        Violation{"FluxDelta", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 3.05; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.5; }, "DeltaEFluxAna"},
        Violation{"FluxROpt", [](Json& j) { j["doFluxAnalysis"] = true; j["rFluxAnalysis"] = 5.0; j["EMinFluxAna"] = 0.2; j["EMaxFluxAna"] = 0.6; j["DeltaEFluxAna"] = 0.01; }, "rFluxAnalysis"},
        Violation{"EMinWP", [](Json& j) { j["EMinWP"] = 0.0; }, "EMinWP"},
        Violation{"EMinGeEMax", [](Json& j) { j["EMinWP"] = 0.6; }, "EMinWP"},
        Violation{"Z0LtZasymp", [](Json& j) { j["Z0"] = 6.0; }, "Z0"},
        Violation{"Z0GtZSpOpt", [](Json& j) { j["Z0"] = 17.0; }, "Z0"},
        Violation{"VIn", [](Json& j) { j["VIn"] = -1; }, "VIn"},
        Violation{"VInGtNr", [](Json& j) { j["VIn"] = 56; }, "VIn"},
        Violation{"VInGtVAna", [](Json& j) { j["VIn"] = 2; j["VAna"] = 1; }, "VIn"},
        Violation{"JInGtJAna", [](Json& j) { j["JIn"] = 2; j["JAna"] = 0; }, "JIn"},
        Violation{"VMinGtVMax", [](Json& j) { j["ImposeVMin"] = true; j["ImposeVMax"] = true; j["VMin"] = 2.0; j["VMax"] = 1.0; }, "VMin"},
        Violation{"ZOptOrder", [](Json& j) { j["ZOptEnd"] = 8.5; }, "ZOptStart"},
        Violation{"ZOptEnd", [](Json& j) { j["ZOptEnd"] = 13.0; }, "ZOptEnd"},
        Violation{"ROptOrder", [](Json& j) { j["ROptEnd"] = 4.0; }, "ROptStart"},
        Violation{"ROptEnd", [](Json& j) { j["ROptEnd"] = 9.2; }, "ROptEnd"},
        Violation{"ZSpOptOrder", [](Json& j) { j["ZSpOptEnd"] = 16.0; }, "ZSpOptStart"},
        Violation{"ZSpOptEnd", [](Json& j) { j["ZSpOptEnd"] = 21.0; }, "ZSpOptEnd"},
        // ZSpOptStart < ZAsymptotic cannot be reached on its own: Z0 >= ZAsymptotic and
        // Z0 <= ZSpOptStart are checked first.
        Violation{"ZSpOptStart", [](Json& j) { j["ZSpOptStart"] = 6.0; }, "Z0"},
        Violation{"EParallel", [](Json& j) { j["EParallel"] = -0.1; }, "EParallel"},
        Violation{"PhiParallelNeg", [](Json& j) { j["PhiParallel"] = -1.0; }, "PhiParallel"},
        Violation{"PhiParallelBig", [](Json& j) { j["PhiParallel"] = 361.0; }, "PhiParallel"},
        Violation{"ZAnaSpDynamic", [](Json& j) { j["ZAnaSp"] = 12.0; }, "ZAnaSp"},
        Violation{"ZAnaSpGtOpt", [](Json& j) { j["DynamZSpAna"] = false; j["ZAnaSp"] = 17.0; }, "ZAnaSp"},
        Violation{"ZAnaSpLtZ0", [](Json& j) { j["DynamZSpAna"] = false; j["ZAnaSp"] = 11.0; }, "ZAnaSp"},
        // Low-Z absorbers. low-Z D3: an incomplete set of keys (the Fortran ignores it).
        Violation{"LowZD3OnlyStart", [](Json& j) { j["ZOptLowStart"] = -2.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowEnd"},
        Violation{"LowZD3OnlyEnd", [](Json& j) { j["ZOptLowEnd"] = 8.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowStart"},
        Violation{"LowZD3OnlyStrength", [](Json& j) { j["ZSpOptLowFactor"] = 0.01; }, "ZSpOptLowStart"},
        // low-Z D1 (Start > End, ignored by the Fortran) and low-Z D2 (Start = End, NaN there).
        Violation{"LowZD1StartAboveEnd", [](Json& j) { j["ZOptLowStart"] = 4.0; j["ZOptLowEnd"] = 3.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowStart"},
        Violation{"LowZD1SpecularStartAboveEnd", [](Json& j) { j["ZSpOptLowStart"] = 4.0; j["ZSpOptLowEnd"] = 3.0; j["ZSpOptLowEnergy"] = 0.3; }, "ZSpOptLowStart"},
        Violation{"LowZD2StartEqualsEnd", [](Json& j) { j["ZOptLowStart"] = 3.0; j["ZOptLowEnd"] = 3.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowStart"},
        Violation{"LowZD2SpecularStartEqualsEnd", [](Json& j) { j["ZSpOptLowStart"] = 3.0; j["ZSpOptLowEnd"] = 3.0; j["ZSpOptLowEnergy"] = 0.3; }, "ZSpOptLowStart"},
        // The Fortran's plain < and >: 1e-13 beyond ZStart or Zana is rejected (.RLT. would accept it).
        Violation{"LowZStartBelowZStart", [](Json& j) { j["ZOptLowStart"] = -2.0 - 1.0e-13; j["ZOptLowEnd"] = 8.0; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowStart"},
        Violation{"LowZEndAboveZana", [](Json& j) { j["ZOptLowStart"] = -2.0; j["ZOptLowEnd"] = 8.0 + 1.0e-13; j["ZOptLowEnergy"] = 0.3; }, "ZOptLowEnd"},
        Violation{"LowZSpecularStartBelowZStart", [](Json& j) { j["ZSpOptLowStart"] = -2.0 - 1.0e-13; j["ZSpOptLowEnd"] = 8.0; j["ZSpOptLowEnergy"] = 0.3; }, "ZSpOptLowStart"},
        Violation{"LowZSpecularEndAboveZana", [](Json& j) { j["ZSpOptLowStart"] = -2.0; j["ZSpOptLowEnd"] = 8.0 + 1.0e-13; j["ZSpOptLowEnergy"] = 0.3; }, "ZSpOptLowEnd"},
        Violation{"LowZNoEnergyOrFactor", [](Json& j) { j["ZOptLowStart"] = -2.0; j["ZOptLowEnd"] = 8.0; }, "ZOptLowEnergy"},
        Violation{"LowZEnergyAndFactor", [](Json& j) { j["ZSpOptLowStart"] = -2.0; j["ZSpOptLowEnd"] = 8.0; j["ZSpOptLowEnergy"] = 0.3; j["ZSpOptLowFactor"] = 0.01; }, "ZSpOptLowEnergy"},
        Violation{"LowZEnergyNegative", [](Json& j) { j["ZOptLowStart"] = -2.0; j["ZOptLowEnd"] = 8.0; j["ZOptLowEnergy"] = -0.1; }, "ZOptLowEnergy"},
        Violation{"LowZSpecularEnergyNegative", [](Json& j) { j["ZSpOptLowStart"] = -2.0; j["ZSpOptLowEnd"] = 8.0; j["ZSpOptLowEnergy"] = -0.1; }, "ZSpOptLowEnergy"}),
    [](const ::testing::TestParamInfo<Violation>& param_info) { return std::string(param_info.param.name); });

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

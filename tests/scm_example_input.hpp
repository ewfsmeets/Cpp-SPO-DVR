#pragma once
// The SCM example of the repository as a program input: the SRP48 data set with its example slab and
// parameter files (potentials/H2Cu111SRP_C3v.dir), on the small grid of crp_cu111_srp48_small from
// Z = 0.5 bohr, the first point above the CRP's Zmin = 0.47 bohr. Paths are absolute.

#include <filesystem>
#include <fstream>
#include <string>

#include <nlohmann/json.hpp>

namespace test_support {

inline const std::filesystem::path scm_example_directory =
    std::filesystem::path(SPODVR_REFERENCE_DIR).parent_path().parent_path() / "potentials/H2Cu111SRP_C3v.dir";

// The input with this position file: a full run of the CRP case (18000 steps); set "Potential" or the
// times as needed.
inline nlohmann::json scm_example_json(const std::filesystem::path& position_file) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/crp_cu111_srp48_small/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  const std::filesystem::path data_set =
      scm_example_directory.parent_path().parent_path() / json["CRPlibDataSet"].get<std::string>();
  json.update({{"SystemID", "CRPDistrtSurface"},
               {"CRPlibDataSet", data_set.string()},
               {"ZStart", 0.5},
               {"nZ", 110},
               {"nZSpecular", 198},
               {"ZSpOptEnd", 30.05},
               {"SCMParFilepath", (scm_example_directory / "params.dat").string()},
               {"PosFilepath", position_file.string()}});
  return json;
}

}  // namespace test_support

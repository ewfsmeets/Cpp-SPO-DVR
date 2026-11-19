// Tests: potential/ (writePotential, report/potential_files.*): the files hold V as the propagation reads it
// (bitwise), the symmetric case agrees with the potential computed on the full grid, the minimum's orientation and
// the sites agree with each other, a NaN propagates, and a potential-only run writes the files and stops.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "app/cli.hpp"
#include "app/run.hpp"
#include "comparison.hpp"
#include "core/error.hpp"
#include "grids/angular_grid.hpp"
#include "grids/spatial_grids.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "input/validate_outputs.hpp"
#include "io/binary_file.hpp"
#include "potential/dai_light_cu111_potential.hpp"
#include "potential/potential_grid.hpp"
#include "report/plot_data.hpp"
#include "report/potential_files.hpp"
#include "temporary_directory.hpp"

using namespace spodvr;
using ::test_support::TemporaryDirectory;
namespace fs = std::filesystem;

namespace {

using Json = nlohmann::json;

Json case_json(const std::string& name) {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/" + name + "/input.json");
  return Json::parse(file);
}

// The Dai-Light C6v case (UseSymmetry, 8 x 8) with writePotential and three sites on its grid.
Json symmetric_input() {
  Json json = case_json("cu111_dailight_short_sym");
  json["writePotential"] = true;
  json["PlotSites"] = Json::array({{{"name", "top"}, {"x", 0.0}, {"y", 0.0}},
                                   {{"name", "bridge"}, {"x", 0.5}, {"y", 0.0}},
                                   {{"name", "off"}, {"x", 0.375}, {"y", 0.125}}});
  return json;
}

struct PotentialSetup {
  Config config;
  SpatialGrids grids;
  AngularGrid angular;
  PotentialGrid reduced;  // with the run's surface symmetry
  PotentialGrid full;     // the same grid without symmetry: every (x, y) point computed
  explicit PotentialSetup(const Json& json)
      : config(build_config(parse_input(json))),
        grids(make_spatial_grids(config)),
        angular(make_angular_grid(config)),
        reduced(compute_potential_grid(DaiLightCu111Potential(), config, grids, angular)),
        full(compute_potential_grid(DaiLightCu111Potential(), config, grids, angular,
                                    GridSymmetry(config.input.nx, config.input.ny, angular.n_phi()))) {}
};

const PotentialSetup& symmetric() {
  static const PotentialSetup s(symmetric_input());
  return s;
}

// Every value of a little-endian binary file of T.
template <typename T>
std::vector<T> read_values(const fs::path& path) {
  std::vector<T> values(fs::file_size(path) / sizeof(T));
  BinaryReader reader(path);
  reader.read_array(std::span(values));
  return values;
}

// The written files of the symmetric case (written once).
struct Written {
  std::vector<double> v_min;
  std::vector<std::int32_t> v_min_node;
  std::vector<double> v_sites;
};

const Written& written() {
  static const Written w = [] {
    TemporaryDirectory dir("potential_files");
    const PotentialSetup& s = symmetric();
    write_potential_files(dir.path(), s.config, s.grids, s.angular, s.reduced);
    return Written{read_values<double>(dir.path() / "potential" / "v_min.bin"),
                   read_values<std::int32_t>(dir.path() / "potential" / "v_min_node.bin"),
                   read_values<double>(dir.path() / "potential" / "v_sites.bin")};
  }();
  return w;
}

// |a - b| relative to max(|a|, |b|, 1 hartree) (NaN if either is).
double scaled_difference(double a, double b) { return std::abs(a - b) / std::max({std::abs(a), std::abs(b), 1.0}); }

// The minimum over the orientations of a grid at (ixy, izr), NaN-propagating, by a plain loop.
double minimum_by_loop(const PotentialGrid& grid, std::size_t ixy, std::size_t izr) {
  double lowest = std::numeric_limits<double>::infinity();
  for (std::size_t itp = 0; itp < grid.index.n0; ++itp) {
    const double v = grid(itp, ixy, izr);
    if (std::isnan(v)) return v;
    lowest = std::min(lowest, v);
  }
  return lowest;
}

// A porting slip: the gather through the symmetry maps without the phi image of the operation.
double gather_without_phi_map(const PotentialGrid& grid, std::size_t itp, std::size_t ixy, std::size_t izr) {
  return grid.values[grid.index(itp, grid.symmetry.representative_of(static_cast<int>(ixy)), izr)];
}

}  // namespace

// v_min and v_sites equal V as the propagation reads it (PotentialGrid::operator()), bitwise; v_min_node is
// the first orientation of the minimum in itp order.
TEST(PotentialFiles, HoldThePotentialAsThePropagationReadsIt) {
  const PotentialSetup& s = symmetric();
  const Written& w = written();
  const std::size_t n_xy = s.reduced.symmetry.n_xy();
  const std::size_t n_zr = s.reduced.index.n2;
  const std::size_t n_tp = s.reduced.index.n0;
  const int n_phi = s.angular.n_phi();
  ASSERT_EQ(w.v_min.size(), n_zr * n_xy);
  ASSERT_EQ(w.v_min_node.size(), 2 * n_zr * n_xy);
  int differing = 0;
  int wrong_nodes = 0;
  for (std::size_t izr = 0; izr < n_zr; ++izr)
    for (std::size_t ixy = 0; ixy < n_xy; ++ixy) {
      const double expected = minimum_by_loop(s.reduced, ixy, izr);
      const std::size_t k = ixy + n_xy * izr;
      differing += (w.v_min[k] == expected) ? 0 : 1;
      const std::size_t itp = w.v_min_node[2 * k + 1] + n_phi * w.v_min_node[2 * k];
      std::size_t first = 0;
      while (first < n_tp && s.reduced(first, ixy, izr) != expected) ++first;
      wrong_nodes += (itp == first) ? 0 : 1;
    }
  EXPECT_EQ(differing, 0);
  EXPECT_EQ(wrong_nodes, 0);

  const auto& sites = s.config.input.plot_sites;
  ASSERT_EQ(w.v_sites.size(), sites.size() * n_zr * n_tp);
  int differing_sites = 0;
  for (std::size_t isite = 0; isite < sites.size(); ++isite) {
    const int ixy = plot_site_ixy(sites[isite], s.config.input.nx, s.config.input.ny);
    for (std::size_t izr = 0; izr < n_zr; ++izr)
      for (std::size_t itp = 0; itp < n_tp; ++itp)
        differing_sites += (w.v_sites[itp + n_tp * (izr + n_zr * isite)] == s.reduced(itp, ixy, izr)) ? 0 : 1;
  }
  EXPECT_EQ(differing_sites, 0);
}

// At each site, v_min is the minimum of its v_sites and v_sites at the stored orientation equals v_min (bitwise).
TEST(PotentialFiles, SitesAndMinimumAgree) {
  const PotentialSetup& s = symmetric();
  const Written& w = written();
  const std::size_t n_xy = s.reduced.symmetry.n_xy();
  const std::size_t n_zr = s.reduced.index.n2;
  const std::size_t n_tp = s.reduced.index.n0;
  const int n_phi = s.angular.n_phi();
  const auto& sites = s.config.input.plot_sites;
  int disagreeing = 0;
  for (std::size_t isite = 0; isite < sites.size(); ++isite) {
    const std::size_t ixy = plot_site_ixy(sites[isite], s.config.input.nx, s.config.input.ny);
    for (std::size_t izr = 0; izr < n_zr; ++izr) {
      const auto row = w.v_sites.begin() + static_cast<std::ptrdiff_t>(n_tp * (izr + n_zr * isite));
      const std::size_t k = ixy + n_xy * izr;
      const double lowest = *std::min_element(row, row + static_cast<std::ptrdiff_t>(n_tp));
      const std::size_t itp = w.v_min_node[2 * k + 1] + n_phi * w.v_min_node[2 * k];
      disagreeing += (lowest == w.v_min[k] && row[static_cast<std::ptrdiff_t>(itp)] == w.v_min[k]) ? 0 : 1;
    }
  }
  EXPECT_EQ(disagreeing, 0);
}

// The symmetric case's v_min equals the minimum of the potential computed at every (x, y) point: bitwise at the
// representatives, to rounding elsewhere (1e-12 of max(|V|, 1 hartree)).
TEST(PotentialFiles, SymmetricCaseMatchesTheFullGrid) {
  const PotentialSetup& s = symmetric();
  const Written& w = written();
  const std::size_t n_xy = s.full.symmetry.n_xy();
  const std::size_t n_zr = s.full.index.n2;
  int representatives_differing = 0;
  double worst = 0.0;
  for (std::size_t izr = 0; izr < n_zr; ++izr)
    for (std::size_t ixy = 0; ixy < n_xy; ++ixy) {
      const double expected = minimum_by_loop(s.full, ixy, izr);
      const double got = w.v_min[ixy + n_xy * izr];
      const int irep = s.reduced.symmetry.representative_of(static_cast<int>(ixy));
      if (s.reduced.symmetry.representative_xy(irep) == static_cast<int>(ixy))
        representatives_differing += (got == expected) ? 0 : 1;
      worst = spodvr::testing::nan_aware_max(worst, scaled_difference(got, expected));
    }
  EXPECT_EQ(representatives_differing, 0);
  EXPECT_LE(worst, 1e-12);
}

// The sites' V at every orientation equals the potential computed at their (x, y) point (the bridge and the third
// site are not representatives: their values come through the symmetry maps). Negative test: the gather without
// the phi image of the symmetry operation fails exactly at the orientations where that image changes V.
TEST(PotentialFiles, SitesMatchTheFullGrid) {
  const PotentialSetup& s = symmetric();
  const Written& w = written();
  const std::size_t n_zr = s.full.index.n2;
  const std::size_t n_tp = s.full.index.n0;
  const int n_phi = s.angular.n_phi();
  const auto& sites = s.config.input.plot_sites;
  int non_representatives = 0;
  double worst = 0.0;
  int slip_fails = 0;
  int image_changes = 0;
  for (std::size_t isite = 0; isite < sites.size(); ++isite) {
    const int ixy = plot_site_ixy(sites[isite], s.config.input.nx, s.config.input.ny);
    const int irep = s.reduced.symmetry.representative_of(ixy);
    non_representatives += (s.reduced.symmetry.representative_xy(irep) == ixy) ? 0 : 1;
    const std::vector<int>& image = s.reduced.symmetry.phi_index(s.reduced.symmetry.operation_of(ixy));
    for (std::size_t izr = 0; izr < n_zr; ++izr)
      for (std::size_t itp = 0; itp < n_tp; ++itp) {
        const double full = s.full(itp, static_cast<std::size_t>(ixy), izr);
        worst = spodvr::testing::nan_aware_max(worst,
                                               scaled_difference(w.v_sites[itp + n_tp * (izr + n_zr * isite)], full));
        const double slip = gather_without_phi_map(s.reduced, itp, static_cast<std::size_t>(ixy), izr);
        slip_fails += (scaled_difference(slip, full) <= 1e-12) ? 0 : 1;
        const std::size_t itheta = itp / n_phi;
        const std::size_t imaged = image[itp % n_phi] + n_phi * itheta;
        const double stored = s.reduced.values[s.reduced.index(imaged, irep, izr)];
        const double unimaged = s.reduced.values[s.reduced.index(itp, irep, izr)];
        image_changes += (scaled_difference(unimaged, stored) <= 1e-12) ? 0 : 1;
      }
  }
  EXPECT_EQ(non_representatives, 2);
  EXPECT_LE(worst, 1e-12);
  EXPECT_GT(image_changes, 0);
  EXPECT_EQ(slip_fails, image_changes);
}

// A directory named like the plot data but not written by this program is refused before the run (it would be
// removed); an earlier plot-data directory of the same kind is replaced.
TEST(PotentialFiles, RefusesToReplaceADirectoryItDidNotWrite) {
  TemporaryDirectory dir("potential_replace");
  const std::vector<PlotDataSize> sizes = {{"potential/", "writePotential", 10, ""},
                                           {"densities/", "writeDensities", 10, ""}};
  EXPECT_NO_THROW(check_replaceable(dir.path(), sizes));
  fs::create_directories(dir.path() / "potential");
  std::ofstream(dir.path() / "potential" / "my_pes.dat") << "1 2 3\n";
  try {
    check_replaceable(dir.path(), sizes);
    FAIL() << "no InputError";
  } catch (const InputError& e) {
    EXPECT_NE(std::string(e.what()).find("writePotential"), std::string::npos) << e.what();
  }
  std::ofstream(dir.path() / "potential" / "manifest.json") << R"({"kind": "potential"})";
  EXPECT_NO_THROW(check_replaceable(dir.path(), sizes));
  fs::create_directories(dir.path() / "densities");
  std::ofstream(dir.path() / "densities" / "manifest.json") << R"({"kind": "potential"})";
  EXPECT_THROW(check_replaceable(dir.path(), sizes), InputError);
}

// A NaN at one stored value makes the minimum NaN at every (x, y) point that reads it (all points of its
// representative), and the orientation is where the NaN sits.
TEST(PotentialFiles, NaNPropagatesIntoTheMinimum) {
  const PotentialSetup& s = symmetric();
  PotentialGrid planted = s.reduced;
  const std::size_t izr = 7;
  const int irep = 2;
  planted.values[planted.index(5, static_cast<std::size_t>(irep), izr)] = std::nan("");
  const std::size_t n_xy = planted.symmetry.n_xy();
  std::vector<double> minimum(n_xy);
  std::vector<int> orientation(n_xy);
  minimum_over_orientations(planted, izr, minimum, orientation);
  int nan_points = 0;
  int members = 0;
  for (std::size_t ixy = 0; ixy < n_xy; ++ixy) {
    const bool member = planted.symmetry.representative_of(static_cast<int>(ixy)) == irep;
    members += member ? 1 : 0;
    if (std::isnan(minimum[ixy])) {
      ++nan_points;
      EXPECT_TRUE(member) << ixy;
      EXPECT_TRUE(std::isnan(planted(static_cast<std::size_t>(orientation[ixy]), ixy, izr))) << ixy;
    }
  }
  EXPECT_GT(members, 1);
  EXPECT_EQ(nan_points, members);
}

// A potential-only run of the hard wall with writePotential and a site writes potential/ (the sizes of the
// manifest's datasets), prints its size, and stops.
TEST(PotentialFiles, PotentialOnlyRunWritesThemAndStops) {
  TemporaryDirectory dir("potential_run");
  Json json = case_json("hardwall_small");
  json["Potential"] = true;
  json["writePotential"] = true;
  json["PlotSites"] = Json::array({{{"name", "origin"}, {"x", 0.0}, {"y", 0.0}}});
  std::ofstream(dir.path() / "input.json") << json.dump(1);
  CommandLine cl;
  cl.input_file = dir.path() / "input.json";
  cl.output_dir = dir.path() / "out";
  std::ostringstream out;
  std::ostringstream err;
  ASSERT_EQ(run_program(cl, out, err, false), 0) << err.str();
  EXPECT_NE(out.str().find("PLOT DATA"), std::string::npos);
  EXPECT_FALSE(fs::exists(cl.output_dir / "averages.dat"));
  const fs::path potential = cl.output_dir / "potential";
  std::ifstream manifest_file(potential / "manifest.json");
  const Json manifest = Json::parse(manifest_file);
  EXPECT_FALSE(manifest["time_dependent"].get<bool>());
  std::uintmax_t total = 0;
  for (const auto& [name, spec] : manifest["datasets"].items()) {
    std::uintmax_t bytes = spec["dtype"] == "<i4" ? 4 : 8;
    for (const int extent : spec["shape"]) bytes *= static_cast<std::uintmax_t>(extent);
    EXPECT_EQ(fs::file_size(potential / spec["file"].get<std::string>()), bytes) << name;
    total += bytes;
  }
  const Config config = build_config(parse_input(json));
  const SpatialGrids grids = make_spatial_grids(config);
  const AngularGrid angular = make_angular_grid(config);
  const std::vector<PlotDataSize> sizes = plot_data_sizes(config, grids, angular);
  ASSERT_EQ(sizes.size(), 1u);
  EXPECT_EQ(sizes[0].bytes, total);
  EXPECT_EQ(sizes[0].key, "writePotential");
  std::ifstream run_info_file(cl.output_dir / "run_info.json");
  const Json run_info = Json::parse(run_info_file);
  EXPECT_EQ(run_info["plot_data"]["sites"][0]["name"], "origin");
}

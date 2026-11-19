#pragma once

// A propagation backend that counts the calls of the time loop (rule "time-loop additions are switchable
// and free when off", docs/CONTRIBUTING.md "Additions to the time loop"): it owns the real backend,
// forwards every call and counts it in a test-owned BackendCallCounts. The tests run the program through run_program
// with a factory that wraps the real one, so the counts come from the program's own key handling and wiring.

#include <filesystem>
#include <fstream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>

#include <nlohmann/json.hpp>

#include "app/cli.hpp"
#include "app/run.hpp"
#include "propagation/propagation_backend.hpp"

namespace spodvr::test_support {

// The calls of one program run.
struct BackendCallCounts {
  int propagate_analysis_step = 0;
  int z_plane = 0;
  int block_densities = 0;
  int specular = 0;  // for information: the specular wave is held on the host anyway
  int wave_packet_on_host = 0;
  int read_wave_packet = 0;
};

class CountingBackend : public PropagationBackend {
 public:
  CountingBackend(std::unique_ptr<PropagationBackend> inner, BackendCallCounts& counts)
      : inner_(std::move(inner)), counts_(counts) {}

  void read_wave_packet(const std::function<void(WavePacket&)>& read) override {
    ++counts_.read_wave_packet;
    inner_->read_wave_packet(read);
  }
  [[nodiscard]] const WavePacket& wave_packet_on_host() override {
    ++counts_.wave_packet_on_host;
    return inner_->wave_packet_on_host();
  }
  void propagate_analysis_step() override {
    ++counts_.propagate_analysis_step;
    inner_->propagate_analysis_step();
  }
  [[nodiscard]] const SpecularWave& specular() const override {
    ++counts_.specular;
    return inner_->specular();
  }
  [[nodiscard]] std::span<const std::complex<double>> z_plane(int iz) const override {
    ++counts_.z_plane;
    return inner_->z_plane(iz);
  }
  [[nodiscard]] std::vector<double> block_densities() const override {
    ++counts_.block_densities;
    return inner_->block_densities();
  }

 private:
  std::unique_ptr<PropagationBackend> inner_;
  BackendCallCounts& counts_;
};

// A factory that wraps the backends `inner` makes in a CountingBackend counting into `counts`.
inline BackendFactory counting_factory(BackendFactory inner, BackendCallCounts& counts) {
  return
      [inner = std::move(inner), &counts](const BackendInput& input, PotentialGrid potential, SpecularWave specular) {
        return std::make_unique<CountingBackend>(inner(input, std::move(potential), std::move(specular)), counts);
      };
}

// The hard-wall input of the call tests: tests/reference/hardwall_small on a 2 x 2 (x, y) grid with
// JMax = MJMax = 2, tEnd 1000 a.u.: 100 analysis steps (dtAna 10, dtStep 2) and 11 averages rows (dtAve 100).
inline nlohmann::json small_hardwall_input() {
  std::ifstream file(std::string(SPODVR_REFERENCE_DIR) + "/hardwall_small/input.json");
  nlohmann::json json = nlohmann::json::parse(file);
  json["nX"] = 2;
  json["nY"] = 2;
  json["JMax"] = 2;
  json["MJMax"] = 2;
  json["ImposeVMax"] = false;
  json["tEnd"] = 1000.0;
  json["nStep"] = 500;
  return json;
}

constexpr int small_hardwall_analysis_steps = 100;
constexpr int small_hardwall_averages_rows = 11;  // tStart and every dtAve

// Runs the program on `json` into `dir` with the backends `inner` makes, counted; returns the exit code.
inline int run_counted(const nlohmann::json& json, const std::filesystem::path& dir, BackendFactory inner,
                       BackendCallCounts& counts) {
  std::filesystem::create_directories(dir);
  std::ofstream(dir / "input.json") << json.dump(1);
  CommandLine cl;
  cl.input_file = dir / "input.json";
  cl.output_dir = dir;
  std::ostringstream out;
  std::ostringstream err;
  const int code = run_program(cl, out, err, false, {}, counting_factory(std::move(inner), counts));
  if (code != 0) ADD_FAILURE() << err.str();
  return code;
}

}  // namespace spodvr::test_support

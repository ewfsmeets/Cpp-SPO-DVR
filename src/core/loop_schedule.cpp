// The OpenMP schedule of the (Z, r) group loops, chosen once at start-up (see loop_schedule.hpp).
#include "core/loop_schedule.hpp"

#include <omp.h>

#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <system_error>

namespace spodvr {

int numa_node_count() {
  std::error_code error;
  std::filesystem::directory_iterator nodes("/sys/devices/system/node", error);
  if (error) return 1;
  int count = 0;
  for (const auto& entry : nodes) {
    const std::string name = entry.path().filename().string();
    if (name.starts_with("node") && name.size() > 4 && std::isdigit(static_cast<unsigned char>(name[4]))) ++count;
  }
  return count > 0 ? count : 1;
}

LoopSchedule default_loop_schedule(int numa_nodes) {
  return numa_nodes > 1 ? LoopSchedule::RoundRobin : LoopSchedule::Dynamic;
}

LoopSchedule parse_loop_schedule(const std::string& text) {
  if (text == "dynamic") return LoopSchedule::Dynamic;
  if (text == "static") return LoopSchedule::RoundRobin;
  throw std::invalid_argument("--schedule must be 'dynamic' or 'static' (got '" + text + "')");
}

std::string to_string(LoopSchedule schedule) { return schedule == LoopSchedule::Dynamic ? "dynamic,1" : "static,1"; }

void apply_loop_schedule(LoopSchedule schedule) {
  omp_set_schedule(schedule == LoopSchedule::Dynamic ? omp_sched_dynamic : omp_sched_static, 1);
}

}  // namespace spodvr

#pragma once

#include <string>

namespace spodvr {

// The OpenMP schedule of the loops over (Z, r) work groups: every loop that computes on, or first
// writes, the large per-(Z, r) arrays uses `schedule(runtime)` and iterates over the same groups.
// It is set once at start-up (apply_loop_schedule); results do not depend on it.
// - Dynamic (dynamic,1): a thread takes the next group when it is done. Balances the uneven cost
//   per group (below Z_asymptotic: angular transform + exp(-iV dt); beyond: one multiply). For
//   machines with one NUMA node.
// - RoundRobin (static,1): thread t always gets groups t, t+T, t+2T, ... Balanced to within one
//   group, and every thread works on the memory it first wrote, so that memory stays on its own
//   socket. For machines with several NUMA nodes.
enum class LoopSchedule { Dynamic, RoundRobin };

// Number of NUMA nodes (/sys/devices/system/node/node*); 1 if the system does not report them.
[[nodiscard]] int numa_node_count();

// RoundRobin with several NUMA nodes, Dynamic otherwise.
[[nodiscard]] LoopSchedule default_loop_schedule(int numa_nodes);

// "dynamic" or "static" (both with chunk 1). Throws std::invalid_argument otherwise.
[[nodiscard]] LoopSchedule parse_loop_schedule(const std::string& text);

// "dynamic,1" or "static,1".
[[nodiscard]] std::string to_string(LoopSchedule schedule);

// Sets the schedule of all `schedule(runtime)` loops (omp_set_schedule; overrides OMP_SCHEDULE).
void apply_loop_schedule(LoopSchedule schedule);

}  // namespace spodvr

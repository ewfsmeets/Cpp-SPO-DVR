// Tests: the OpenMP schedule of the (Z, r) group loops (core/loop_schedule).
#include <gtest/gtest.h>
#include <omp.h>

#include <stdexcept>

#include "core/loop_schedule.hpp"

using namespace spodvr;

TEST(LoopSchedule, ParsesTheTwoSchedules) {
  EXPECT_EQ(parse_loop_schedule("dynamic"), LoopSchedule::Dynamic);
  EXPECT_EQ(parse_loop_schedule("static"), LoopSchedule::RoundRobin);
  EXPECT_THROW((void)parse_loop_schedule("guided"), std::invalid_argument);
  EXPECT_THROW((void)parse_loop_schedule("dynamic,4"), std::invalid_argument);
  EXPECT_EQ(to_string(LoopSchedule::Dynamic), "dynamic,1");
  EXPECT_EQ(to_string(LoopSchedule::RoundRobin), "static,1");
}

TEST(LoopSchedule, DefaultFollowsTheNumberOfNumaNodes) {
  EXPECT_EQ(default_loop_schedule(1), LoopSchedule::Dynamic);
  EXPECT_EQ(default_loop_schedule(2), LoopSchedule::RoundRobin);
  EXPECT_GE(numa_node_count(), 1);
}

TEST(LoopSchedule, SetsTheRuntimeScheduleWithChunkOne) {
  omp_sched_t kind{};
  int chunk = 0;
  apply_loop_schedule(LoopSchedule::RoundRobin);
  omp_get_schedule(&kind, &chunk);
  EXPECT_EQ(kind, omp_sched_static);
  EXPECT_EQ(chunk, 1);
  apply_loop_schedule(LoopSchedule::Dynamic);
  omp_get_schedule(&kind, &chunk);
  EXPECT_EQ(kind, omp_sched_dynamic);
  EXPECT_EQ(chunk, 1);
}

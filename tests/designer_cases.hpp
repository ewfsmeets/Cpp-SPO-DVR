#pragma once
// The shared cases of the wave-packet designer's drift guard (tools/wave_packet_designer/expected_cases.json):
// the designer's copies of the grid indices, angular sizes and memory estimates against the program's. When a
// check here fails, the program changed something the designer copies: update the copy in the designer's
// program/ and regenerate the file (tools/wave_packet_designer/README.md, "Keeping up with the C++").

#include <gtest/gtest.h>

#include <cstddef>
#include <fstream>
#include <optional>
#include <string>

#include <nlohmann/json.hpp>

#include "core/index_range.hpp"
#include "propagation/memory_estimate.hpp"

namespace spodvr::testing {

// A file of the repository (SPODVR_SOURCE_DIR: its root, compiled in).
inline nlohmann::json read_repository_json(const std::string& relative_path) {
  std::ifstream file(std::string(SPODVR_SOURCE_DIR) + "/" + relative_path);
  return nlohmann::json::parse(file);
}

// The cases of the file.
inline nlohmann::json designer_cases() {
  return read_repository_json("tools/wave_packet_designer/expected_cases.json").at("cases");
}

// The input of a case: its base input file with the case's key changes (null removes the key).
inline nlohmann::json designer_case_input(const nlohmann::json& designer_case) {
  nlohmann::json input = read_repository_json(designer_case.at("input").get<std::string>());
  for (const auto& [key, value] : designer_case.at("changes").items()) {
    if (value.is_null())
      input.erase(key);
    else
      input[key] = value;
  }
  return input;
}

// An index range against the file's [first, last], or null for an absorber that is off.
inline void expect_designer_range(const std::optional<IndexRange>& range, const nlohmann::json& expected,
                                  const std::string& name) {
  SCOPED_TRACE(name);
  ASSERT_EQ(range.has_value(), !expected.is_null());
  if (!range) return;
  EXPECT_EQ(range->first, expected.at(0).get<int>());
  EXPECT_EQ(range->last, expected.at(1).get<int>());
}

// A memory estimate against the file's: every item (name, values, bytes, when allocated) and the totals.
inline void expect_designer_estimate(const MemoryEstimate& estimate, const nlohmann::json& expected) {
  const nlohmann::json& items = expected.at("items");
  ASSERT_EQ(estimate.items.size(), items.size());
  for (std::size_t i = 0; i < items.size(); ++i) {
    const MemoryItem& item = estimate.items[i];
    SCOPED_TRACE("memory item " + item.name);
    EXPECT_EQ(item.name, items[i].at("name").get<std::string>());
    EXPECT_EQ(item.values, items[i].at("values").get<std::size_t>());
    EXPECT_EQ(item.bytes, items[i].at("bytes").get<std::size_t>());
    EXPECT_EQ(item.during_setup, items[i].at("during_setup").get<bool>());
    EXPECT_EQ(item.during_propagation, items[i].at("during_propagation").get<bool>());
  }
  EXPECT_EQ(estimate.setup_bytes, expected.at("setup_bytes").get<std::size_t>());
  EXPECT_EQ(estimate.propagation_bytes, expected.at("propagation_bytes").get<std::size_t>());
  EXPECT_EQ(estimate.peak_bytes, expected.at("peak_bytes").get<std::size_t>());
}

}  // namespace spodvr::testing

// Unit tests of the screen formatting (report/format, report/propagation_report): durations and
// the live progress line, which appears only with colour so that logs stay plain text.
#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "report/format.hpp"
#include "report/propagation_report.hpp"

namespace spodvr {
namespace {

bool has_control_codes(const std::string& text) {
  return text.find('\r') != std::string::npos || text.find('\033') != std::string::npos;
}

TEST(Format, DurationUnderAnHourIsMinutesAndSeconds) {
  EXPECT_EQ(format_duration(0.0), "0:00");
  EXPECT_EQ(format_duration(59.4), "0:59");
  EXPECT_EQ(format_duration(61.0), "1:01");
}

TEST(Format, DurationFromAnHourHasHours) {
  EXPECT_EQ(format_duration(3600.0), "1:00:00");
  EXPECT_EQ(format_duration(5532.1), "1:32:12");
}

TEST(AveragesTable, WithoutColourThereIsNoProgressLineAndNoControlCodes) {
  std::ostringstream out;
  const Style style{false};
  AveragesTable table(out, style, 8, 2.5, 0.0);
  table.print_header();
  table.print_row(Averages{0.0, 1.0, 17.0, 1.4, 1.0});
  table.print_status(10.0);
  table.print_note("analysis surface fixed", 10.0);
  table.finish();
  const std::string text = out.str();
  EXPECT_FALSE(has_control_codes(text));
  EXPECT_EQ(text.find("step 4/8"), std::string::npos);
  EXPECT_NE(text.find("Propagation time: "), std::string::npos);
}

TEST(AveragesTable, WithColourTheProgressLineIsDrawnAndClearedBeforeTheNextRow) {
  std::ostringstream out;
  const Style style{true};
  AveragesTable table(out, style, 8, 2.5, 0.0);
  table.print_status(10.0);
  const std::string status = out.str();
  EXPECT_EQ(status.rfind("\r\033[K", 0), 0u);
  EXPECT_NE(status.find("step 4/8"), std::string::npos);
  EXPECT_NE(status.find(" 50%"), std::string::npos);
  table.print_row(Averages{10.0, 1.0, 17.0, 1.4, 1.0});
  EXPECT_EQ(out.str().substr(status.size(), 4), "\r\033[K");
}

}  // namespace
}  // namespace spodvr

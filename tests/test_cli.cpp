#include <gtest/gtest.h>

#include <stdexcept>

#include "app/cli.hpp"

using namespace spodvr;

TEST(Cli, ParsesAllOptions) {
  const auto cl = parse_command_line({"in.json", "--threads", "4", "--output-dir", "out", "--no-color"});
  EXPECT_EQ(cl.input_file, "in.json");
  EXPECT_EQ(cl.threads, 4);
  EXPECT_EQ(cl.output_dir, "out");
  EXPECT_TRUE(cl.no_color);
}

TEST(Cli, HelpAndVersionNeedNoInput) {
  EXPECT_TRUE(parse_command_line({"--help"}).help);
  EXPECT_TRUE(parse_command_line({"--version"}).version);
}

TEST(Cli, Errors) {
  EXPECT_THROW((void)parse_command_line({}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "b.json"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--threads"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--threads", "0"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--threads", "4x"}), std::invalid_argument);
  EXPECT_THROW((void)parse_command_line({"a.json", "--bogus"}), std::invalid_argument);
}

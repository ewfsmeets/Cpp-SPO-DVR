// Program entry point: parses the command line and runs the calculation.
#include <unistd.h>

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "app/cli.hpp"
#include "app/run.hpp"

int main(int argc, char** argv) {
  const std::vector<std::string> args(argv + 1, argv + argc);
  spodvr::CommandLine command_line;
  try {
    command_line = spodvr::parse_command_line(args);
  } catch (const std::invalid_argument& e) {
    std::cerr << "spodvr: " << e.what() << "\n\n" << spodvr::usage();
    return 2;
  }
  if (command_line.help) {
    std::cout << spodvr::usage();
    return 0;
  }
  if (command_line.version) {
    std::cout << "spodvr " SPODVR_VERSION "\n";
    return 0;
  }
  const bool color = !command_line.no_color && isatty(STDOUT_FILENO) != 0;
  return spodvr::run_program(command_line, std::cout, std::cerr, color);
}

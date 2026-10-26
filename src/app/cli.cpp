#include "app/cli.hpp"

#include <stdexcept>

namespace spodvr {

std::string usage() {
  return "Usage: spodvr <input.json> [--threads N] [--output-dir DIR] [--no-color] [--help] [--version]\n"
         "\n"
         "  --threads N       number of OpenMP threads (default: OpenMP default)\n"
         "  --output-dir DIR  directory for the output data files (default: current directory)\n"
         "  --no-color        disable ANSI colours (also disabled when stdout is not a terminal)\n"
         "  --help            show this help\n"
         "  --version         show the program version\n";
}

CommandLine parse_command_line(const std::vector<std::string>& args) {
  CommandLine cl;
  bool have_input = false;
  for (std::size_t i = 0; i < args.size(); ++i) {
    const std::string& a = args[i];
    auto value = [&](const std::string& option) -> const std::string& {
      if (i + 1 >= args.size()) throw std::invalid_argument(option + " requires a value");
      return args[++i];
    };
    if (a == "--help" || a == "-h") {
      cl.help = true;
    } else if (a == "--version") {
      cl.version = true;
    } else if (a == "--no-color") {
      cl.no_color = true;
    } else if (a == "--threads") {
      const std::string& v = value(a);
      std::size_t used = 0;
      int n = 0;
      try {
        n = std::stoi(v, &used);
      } catch (const std::exception&) {
        used = 0;
      }
      if (used != v.size() || n < 1) throw std::invalid_argument("--threads requires a positive integer");
      cl.threads = n;
    } else if (a == "--output-dir") {
      cl.output_dir = value(a);
    } else if (!a.empty() && a[0] == '-') {
      throw std::invalid_argument("unknown option '" + a + "'");
    } else {
      if (have_input) throw std::invalid_argument("more than one input file given");
      cl.input_file = a;
      have_input = true;
    }
  }
  if (!have_input && !cl.help && !cl.version) throw std::invalid_argument("no input file given");
  return cl;
}

}  // namespace spodvr

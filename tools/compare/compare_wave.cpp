// Compares two wave packets written by spodvr_dump_wave: bitwise identical, and the largest
// difference relative to the largest amplitude. Exit code 0 if the difference is within the
// tolerance (default 0: bitwise identical).
// Usage: spodvr_compare_wave <a> <b> [tolerance]
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace {

std::vector<double> read_values(const char* path) {
  std::ifstream in(path, std::ios::binary);
  std::vector<double> values;
  double x = 0.0;
  while (in.read(reinterpret_cast<char*>(&x), sizeof(x))) values.push_back(x);
  return values;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <a> <b> [tolerance]\n", argv[0]);
    return 2;
  }
  const double tolerance = argc > 3 ? std::stod(argv[3]) : 0.0;
  const std::vector<double> a = read_values(argv[1]);
  const std::vector<double> b = read_values(argv[2]);
  if (a.empty() || a.size() != b.size()) {
    std::printf("DIFFERENT: sizes %zu and %zu\n", a.size(), b.size());
    return 1;
  }

  bool identical = true;
  double max_difference = 0.0;
  double max_amplitude = 0.0;
  for (std::size_t i = 0; i + 1 < a.size(); i += 2) {
    const double dre = a[i] - b[i];
    const double dim = a[i + 1] - b[i + 1];
    max_difference = std::max(max_difference, std::hypot(dre, dim));
    max_amplitude = std::max(max_amplitude, std::hypot(a[i], a[i + 1]));
    if (a[i] != b[i] || a[i + 1] != b[i + 1]) identical = false;
  }
  const double relative = max_amplitude > 0.0 ? max_difference / max_amplitude : max_difference;
  std::printf("%s: max |a - b| / max |a| = %.3e\n", identical ? "bitwise identical" : "different", relative);
  return relative <= tolerance ? 0 : 1;
}

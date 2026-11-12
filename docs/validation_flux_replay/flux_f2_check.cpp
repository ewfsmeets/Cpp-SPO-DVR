// Flux check 2 (docs/validation_flux_replay/README.md): replays the 5 analysis steps of
// flux_f2_check.f90 (stdin: hex bits) through the C++ FluxAccumulator and flux_results; phi and phi'
// must be bitwise equal to the Fortran's flux coefficient files, P and P(x, y) within half a printed
// unit of its E15.8 files. Usage: flux_f2_check <Fortran run dir> <tests/reference> < <Fortran output>
// Exit status: 0 all checks pass, 1 otherwise. Built and run by run_checks.sh.
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <nlohmann/json.hpp>

#include "analysis/flux_analysis.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"

using namespace spodvr;
using cplx = std::complex<double>;

double bits(const std::string& h) { return std::bit_cast<double>(std::stoull(h, nullptr, 16)); }
bool same(double a, const std::string& h) { return std::bit_cast<std::uint64_t>(a) == std::stoull(h, nullptr, 16); }

// |value - printed| <= half a unit of the last printed digit of an E15.8 number.
bool within_half_unit(double value, const std::string& printed) {
  const auto e = printed.find('E');
  const int exponent = std::stoi(printed.substr(e + 1));
  return std::abs(value - std::stod(printed)) <= 0.5 * std::pow(10.0, exponent - 8) * (1.0 + 1e-9);
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: flux_f2_check <Fortran run dir> <tests/reference> < <Fortran output>\n");
    return 2;
  }
  const std::string run_dir = argv[1];
  const std::string reference_dir = argv[2];
  const int n_rot = 2, nx = 1, ny = 2, nr = 8, n_planes = 2, ir_flux = 2, steps = 5;
  std::map<std::tuple<int, int, int, int, int>, cplx> psi;  // (step, irot, ixy, ir, iz), Fortran 1-based except step
  std::vector<double> times(steps);
  std::map<std::tuple<char, int, int, int>, std::pair<std::string, std::string>> fortran;  // (F|D, ien, iz, iredist)
  std::string g_mass, g_reduced, g_skew, g_dt;
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string tag;
    in >> tag;
    if (tag == "G") {
      in >> g_mass >> g_reduced >> g_skew >> g_dt;
    } else if (tag == "P") {
      int s, irot, ixy, ir, iz;
      std::string re, im;
      in >> s >> irot >> ixy >> ir >> iz >> re >> im;
      psi[{s, irot, ixy, ir, iz}] = cplx(bits(re), bits(im));
    } else if (tag == "T") {
      int s;
      std::string t;
      in >> s >> t;
      times[s] = bits(t);
    } else if (tag == "F" || tag == "D") {
      int ien, iz, iredist;
      std::string re, im;
      in >> ien >> iz >> iredist >> re >> im;
      fortran[{tag[0], ien, iz, iredist}] = {re, im};
    }
  }

  std::ifstream json_file(reference_dir + "/hardwall_small/input.json");
  nlohmann::json json = nlohmann::json::parse(json_file);
  json["doFluxAnalysis"] = true;
  json["rFluxAnalysis"] = 2.0;
  json["EMinFluxAna"] = 0.2;
  json["EMaxFluxAna"] = 0.6;
  json["DeltaEFluxAna"] = 0.15;
  const Config config = build_config(parse_input(json));
  std::printf("total mass bitwise the Fortran's rmass: %s\n", same(config.derived.total_mass, g_mass) ? "yes" : "NO");
  const FluxEnergies energies = make_flux_energies(config, 0.0123, 0.0098765);

  const std::vector<double> weights = flux_derivative_weights(nr, 0.15, ir_flux);
  FluxAccumulator flux(n_rot * nx * ny, n_planes, energies, bits(g_dt), 2);
  for (int s = 0; s < steps; ++s)
    flux.add_step(times[s], [&](FluxSurfaceWave& slot) {
      for (int iz = 0; iz < n_planes; ++iz)
        for (int ixy = 0; ixy < nx * ny; ++ixy)
          for (int irot = 0; irot < n_rot; ++irot) {
            std::vector<cplx> column(nr);
            for (int ir = 0; ir < nr; ++ir) column[ir] = psi.at({s, irot + 1, ixy + 1, ir + 1, iz + 1});
            const int block = irot + n_rot * ixy;
            slot.psi[slot.index(block, iz)] = column[ir_flux];
            slot.derivative[slot.index(block, iz)] = flux_r_derivative(weights, column);
          }
    });
  flux.flush();
  const FluxCoefficients& c = flux.coefficients();

  int total = 0, equal = 0;
  for (const auto& [key, value] : fortran) {
    const auto [tag, ien, iz, iredist] = key;
    const std::size_t k = c.index(iredist - 1, iz - 1, ien - 1);
    const cplx v = (tag == 'F') ? c.phi[k] : c.phi_derivative[k];
    ++total;
    if (same(v.real(), value.first) && same(v.imag(), value.second))
      ++equal;
    else
      std::printf("%c ien %d iz %d iredist %d differs: %.17g %.17g vs %.17g %.17g\n", tag, ien, iz, iredist, v.real(),
                  v.imag(), bits(value.first), bits(value.second));
  }
  std::printf("phi and phi': %d of %d values bitwise equal\n", equal, total);

  const FluxGrid grid{.n_rot = n_rot,
                      .nx = nx,
                      .ny = ny,
                      .dz = 0.1,
                      .dx = 0.8,
                      .dy = 0.9,
                      .skewing_angle = bits(g_skew),
                      .total_mass = bits(g_mass),
                      .reduced_mass = bits(g_reduced)};
  const FluxResults r = flux_results(c, energies, grid);
  std::ifstream reaction(run_dir + "/flux_reaction.dat");
  int p_ok = 0, p_total = 0;
  for (int n = 0; n < energies.size(); ++n) {
    std::string e_text, p_text;
    reaction >> e_text >> p_text;
    ++p_total;
    if (within_half_unit(r.probability[n], p_text) && within_half_unit(r.energy[n], e_text)) ++p_ok;
    std::printf("E %s P %s   C++ %.16e %.16e\n", e_text.c_str(), p_text.c_str(), r.energy[n], r.probability[n]);
  }
  std::ifstream density(run_dir + "/flux_density.dat");
  std::string skip;
  for (int i = 0; i < 4; ++i) std::getline(density, skip);
  int d_ok = 0, d_total = 0;
  for (int n = 0; n < energies.size(); ++n)
    for (int ix = 0; ix < nx; ++ix)
      for (int iy = 0; iy < ny; ++iy) {
        std::string e_text, x_text, y_text, d_text;
        density >> e_text >> x_text >> y_text >> d_text;
        ++d_total;
        if (within_half_unit(r.density[r.density_index(ix, iy, n)], d_text)) ++d_ok;
      }
  std::printf("P(E): %d of %d and P(x, y; E): %d of %d within half a printed unit\n", p_ok, p_total, d_ok, d_total);
  return (equal == total && total == 48 && p_ok == p_total && d_ok == d_total) ? 0 : 1;
}

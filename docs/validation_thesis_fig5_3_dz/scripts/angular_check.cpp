// Validation helper (angular matrices, docs/validation_thesis_fig5_3_dz): the angular DVR<->FBR
// transform of a run, built with the library of the build that ran, and its round trip.
// Usage: angular_check <input.json> <output dir>
// Writes into the output directory (raw little-endian float64 / complex128, no headers):
//   thetas.bin, weights.bin  the theta nodes and weights (AngularGrid)
//   matrix.bin               AngularTransform::matrix(): sqrt(w_t / nphi) P_J^m(theta_t), ordered
//                            (mJ block, theta, J)
//   basis.txt                n_rot, nphi, ntheta; then per FBR index: J mJ
//   roundtrip.bin            for every mJ = 0 basis function e_J (J ascending): the FBR vector
//                            dvr_to_fbr(fbr_to_dvr(e_J)) of column 0 (n_rot complex values); all
//                            nx*ny columns carry e_J, as the uniform (n, m) = (0, 0) component in
//                            the potential operator (the run's column count)
//   columns_identical.txt    1 if every column's result equals column 0 bitwise
#include <algorithm>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <span>
#include <vector>

#include "fft/fftw_fft.hpp"
#include "grids/angular_grid.hpp"
#include "grids/rotational_basis.hpp"
#include "input/json_input.hpp"
#include "input/process.hpp"
#include "io/binary_file.hpp"
#include "propagation/angular_transform.hpp"

using namespace spodvr;
namespace fs = std::filesystem;
using cplx = std::complex<double>;

namespace {

template <typename T>
void write_binary(const fs::path& path, std::span<const T> values) {
  BinaryWriter file(path);
  file.write_array(values);
  file.close();
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "usage: angular_check <input.json> <output dir>\n";
    return 2;
  }
  const fs::path out_dir = argv[2];
  fs::create_directories(out_dir);
  const Config config = build_config(read_input_file(argv[1]));
  const AngularGrid angular = make_angular_grid(config);
  const RotationalBasis basis(angular, config.input.j_in, config.input.mj_in);
  const int n_columns = config.input.nx * config.input.ny;  // as CpuBackend
  const AngularTransform transform(angular, basis, FftwFft(), n_columns);

  write_binary(out_dir / "thetas.bin", angular.thetas());
  write_binary(out_dir / "weights.bin", angular.theta_weights());
  write_binary(out_dir / "matrix.bin", std::span<const double>(transform.matrix()));
  {
    std::ofstream out(out_dir / "basis.txt");
    out << basis.size() << ' ' << angular.n_phi() << ' ' << angular.n_theta() << '\n';
    for (int i = 0; i < basis.size(); ++i) out << basis.j(i) << ' ' << basis.mj(i) << '\n';
  }

  const int n_rot = basis.size();
  std::vector<cplx> fbr(std::size_t(n_rot) * n_columns);
  std::vector<cplx> dvr(std::size_t(transform.dvr_block_size()) * n_columns);
  std::vector<cplx> roundtrip;
  bool columns_identical = true;
  const IndexRange m0 = basis.mj_block(0);
  for (int irot = m0.first; irot < m0.first + m0.size(); ++irot) {
    std::fill(fbr.begin(), fbr.end(), cplx(0.0, 0.0));
    for (int c = 0; c < n_columns; ++c) fbr[irot + std::size_t(n_rot) * c] = 1.0;
    std::fill(dvr.begin(), dvr.end(), cplx(0.0, 0.0));
    transform.fbr_to_dvr(fbr, dvr);
    transform.dvr_to_fbr(dvr, fbr);
    for (int c = 1; c < n_columns; ++c)
      for (int i = 0; i < n_rot; ++i)
        if (fbr[i + std::size_t(n_rot) * c] != fbr[i]) columns_identical = false;
    roundtrip.insert(roundtrip.end(), fbr.begin(), fbr.begin() + n_rot);
  }
  write_binary(out_dir / "roundtrip.bin", std::span<const cplx>(roundtrip));
  std::ofstream(out_dir / "columns_identical.txt") << (columns_identical ? 1 : 0) << '\n';
  std::cout << "angular_check: " << n_rot << " basis functions, " << m0.size() << " round trips, ntheta "
            << angular.n_theta() << ", nphi " << angular.n_phi() << '\n';
  return 0;
}

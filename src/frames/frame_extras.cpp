// The optional datasets of a frame (FrameExtras): the angular, diffraction, vibrational, rotational-moment and
// momentum reductions of psi_full, one Z plane at a time. C++ only.
#include "frames/frame_extras.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "core/aligned_vector.hpp"
#include "core/array_index.hpp"

namespace spodvr {

namespace {

using Complex = std::complex<double>;

enum Extra { angular_extra, diffraction_extra, vibration_extra, rotation_extra, momentum_extra, n_extras };

// <J m| cos(theta) |J+1 m> of normalised functions.
double cos_coupling(int j, int m) {
  return std::sqrt(((j + 1.0) * (j + 1.0) - m * m) / ((2.0 * j + 1.0) * (2.0 * j + 3.0)));
}

}  // namespace

// <J m| P2 |J' m> = (3 <J m|cos^2|J' m> - delta) / 2, <cos^2> from the cos couplings (J' = J, J +- 2).
double p2_element(int j, int jp, int m) {
  if (std::abs(m) > std::min(j, jp)) return 0.0;
  if (j == jp) {
    const double up = cos_coupling(j, m);
    const double down = j - 1 >= std::abs(m) ? cos_coupling(j - 1, m) : 0.0;
    return (3.0 * (up * up + down * down) - 1.0) / 2.0;
  }
  if (std::abs(j - jp) == 2) return 1.5 * cos_coupling(std::min(j, jp), m) * cos_coupling(std::min(j, jp) + 1, m);
  return 0.0;
}

namespace {

// Times one extra's work into seconds[extra].
class ExtraTimer {
 public:
  ExtraTimer(std::vector<double>& seconds, Extra extra)
      : seconds_(seconds), extra_(extra), start_(std::chrono::steady_clock::now()) {}
  ~ExtraTimer() {
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_;
    seconds_[extra_] += elapsed.count();
  }
  ExtraTimer(const ExtraTimer&) = delete;
  ExtraTimer& operator=(const ExtraTimer&) = delete;

 private:
  std::vector<double>& seconds_;
  Extra extra_;
  std::chrono::steady_clock::time_point start_;
};

}  // namespace

FrameExtraSwitches FrameExtraSwitches::of(const InputParameters& p) {
  FrameExtraSwitches on;
  for (const std::string& name : p.frame_extras) {
    on.angular = on.angular || name == "angular";
    on.diffraction = on.diffraction || name == "diffraction";
    on.vibration = on.vibration || name == "vibration";
    on.rotation_moments = on.rotation_moments || name == "rotation_moments";
    on.momentum = on.momentum || name == "momentum";
  }
  return on;
}

std::vector<int> basis_j_values(const RotationalBasis& basis) {
  std::vector<int> j;
  for (int irot = 0; irot < basis.size(); ++irot) j.push_back(basis.j(irot));
  std::ranges::sort(j);
  j.erase(std::unique(j.begin(), j.end()), j.end());
  return j;
}

FrameExtraSizes frame_extra_sizes(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                                  const RotationalBasis& basis) {
  const FrameExtraSwitches on = FrameExtraSwitches::of(config.input);
  const std::size_t nz = grids.z_scattering.size();
  const std::size_t nr = grids.r.size();
  FrameExtraSizes sizes;
  if (on.angular) sizes.angle_z = nz * angular.n_theta() * angular.n_phi();
  if (on.diffraction) sizes.nm_z = nz * config.input.nx * config.input.ny;
  if (on.vibration) sizes.vj_z = nz * nr * basis_j_values(basis).size();
  if (on.rotation_moments) sizes.moments_z = nz;
  if (on.momentum) {
    sizes.kr = nr;
    sizes.kz_spec = grids.z_specular.size();
  }
  return sizes;
}

// The tables and transforms of the extras that are on.
struct FrameExtras::Tables {
  // angular: the FBR -> DVR transform of one column pass (columns_per_pass of the nX nY columns of a (Z, r)
  // point, as the potential step)
  std::unique_ptr<AngularTransform> angular;
  // diffraction: the x-y FFTs of the n_rot columns of one (Z, r) point (layout (irot, ix, iy)), the shift
  // exp(-i k_par . R) and the probability per |FFT value|^2: dZ dr dx dy sin(skew) / (nX nY)
  std::unique_ptr<ComplexFftPlan> xy_fft;
  std::span<const Complex> shift_forward;
  double diffraction_scale = 0.0;
  // vibration: phi_vJ(r) dr for each J of the basis (index jv = v + nr j_index), the J index of each state, and the
  // probability per |c_v|^2: dZ dx dy sin(skew)
  std::vector<double> phi_dr;
  std::vector<int> j_index;
  int n_j = 0;
  double vibration_scale = 0.0;
  // rotation_moments: the non-zero P2 elements (same mJ) and J(J+1) / (2 mu r^2) per (irot, ir)
  struct Element {
    int row;
    int col;
    double value;
  };
  std::vector<Element> p2;
  std::span<const double> rotational_energy;  // the propagation's table, (irot, ir)
  // momentum: the r FFT of a whole plane (n_rot nX nY transforms of nR points) and the Z FFT of the specular wave
  // (nR transforms of nZSp points); the probability per |FFT value|^2: the cell volume / nR, dZ dr / nZSp
  std::unique_ptr<ComplexFftPlan> r_fft;
  std::unique_ptr<ComplexFftPlan> z_fft;
  std::vector<Complex> plane_work;
  std::vector<Complex> specular_work;
};

FrameExtras::FrameExtras(const Config& config, const SpatialGrids& grids, const AngularGrid& angular,
                         const RotationalBasis& basis, const GasPhase& gas, const PropagationTables& tables,
                         const Fft& fft, const FrameLayout& layout)
    : layout_(layout),
      on_(FrameExtraSwitches::of(config.input)),
      tables_(std::make_unique<Tables>()),
      seconds_(n_extras, 0.0) {
  const InputParameters& p = config.input;
  const DerivedParameters& d = config.derived;
  const double dz = grids.z_scattering.spacing();
  const double dr = grids.r.spacing();
  Tables& t = *tables_;
  if (on_.angular)
    t.angular = std::make_unique<AngularTransform>(
        angular, basis, fft, AngularTransform::columns_per_pass(layout.nxy, angular.n_theta(), angular.n_phi()));
  if (on_.diffraction) {
    const FftLayout2D xy{.outer_length = p.nx,
                         .outer_stride = layout.n_rot * p.ny,
                         .inner_length = p.ny,
                         .inner_stride = layout.n_rot,
                         .howmany = layout.n_rot,
                         .distance = 1};
    t.xy_fft = fft.plan_complex_2d(xy, FftSign::Minus, 1.0);  // the BK analysis's convention
    t.shift_forward = tables.shift_forward;
    t.diffraction_scale = dz * dr * d.dx * d.dy * std::sin(d.system.skewing_angle) / (p.nx * p.ny);
  }
  if (on_.vibration) {
    const std::vector<int> j_values = basis_j_values(basis);
    t.n_j = static_cast<int>(j_values.size());
    for (int irot = 0; irot < basis.size(); ++irot)
      t.j_index.push_back(static_cast<int>(std::ranges::find(j_values, basis.j(irot)) - j_values.begin()));
    for (const int j : j_values)
      for (int v = 0; v < layout.nr; ++v)
        for (const double phi : gas.states.wave_function(v, j)) t.phi_dr.push_back(phi * dr);
    t.vibration_scale = dz * d.xy_volume_element;
  }
  if (on_.rotation_moments) {
    for (int row = 0; row < basis.size(); ++row)
      for (int col = 0; col < basis.size(); ++col)
        if (basis.mj(row) == basis.mj(col))
          if (const double value = p2_element(basis.j(row), basis.j(col), basis.mj(row)); value != 0.0)
            t.p2.push_back({row, col, value});
    t.rotational_energy = tables.rotational_energies;
  }
  if (on_.momentum) {
    const int columns = layout.n_rot * layout.nxy;
    t.r_fft = fft.plan_complex({.length = layout.nr, .howmany = columns, .stride = columns, .distance = 1},
                               FftSign::Minus, 1.0);
    t.z_fft = fft.plan_complex({.length = layout.nz_specular, .howmany = layout.nr, .stride = layout.nr, .distance = 1},
                               FftSign::Minus, 1.0);
    t.plane_work.resize(std::size_t(columns) * layout.nr);
    t.specular_work.resize(std::size_t(layout.nz_specular) * layout.nr);
  }
}

FrameExtras::~FrameExtras() = default;

std::size_t FrameExtras::memory_values(const FrameLayout& layout, const AngularGrid& angular,
                                       const FrameExtraSwitches& on, int n_j, int n_threads) {
  const std::size_t nz = layout.nz;
  const std::size_t nr = layout.nr;
  const std::size_t nxy = layout.nxy;
  const std::size_t plane = std::size_t(layout.n_rot) * nxy * nr;
  const std::size_t nodes = std::size_t(angular.n_theta()) * angular.n_phi();
  const std::size_t threads = n_threads;
  std::size_t values = 0;
  if (on.any()) values += 2 * plane;  // the plane of psi_full (complex)
  if (on.angular)  // dataset, per-r sums, the DVR blocks of one column pass per thread
    values += nz * nodes + nr * nodes +
              threads * 2 * std::size_t(AngularTransform::dvr_block_size(angular.n_theta(), angular.n_phi())) *
                  AngularTransform::columns_per_pass(layout.nxy, angular.n_theta(), angular.n_phi());
  if (on.diffraction)  // dataset, per-r sums, a (Z, r) block per thread
    values += nz * nxy + nr * nxy + threads * 2 * std::size_t(layout.n_rot) * nxy;
  if (on.vibration)  // dataset, phi_vJ dr, per-(x, y) sums
    values += nz * nr * n_j + std::size_t(n_j) * nr * nr + nxy * nr * n_j;
  if (on.rotation_moments) values += 2 * nz + 2 * nr;  // datasets, per-r sums (the energies are the tables')
  if (on.momentum)  // datasets, the plane and the specular wave to transform
    values += nr + std::size_t(layout.nz_specular) + 2 * plane + 2 * std::size_t(layout.nz_specular) * nr;
  return values;
}

void FrameExtras::reduce_plane(int iz, std::span<const Complex> full_plane, FrameRecord& record) {
  if (on_.angular) reduce_angular(iz, full_plane, record);
  if (on_.diffraction) reduce_diffraction(iz, full_plane, record);
  if (on_.vibration) reduce_vibration(iz, full_plane, record);
  if (on_.rotation_moments) reduce_rotation_moments(iz, full_plane, record);
  if (on_.momentum) {
    if (iz == 0) std::ranges::fill(record.prob_kr, 0.0);  // a sum over the planes
    reduce_momentum(full_plane, record);
  }
}

int FrameExtras::angular_columns() const { return tables_->angular ? tables_->angular->n_columns() : 0; }

// Per (Z, r) point the nX nY columns to the DVR grid, in passes of transform.n_columns() columns (as the potential
// step: a pass's DVR blocks stay in cache); |value|^2 summed over (x, y) per node, then over r in order.
void FrameExtras::reduce_angular(int iz, std::span<const Complex> plane, FrameRecord& record) {
  ExtraTimer timer(seconds_, angular_extra);
  const AngularTransform& transform = *tables_->angular;
  const int nodes = transform.n_theta() * transform.n_phi();
  const int ld = transform.leading_dimension();
  const std::size_t column = std::size_t(layout_.n_rot) * layout_.nxy;
  std::vector<double> per_r(std::size_t(layout_.nr) * nodes, 0.0);
#pragma omp parallel
  {
    const int n_col = transform.n_columns();
    ComplexVector dvr(std::size_t(transform.dvr_block_size()) * n_col);  // aligned: every pass runs the same plan
#pragma omp for schedule(static)
    for (int ir = 0; ir < layout_.nr; ++ir)
      for (int c0 = 0; c0 < layout_.nxy; c0 += n_col) {
        transform.fbr_to_dvr(
            plane.subspan(column * ir + std::size_t(layout_.n_rot) * c0, std::size_t(layout_.n_rot) * n_col), dvr);
        // The same additions in the same order as one pass over all columns: ixy ascending per node.
        for (int ixy = c0; ixy < c0 + n_col; ++ixy)
          for (int it = 0; it < transform.n_theta(); ++it)
            for (int ip = 0; ip < transform.n_phi(); ++ip)
              per_r[std::size_t(ir) * nodes + ip + transform.n_phi() * it] +=
                  std::norm(dvr[std::size_t(transform.dvr_block_size()) * (ixy - c0) + ip + ld * it]);
      }
  }
  for (int node = 0; node < nodes; ++node) {
    double sum = 0.0;
    for (int ir = 0; ir < layout_.nr; ++ir) sum += per_r[std::size_t(ir) * nodes + node];  // in r order
    record.prob_angle_z[std::size_t(iz) * nodes + node] = sum * layout_.cell_volume;
  }
}

// Per (Z, r) point the x-y FFT of exp(-i k_par . R) psi_full of every rotational state; |value|^2 summed over the
// states per (kx, ky) point, then over r in order.
void FrameExtras::reduce_diffraction(int iz, std::span<const Complex> plane, FrameRecord& record) {
  ExtraTimer timer(seconds_, diffraction_extra);
  const Tables& t = *tables_;
  const std::size_t column = std::size_t(layout_.n_rot) * layout_.nxy;
  std::vector<double> per_r(std::size_t(layout_.nr) * layout_.nxy, 0.0);
#pragma omp parallel
  {
    std::vector<Complex> shifted(column);
#pragma omp for schedule(static)
    for (int ir = 0; ir < layout_.nr; ++ir) {
      for (int ixy = 0; ixy < layout_.nxy; ++ixy)
        for (int irot = 0; irot < layout_.n_rot; ++irot) {
          const std::size_t k = irot + std::size_t(layout_.n_rot) * ixy;
          shifted[k] = t.shift_forward[ixy] * plane[column * ir + k];
        }
      t.xy_fft->execute(shifted);
      for (int ixy = 0; ixy < layout_.nxy; ++ixy)
        for (int irot = 0; irot < layout_.n_rot; ++irot)
          per_r[std::size_t(ir) * layout_.nxy + ixy] += std::norm(shifted[irot + std::size_t(layout_.n_rot) * ixy]);
    }
  }
  for (int ixy = 0; ixy < layout_.nxy; ++ixy) {
    double sum = 0.0;
    for (int ir = 0; ir < layout_.nr; ++ir) sum += per_r[std::size_t(ir) * layout_.nxy + ixy];
    record.pop_nm_z[std::size_t(iz) * layout_.nxy + ixy] = sum * t.diffraction_scale;
  }
}

// Per (x, y) point and rotational state the overlaps c_v = sum_r phi_vJ(r) psi(r) dr with all nR states of its J;
// |c_v|^2 summed over the states of each J, then over (x, y) in order.
void FrameExtras::reduce_vibration(int iz, std::span<const Complex> plane, FrameRecord& record) {
  ExtraTimer timer(seconds_, vibration_extra);
  const Tables& t = *tables_;
  const int nr = layout_.nr;
  const std::size_t states = std::size_t(nr) * t.n_j;
  std::vector<double> per_xy(std::size_t(layout_.nxy) * states, 0.0);
  const Index3D in_plane = layout_.plane();
#pragma omp parallel for schedule(static)
  for (int ixy = 0; ixy < layout_.nxy; ++ixy) {
    for (int irot = 0; irot < layout_.n_rot; ++irot) {
      const int j = t.j_index[irot];
      for (int v = 0; v < nr; ++v) {
        Complex overlap(0.0, 0.0);
        const std::size_t phi = (std::size_t(j) * nr + v) * nr;
        for (int ir = 0; ir < nr; ++ir) overlap += t.phi_dr[phi + ir] * plane[in_plane(irot, ixy, ir)];
        per_xy[std::size_t(ixy) * states + j + std::size_t(t.n_j) * v] += std::norm(overlap);  // [v, J]
      }
    }
  }
  for (std::size_t s = 0; s < states; ++s) {
    double sum = 0.0;
    for (int ixy = 0; ixy < layout_.nxy; ++ixy) sum += per_xy[std::size_t(ixy) * states + s];  // in (x, y) order
    record.pop_vj_z[std::size_t(iz) * states + s] = sum * t.vibration_scale;
  }
}

// Per (Z, r) point sum over (x, y) of psi^+ P2 psi and of |psi_irot|^2 J(J+1) / (2 mu r^2), then over r in order.
void FrameExtras::reduce_rotation_moments(int iz, std::span<const Complex> plane, FrameRecord& record) {
  ExtraTimer timer(seconds_, rotation_extra);
  const Tables& t = *tables_;
  const Index3D in_plane = layout_.plane();
  std::vector<double> p2_r(layout_.nr, 0.0);
  std::vector<double> rot_r(layout_.nr, 0.0);
#pragma omp parallel for schedule(static)
  for (int ir = 0; ir < layout_.nr; ++ir) {
    for (int ixy = 0; ixy < layout_.nxy; ++ixy) {
      for (const Tables::Element& e : t.p2)
        p2_r[ir] += e.value * (std::conj(plane[in_plane(e.row, ixy, ir)]) * plane[in_plane(e.col, ixy, ir)]).real();
      for (int irot = 0; irot < layout_.n_rot; ++irot)
        rot_r[ir] +=
            t.rotational_energy[irot + std::size_t(layout_.n_rot) * ir] * std::norm(plane[in_plane(irot, ixy, ir)]);
    }
  }
  double p2 = 0.0;
  double rot = 0.0;
  for (int ir = 0; ir < layout_.nr; ++ir) {
    p2 += p2_r[ir];
    rot += rot_r[ir];
  }
  record.p2_numerator_z[iz] = p2 * layout_.cell_volume;
  record.rot_numerator_z[iz] = rot * layout_.cell_volume;
}

// The r FFT of every (rotational state, x, y) column of the plane; |value|^2 summed per k_r point over the columns,
// added to the frame's sum over the planes.
void FrameExtras::reduce_momentum(std::span<const Complex> plane, FrameRecord& record) {
  ExtraTimer timer(seconds_, momentum_extra);
  Tables& t = *tables_;
  std::ranges::copy(plane, t.plane_work.begin());
  t.r_fft->execute(t.plane_work);
  const std::size_t columns = std::size_t(layout_.n_rot) * layout_.nxy;
  const double scale = layout_.cell_volume / layout_.nr;
#pragma omp parallel for schedule(static)
  for (int k = 0; k < layout_.nr; ++k) {
    double sum = 0.0;
    for (std::size_t c = 0; c < columns; ++c) sum += std::norm(t.plane_work[c + columns * k]);
    record.prob_kr[k] += sum * scale;
  }
}

void FrameExtras::reduce_specular(const SpecularWave& specular, FrameRecord& record) {
  if (!on_.momentum) return;
  ExtraTimer timer(seconds_, momentum_extra);
  Tables& t = *tables_;
  for (int iz = 0; iz < layout_.nz_specular; ++iz)
    for (int ir = 0; ir < layout_.nr; ++ir) t.specular_work[ir + std::size_t(layout_.nr) * iz] = specular(ir, iz);
  t.z_fft->execute(t.specular_work);
  const double scale = layout_.specular_cell_volume / layout_.nz_specular;
  for (int k = 0; k < layout_.nz_specular; ++k) {
    double sum = 0.0;
    for (int ir = 0; ir < layout_.nr; ++ir) sum += std::norm(t.specular_work[ir + std::size_t(layout_.nr) * k]);
    record.prob_kz_spec[k] = sum * scale;
  }
}

}  // namespace spodvr

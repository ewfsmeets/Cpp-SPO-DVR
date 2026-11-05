// The wave packet in host memory, propagated by the OpenMP operator stages. Fortran: SplitOperator.f90
#include "propagation/cpu_backend.hpp"

#include <utility>

#include "propagation/averages.hpp"

namespace spodvr {

namespace {

ScatteringShape scattering_shape(const Config& config, const RotationalBasis& basis) {
  const InputParameters& p = config.input;
  return {basis.size(), p.nx, p.ny, p.nr, p.nz};
}

AlignedVector<double> premultiplied(AlignedVector<double>&& values, double time_step) {
  premultiply_with_time_step(values, time_step);
  return std::move(values);
}

}  // namespace

CpuBackend::CpuBackend(const BackendInput& in, PotentialGrid potential, SpecularWave initial)
    : transform_(in.angular, in.basis, in.fft, in.config.input.nx * in.config.input.ny),
      premultiplied_potential_(premultiplied(std::move(potential.values), in.config.input.time_step)),
      kinetic_(scattering_shape(in.config, in.basis), in.tables, in.fft),
      potential_(scattering_shape(in.config, in.basis), in.grids.nz_potential, in.basis.initial_index(), in.tables,
                 transform_, premultiplied_potential_, potential.symmetry),
      specular_(in.config.input.nr, in.grids.z_specular.size(), in.grids.nz_potential, in.basis.initial_index(),
                *in.grids.z_specular_optical, in.grids.r_optical, in.tables, in.fft),
      propagator_(in.config.derived.steps_per_analysis, kinetic_, potential_, specular_),
      // The initial wave packet lives on the specular grid; the scattering grid starts empty.
      // Zeroed per (Z, r) group of the propagator (first touch, core/loop_schedule.hpp).
      psi_{ScatteringWave(in.basis.size(), in.config.input.nx * in.config.input.ny,
                          in.config.input.nr * in.config.input.nz,
                          xy_group_points(scattering_shape(in.config, in.basis))),
           std::move(initial)} {
  // With precomputed exp(-i V dt) the real potential is no longer needed.
  if (potential_.uses_precomputed_exponentials()) {
    premultiplied_potential_.clear();
    premultiplied_potential_.shrink_to_fit();
  }
}

std::span<const std::complex<double>> CpuBackend::z_plane(int iz) const {
  const std::size_t plane = std::size_t(psi_.scattering.block_size()) * psi_.specular.n_r();
  return std::span<const std::complex<double>>(psi_.scattering.values).subspan(plane * iz, plane);
}

std::vector<double> CpuBackend::block_densities() const { return spodvr::block_densities(psi_.scattering); }

std::unique_ptr<PropagationBackend> make_cpu_backend(const BackendInput& input, PotentialGrid potential,
                                                     SpecularWave initial) {
  return std::make_unique<CpuBackend>(input, std::move(potential), std::move(initial));
}

}  // namespace spodvr

"""The memory of the large arrays of a run, array by array, on the CPU and on the GPU.

Shadows: src/propagation/memory_estimate.cpp (estimate_memory),
src/propagation/kinetic_operators.cpp (zr_chunk_size, xy_group_points, xy_in_place_possible,
xy_momentum_shift_needed),
src/propagation/propagation_tables.cpp (momentum_shift_forward),
src/propagation/angular_transform.hpp (dvr_block_size, cpu_pass_bytes),
src/propagation/angular_transform.cpp (columns_per_pass),
src/propagation/potential_operator.hpp (default_table_limit),
src/cuda/cuda_backend.cpp (estimate_gpu_memory, staging_points, fft_work_limit, staging_bytes),
src/cuda/device_potential_operator.cpp (slab_size),
src/cuda/device_potential_operator.hpp (default_scratch_limit)

Sizes in bytes, of the program's types: double 8, complex<double> 16, int 4. With UseSymmetry the program
stores the potential and its exp table at the (x, y) representatives only. The designer cannot work out
the representatives (that needs the lattice and the grid): it takes their number from the user (the
program's SURFACE SYMMETRY block prints it) and otherwise assumes every (x, y) point, an upper bound.
"""

from dataclasses import dataclass

from program.angular_grid import AngularGrid, rotational_basis_size
from program.config import DerivedParameters, InputParameters
from program.spatial_grids import SpatialGrids

REAL_BYTES = 8  # double
COMPLEX_BYTES = 16  # std::complex<double>
INT_BYTES = 4  # int
TABLE_LIMIT = 4 << 30  # exp(-i V dt) is precomputed up to 4 GiB, else applied on the fly (default_table_limit)
XY_GROUP_LIMIT = 1 << 20  # bytes of the (x, y) blocks of one T_xy group (xy_group_points)
FFT_WORK_LIMIT = 512 << 20  # the GPU's cuFFT work area, shared by all plans (fft_work_limit)
STAGING_BYTES = 64 << 20  # the GPU's buffer for copying the wave function (staging_bytes)
SCRATCH_LIMIT = 1 << 30  # the GPU's DVR scratch of one slab of (Z, r) points (default_scratch_limit)


@dataclass(frozen=True)
class SymmetryCounts:
    """What the estimate needs of the surface symmetry (GridSymmetry): no symmetry is n_xy and 1."""

    n_representatives: int  # (x, y) points at which the potential is stored
    n_operations: int  # operations of the point group usable on the grid


@dataclass(frozen=True)
class MemoryItem:
    """One large array of a run (MemoryItem)."""

    name: str
    values: int  # number of values (per thread for the work arrays)
    bytes: int  # total size (all threads for the work arrays)
    note: str
    during_setup: bool  # allocated while the propagation is set up
    during_propagation: bool  # allocated during the time loop


@dataclass(frozen=True)
class MemoryEstimate:
    """The arrays and the totals (MemoryEstimate)."""

    items: tuple[MemoryItem, ...]
    setup_bytes: int  # while the propagation is set up (potential and exp table together)
    propagation_bytes: int  # during the time loop
    peak_bytes: int  # the larger of the two


def estimate_memory(
    p: InputParameters,
    d: DerivedParameters,
    grids: SpatialGrids,
    angular: AngularGrid,
    symmetry: SymmetryCounts,
    n_threads: int,
) -> MemoryEstimate:
    """The CPU memory of the run (estimate_memory): the potential only for a potential-only run (6), else also
    the exp(-i V dt) table, the wave functions and the work arrays of n_threads OpenMP threads."""
    n_xy = p.nx * p.ny
    n_tp = angular.n_theta * angular.n_phi
    reduced = symmetry.n_representatives < n_xy
    items = []

    # Potential on the interaction region (Z < Z_asymptotic), at the (x, y) representatives. Fortran: vdvr
    potential_values = n_tp * symmetry.n_representatives * grids.nz_potential * p.nr
    precomputed = potential_values * COMPLEX_BYTES <= TABLE_LIMIT
    propagation = d.run_mode != 6
    potential_note = ("freed after setup" if precomputed else "kept") if propagation else ""
    items.append(
        MemoryItem(
            "Potential V",
            potential_values,
            potential_values * REAL_BYTES,
            potential_note,
            during_setup=True,
            during_propagation=propagation and not precomputed,
        )
    )
    if reduced:
        # GridSymmetry: representative and operation per (x, y) point, a phi map per operation.
        map_values = 2 * n_xy + symmetry.n_operations * angular.n_phi
        note = f"{symmetry.n_representatives} of {n_xy} (x, y) points"
        items.append(MemoryItem("Surface symmetry maps", map_values, map_values * INT_BYTES, note, True, propagation))
    if propagation:
        items += propagation_items(p, grids, angular, symmetry, potential_values, precomputed, n_threads)
    return with_totals(items, peak_is_propagation=False)


def propagation_items(
    p: InputParameters,
    grids: SpatialGrids,
    angular: AngularGrid,
    symmetry: SymmetryCounts,
    potential_values: int,
    precomputed: bool,
    n_threads: int,
) -> list[MemoryItem]:
    """The exp(-i V dt) table, the wave functions and the work arrays of a propagating run."""
    n_xy = p.nx * p.ny
    items = []
    if precomputed:
        note = "precomputed (at most 4 GiB)"
        items.append(
            MemoryItem("exp(-i V dt) table", potential_values, potential_values * COMPLEX_BYTES, note, True, True)
        )
    else:
        items.append(MemoryItem("exp(-i V dt) table", 0, 0, "computed in every step (more than 4 GiB)", True, True))

    n_rot = rotational_basis_size(angular)
    scattering_values = n_rot * n_xy * p.nz * p.nr
    items.append(
        MemoryItem(
            "Wave function, scattering grid", scattering_values, scattering_values * COMPLEX_BYTES, "", True, True
        )
    )
    # The propagated specular wave plus the initial wave packet.
    specular_values = p.nr * grids.z_specular.size
    specular_bytes = 2 * specular_values * COMPLEX_BYTES
    items.append(MemoryItem("Wave function, specular grid", specular_values, specular_bytes, "2 copies", True, True))

    work_values = work_values_per_thread(p, angular, symmetry, n_rot, precomputed)
    note = f"per thread, x {n_threads} threads"
    items.append(MemoryItem("Work arrays", work_values, work_values * COMPLEX_BYTES * n_threads, note, False, True))
    return items


def work_values_per_thread(
    p: InputParameters, angular: AngularGrid, symmetry: SymmetryCounts, n_rot: int, precomputed: bool
) -> int:
    """Values of one thread's work arrays: the (r, Z) unit of T_Z + T_r, or the (x, y) group of T_xy together
    with the DVR blocks of the potential step (the stages run one after the other)."""
    n_xy = p.nx * p.ny
    zr_unit = zr_chunk_size(n_rot * p.ny) * p.nr * p.nz
    group_points = xy_group_points(n_rot, p.nx, p.ny, p.nr)
    # T_xy in place (one aligned point per group, no shift) needs no (x, y) group array.
    in_place = xy_in_place_possible(n_rot, p.nx, p.ny, group_points, xy_momentum_shift_needed(p))
    xy_group = 0 if in_place else n_rot * n_xy * group_points
    # The DVR blocks of the columns of one angular-transform pass.
    dvr_blocks = dvr_block_size(angular.n_theta, angular.n_phi) * columns_per_pass(n_xy, angular.n_theta, angular.n_phi)
    # On the fly with symmetry: exp(i A) of the representatives of one (Z, r) point.
    reduced = symmetry.n_representatives < n_xy
    exp_scratch = angular.n_theta * angular.n_phi * symmetry.n_representatives if not precomputed and reduced else 0
    return max(zr_unit, xy_group + dvr_blocks + exp_scratch)


def zr_chunk_size(n_roty: int) -> int:
    """Adjacent (irot, iy) values per (r, Z) work unit: up to 4 (one cache line), a divisor of n_rot * ny."""
    for c in range(min(4, n_roty), 1, -1):
        if n_roty % c == 0:
            return c
    return 1


def xy_group_points(n_rot: int, nx: int, ny: int, nr: int) -> int:
    """r points per (x, y) group: the largest divisor of nR whose (irot, x, y) blocks fit in 1 MiB (at least 1)."""
    block_bytes = COMPLEX_BYTES * n_rot * nx * ny
    for g in range(nr, 1, -1):
        if nr % g == 0 and g * block_bytes <= XY_GROUP_LIMIT:
            return g
    return 1


def xy_in_place_possible(n_rot: int, nx: int, ny: int, group_points: int, shift: bool) -> bool:
    """Whether T_xy transforms a group in the wave function itself: one point per group, no momentum shift, and a
    block of n_rot * nx * ny complex values that is a multiple of 64 bytes."""
    return group_points == 1 and not shift and (COMPLEX_BYTES * n_rot * nx * ny) % 64 == 0


def xy_momentum_shift_needed(p: InputParameters) -> bool:
    """Whether T_xy shifts to the frame of the initial parallel momentum. The C++ decides from the shift factors
    (some factor not exactly 1), which for a valid input is EParallel != 0; the designer has no lattice to compute
    the factors."""
    return p.e_parallel != 0.0


def dvr_block_size(n_theta: int, n_phi: int) -> int:
    """Values of one DVR block of the angular transform, padding included: (n_phi + 1) * n_theta."""
    return (n_phi + 1) * n_theta


CPU_PASS_BYTES = 58 << 10  # AngularTransform::cpu_pass_bytes


def columns_per_pass(n_xy: int, n_theta: int, n_phi: int) -> int:
    """Columns per CPU angular-transform pass: the largest divisor of n_xy whose DVR blocks fit CPU_PASS_BYTES,
    at least 1."""
    block_bytes = COMPLEX_BYTES * dvr_block_size(n_theta, n_phi)
    for n in range(n_xy, 1, -1):
        if n_xy % n == 0 and n * block_bytes <= CPU_PASS_BYTES:
            return n
    return 1


def estimate_gpu_memory(
    p: InputParameters, grids: SpatialGrids, angular: AngularGrid, n_representatives: int
) -> MemoryEstimate:
    """The GPU memory of a propagating run (estimate_gpu_memory), the cuFFT work area at its limit; the CUDA
    context itself (a few hundred MiB) is not included."""
    n_rot = rotational_basis_size(angular)
    n_xy = p.nx * p.ny
    n_tp = angular.n_theta * angular.n_phi
    table = n_tp * n_representatives * grids.nz_potential * p.nr
    slab = slab_size(p, grids.nz_potential, angular)
    scratch = n_tp * n_xy * slab
    scattering = p.nr * p.nz * n_rot * n_xy
    specular = p.nr * grids.z_specular.size
    plane = n_rot * n_xy * p.nr
    staging = staging_points(n_rot, n_xy) * n_rot * n_xy
    c = COMPLEX_BYTES
    items = [
        MemoryItem("Wave function, scattering grid", scattering, c * scattering, "", True, True),
        MemoryItem("Wave function, specular grid", specular, c * specular, "", True, True),
        MemoryItem("exp(-i V dt) table", table, c * table, "computed on the host", True, True),
        MemoryItem("DVR scratch (one slab)", scratch, c * scratch, f"{slab} (Z, r) points per slab", True, True),
        MemoryItem("cuFFT work area", FFT_WORK_LIMIT // c, FFT_WORK_LIMIT, "at most", True, True),
        MemoryItem("BK plane", plane, c * plane, "", True, True),
        MemoryItem("Copy buffer", staging, c * staging, "only while copying (restart, checkpoint)", True, True),
    ]
    return with_totals(items, peak_is_propagation=True)


def slab_size(p: InputParameters, nz_potential: int, angular: AngularGrid) -> int:
    """(Z, r) points per slab of the GPU's potential step: the largest divisor of nzpot * nR whose DVR scratch
    fits in 1 GiB (at least 1)."""
    points = nz_potential * p.nr
    bytes_per_point = COMPLEX_BYTES * p.nx * p.ny * angular.n_theta * angular.n_phi
    for n in range(points, 1, -1):
        if points % n == 0 and n * bytes_per_point <= SCRATCH_LIMIT:
            return n
    return 1


def staging_points(n_rot: int, n_xy: int) -> int:
    """(Z, r) points per chunk of the GPU's copy buffer (at least one)."""
    return max(1, STAGING_BYTES // (COMPLEX_BYTES * n_rot * n_xy))


def with_totals(items: list[MemoryItem], peak_is_propagation: bool) -> MemoryEstimate:
    """The estimate with its totals; on the GPU the peak is the propagation total (everything stays)."""
    setup = sum(item.bytes for item in items if item.during_setup)
    propagation = sum(item.bytes for item in items if item.during_propagation)
    peak = propagation if peak_is_propagation else max(setup, propagation)
    return MemoryEstimate(tuple(items), setup, propagation, peak)

"""The sizes of the angular grid (theta, phi) and of the rotational basis.

Shadows: src/grids/angular_grid.cpp (AngularGrid, make_angular_grid),
src/grids/rotational_basis.cpp (RotationalBasis)

Only the sizes are copied (they set the memory of the potential and of the wave function), not the
Gauss-Legendre points. The checks of the AngularGrid constructor are not copied either: they repeat checks
of validate_input or hold by construction (README, "Limits").
"""

from dataclasses import dataclass

from program.config import DerivedParameters, InputParameters
from program.fft_sizes import nearest_fft_size


@dataclass(frozen=True)
class AngularGrid:
    """The (theta, phi) grid sizes and the J range of the rotational basis (AngularGrid)."""

    n_theta: int  # Gauss-Legendre points (half of them for a homonuclear molecule)
    n_phi: int  # equidistant points: a valid FFT size >= 2 MJMax + 1
    j_min: int  # J range of the basis, after the parity rule of a homonuclear molecule
    j_max: int
    delta_j: int  # 2 for a homonuclear molecule (only J of JIn's parity), else 1
    mj_max: int  # MJMax of the basis


def make_angular_grid(p: InputParameters, d: DerivedParameters) -> AngularGrid:
    """The angular grid of a validated input (make_angular_grid)."""
    j_min, j_max = 0, p.j_max
    if p.homonuclear:
        # Only the J of JIn's parity are coupled: the J range starts and ends on that parity.
        if p.j_in % 2 == 0:
            if j_max % 2 != 0:
                j_max -= 1
            if j_min % 2 != 0:
                j_min += 1
        else:
            if j_max % 2 == 0:
                j_max -= 1
            if j_min % 2 == 0:
                j_min += 1
        n_theta = (j_max - j_min) // 2 + 1
        delta_j = 2
    else:
        n_theta = j_max - j_min + 1
        delta_j = 1
    n_phi = nearest_fft_size(2 * d.mj_max_basis + 1, d.phi_points_multiple)
    return AngularGrid(n_theta, n_phi, j_min, j_max, delta_j, d.mj_max_basis)


def rotational_basis_size(grid: AngularGrid) -> int:
    """The number of (J, mJ) functions: |mJ| <= min(MJMax, JMax), J = j_min .. j_max in steps delta_j, |mJ| <= J
    (RotationalBasis)."""
    mj_max = min(grid.mj_max, grid.j_max)
    size = 0
    for m in range(-mj_max, mj_max + 1):
        for j in range(grid.j_min, grid.j_max + 1, grid.delta_j):
            if abs(m) <= j:
                size += 1
    return size

"""The Z and r grids, the analysis point, the interaction region and the absorber index ranges.

Shadows: src/grids/spatial_grids.cpp (make_spatial_grids, optical_range, fortran_grid_point, add_split_points,
absorber_split_points), src/grids/uniform_grid.hpp (UniformGrid), src/grids/momentum_grid.cpp (MomentumGrid)

Grid points are `start + i * spacing` in double precision, as in the C++. The Fortran computes them in quad
precision and rounds once, so a point can differ by one ulp; `fortran_grid_point` emulates that with exact
fractions to find the points that the two programs put on different sides of an absorber's strongest end
(inputs/README.md, "The grid-point trap").
"""

from collections.abc import Callable
from dataclasses import dataclass, replace
from fractions import Fraction

from program.config import DerivedParameters, InputParameters, low_z_absorber_enabled
from program.findings import Finding, FindingList, has_errors
from program.real_compare import rge, rgt, rlt
from program.unit_conversion import PI


@dataclass(frozen=True)
class UniformGrid:
    """Points start + i * spacing, i = 0 .. size-1 (UniformGrid)."""

    start: float
    spacing: float
    size: int

    def point(self, i: int) -> float:
        return self.start + i * self.spacing

    def points(self) -> list[float]:
        return [self.point(i) for i in range(self.size)]

    @property
    def last(self) -> float:
        return self.point(self.size - 1)


@dataclass(frozen=True)
class MomentumGrid:
    """The FFT momentum grid of a coordinate grid: spacing 2 pi / (n d), largest |k| = (n/2) dk (MomentumGrid)."""

    size: int
    coordinate_spacing: float
    mass: float  # electron masses

    @property
    def spacing(self) -> float:
        return 2.0 * PI / (self.size * self.coordinate_spacing)

    @property
    def max_kinetic_energy(self) -> float:
        """Kinetic energy (hartree) at the largest momentum of the grid, k = pi / d."""
        return 0.5 / self.mass * ((self.size // 2) * self.spacing) ** 2


@dataclass(frozen=True)
class IndexRange:
    """Grid indices first .. last (inclusive) inside an absorber; empty if first > last."""

    first: int
    last: int


@dataclass(frozen=True)
class SplitGridPoint:
    """A grid point that this program and the Fortran put on different sides of an absorber's strongest end."""

    key: str  # the input key of that end, e.g. "ZOptEnd"
    strongest_end: float
    index: int
    point: float  # this program's point
    fortran_point: float  # the Fortran's point (from quad precision)
    inside_here: bool  # inside the absorber in this program


@dataclass(frozen=True)
class SpatialGrids:
    """The grids and indices of make_spatial_grids (Z and r only)."""

    z_scattering: UniformGrid
    z_specular: UniformGrid
    r: UniformGrid
    kz_scattering: MomentumGrid
    kr: MomentumGrid
    iz_analysis: int  # index of the analysis line: first point >= Zana
    z_analysis: float  # Zana moved onto that point
    nz_potential: int  # points with Z < Z_asymptotic: the 6D potential is stored for these
    z_optical: IndexRange
    z_specular_optical: IndexRange
    r_optical: IndexRange
    z_optical_low: IndexRange | None
    z_specular_optical_low: IndexRange | None
    split_points: tuple[SplitGridPoint, ...]


def make_spatial_grids(p: InputParameters, d: DerivedParameters) -> tuple[SpatialGrids | None, list[Finding]]:
    """The grids of a validated input, or None and the errors at which the program stops."""
    checks = FindingList("make_spatial_grids")
    checks.error_if(
        p.nz_specular < p.nz, "nZSpecular", "nZSpecular < nZ (the specular grid must contain the scattering grid)"
    )
    if has_errors(checks.findings):
        return None, checks.findings

    z_scattering = UniformGrid(p.z_start, p.dz, p.nz)
    z_specular = UniformGrid(p.z_start, p.dz, p.nz_specular)
    r = UniformGrid(p.r_start, p.dr, p.nr)

    # Analysis point: the first Z point >= Zana; Zana is moved onto it. validate_input already put Zana within
    # 1e-12 of a grid point, so the move is at most that. The C++ repeats the checks of validate_input here
    # (ZAna > ZOptStart is the third time), as the Fortran's Setup.f90 does.
    iz_analysis = first_index(z_scattering, lambda z: rge(z, p.z_analysis))
    if iz_analysis is None:
        checks.error("Zana", "analysis value ZAnalysis lies outside grid")
        return None, checks.findings
    z_analysis = z_scattering.point(iz_analysis)
    checks.error_if(rgt(z_analysis, p.z_opt_start), "Zana", "analysis value ZAnalysis lies inside optical potential")

    # Interaction region: Z < Z_asymptotic (a point within 1e-12 of it counts as gas phase).
    nz_potential = first_index(z_scattering, lambda z: rge(z, d.z_asymptotic))
    checks.error_if(nz_potential is None, "nZ", "the scattering grid in Z does not reach ZAsymptotic")
    if has_errors(checks.findings):
        return None, checks.findings

    grids = SpatialGrids(
        z_scattering=z_scattering,
        z_specular=z_specular,
        r=r,
        kz_scattering=MomentumGrid(p.nz, p.dz, d.total_mass),
        kr=MomentumGrid(p.nr, p.dr, d.reduced_mass),
        iz_analysis=iz_analysis,
        z_analysis=z_analysis,
        nz_potential=nz_potential,
        z_optical=optical_range(z_scattering, p.z_opt_start, p.z_opt_end),
        z_specular_optical=optical_range(z_specular, p.zsp_opt_start, p.zsp_opt_end),
        r_optical=IndexRange(0, 0) if p.nr == 1 else optical_range(r, p.r_opt_start, p.r_opt_end),
        z_optical_low=(
            optical_range(z_scattering, p.z_opt_low_start, p.z_opt_low_end)
            if low_z_absorber_enabled(p.z_opt_low_start, p.z_opt_low_end)
            else None
        ),
        z_specular_optical_low=(
            optical_range(z_specular, p.zsp_opt_low_start, p.zsp_opt_low_end)
            if low_z_absorber_enabled(p.zsp_opt_low_start, p.zsp_opt_low_end)
            else None
        ),
        split_points=(),
    )
    split = absorber_split_points(p, grids)
    for s in split:
        side = "inside" if s.inside_here else "outside"
        checks.warning(
            s.key,
            f"grid point {s.index} (Z = {s.point!r}) is {side} the absorber here but not in the "
            f"Fortran (its point is {s.fortran_point!r}); move {s.key} by about 1e-13 bohr",
        )
    return replace(grids, split_points=tuple(split)), checks.findings


def first_index(grid: UniformGrid, condition: Callable[[float], bool]) -> int | None:
    """The first index whose point satisfies `condition`, or None."""
    for i in range(grid.size):
        if condition(grid.point(i)):
            return i
    return None


def optical_range(grid: UniformGrid, start: float, end: float) -> IndexRange:
    """The points inside [start, end] with the tolerant comparisons (optical_range)."""
    first = 0
    while first < grid.size and rlt(grid.point(first), start):
        first += 1
    last = grid.size - 1
    while last >= 0 and rgt(grid.point(last), end):
        last -= 1
    return IndexRange(first, last)


QUAD_SIGNIFICAND_BITS = 113  # IEEE binary128, gfortran's REAL(16) and GCC's __float128


def round_to_quad(x: Fraction) -> Fraction:
    """x rounded to the nearest quad-precision number (113 significant bits, ties to even)."""
    if x == 0:
        return x
    exponent = floor_log2(abs(x))  # 2**exponent <= |x| < 2**(exponent + 1)
    scale = Fraction(2) ** (QUAD_SIGNIFICAND_BITS - 1 - exponent)
    return Fraction(round(x * scale)) / scale  # round() of a Fraction rounds ties to even


def floor_log2(x: Fraction) -> int:
    """The largest e with 2**e <= x, for x > 0."""
    e = x.numerator.bit_length() - x.denominator.bit_length()
    return e - 1 if Fraction(2) ** e > x else e


def fortran_grid_point(grid: UniformGrid, i: int) -> float:
    """Point i as the Fortran computes it: start + i * spacing in quad precision, rounded once to double."""
    offset = round_to_quad(Fraction(i) * Fraction(grid.spacing))
    return float(round_to_quad(Fraction(grid.start) + offset))  # float() of a Fraction rounds correctly


def add_split_points(
    split: list[SplitGridPoint], grid: UniformGrid, key: str, strongest_end: float, high_side: bool
) -> None:
    """Adds the points that the two programs put on different sides of `strongest_end`.

    The inside tests are the Fortran's plain <= (absorber strongest at its end, high_side) or >= (low-Z
    absorber, strongest at its start).
    """

    def inside(q: float) -> bool:
        return q <= strongest_end if high_side else q >= strongest_end

    for i in range(grid.size):
        point, fortran_point = grid.point(i), fortran_grid_point(grid, i)
        if inside(point) != inside(fortran_point):
            split.append(SplitGridPoint(key, strongest_end, i, point, fortran_point, inside(point)))


def absorber_split_points(p: InputParameters, g: SpatialGrids) -> list[SplitGridPoint]:
    """The split points of every absorber that is on (absorber_split_points)."""
    split: list[SplitGridPoint] = []
    add_split_points(split, g.z_scattering, "ZOptEnd", p.z_opt_end, high_side=True)
    if p.nr > 1:
        add_split_points(split, g.r, "ROptEnd", p.r_opt_end, high_side=True)
    add_split_points(split, g.z_specular, "ZSpOptEnd", p.zsp_opt_end, high_side=True)
    if g.z_optical_low is not None:
        add_split_points(split, g.z_scattering, "ZOptLowStart", p.z_opt_low_start, high_side=False)
    if g.z_specular_optical_low is not None:
        add_split_points(split, g.z_specular, "ZSpOptLowStart", p.zsp_opt_low_start, high_side=False)
    return split

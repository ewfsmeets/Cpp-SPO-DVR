"""The program's surface point groups: each group's generators, their closure, and how an operation moves a
Cartesian point (x along the first lattice vector, about the cell's origin, the top site). An operation is a
rotation by a multiple of 30 degrees, after a reflection y -> -y when it is a mirror.

Shadows: src/symmetry/point_group.cpp (PointGroup::PointGroup, group_closure), src/symmetry/surface_operation.cpp
(compose, apply_to_cartesian)
"""

import math
from dataclasses import dataclass


@dataclass(frozen=True, kw_only=True)
class SurfaceOperation:
    """A rotation by angle_twelfths x 30 degrees, after the reflection y -> -y when mirror."""

    mirror: bool
    angle_twelfths: int


def rotation(twelfths: int) -> SurfaceOperation:
    """The rotation by twelfths x 30 degrees."""
    return SurfaceOperation(mirror=False, angle_twelfths=twelfths)


def mirror(twelfths: int) -> SurfaceOperation:
    """The mirror operation of angle twelfths x 30 degrees (the reflection, then that rotation)."""
    return SurfaceOperation(mirror=True, angle_twelfths=twelfths)


# The generators of each group the program knows (C3v: a mirror at 90 degrees, p3m1 about a top site).
GENERATORS = {
    "C2v": [rotation(6), mirror(0)],
    "C4v": [rotation(3), mirror(0)],
    "C6v": [rotation(2), mirror(0)],
    "C3v": [rotation(4), mirror(6)],
    "C3": [rotation(4)],
    "none": [],
}


def compose(a: SurfaceOperation, b: SurfaceOperation) -> SurfaceOperation:
    """The operation a after b."""
    sign_a = -1 if a.mirror else 1
    return SurfaceOperation(
        mirror=a.mirror != b.mirror, angle_twelfths=(a.angle_twelfths + sign_a * b.angle_twelfths) % 12
    )


def group_closure(generators: list[SurfaceOperation]) -> list[SurfaceOperation]:
    """Every operation the generators make, the identity included."""
    operations = [rotation(0)]
    grown = True
    while grown:
        grown = False
        for g in generators:
            for op in list(operations):
                product = compose(g, op)
                if product not in operations:
                    operations.append(product)
                    grown = True
    return operations


def apply_to_cartesian(g: SurfaceOperation, x: float, y: float) -> tuple[float, float]:
    """The Cartesian point (x, y) moved by g."""
    beta = 2.0 * math.pi * g.angle_twelfths / 12.0
    c, s = math.cos(beta), math.sin(beta)
    y = -y if g.mirror else y
    return c * x - s * y, s * x + c * y

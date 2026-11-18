"""The surface point groups: their order and the multiple the number of phi points must have.

Shadows: src/symmetry/point_group.cpp (PointGroup, group_closure, phi_points_multiple)

The C++ builds each group from its generators (rotations and mirrors in steps of 30 degrees, "twelfths")
and takes phi_points_multiple = 12 / gcd(12, the angles of all its operations): then every operation maps
the phi grid onto itself. An operation of angle b maps phi to b + phi (rotation) or b - phi (mirror); a mirror
of angle b reflects in the line at b/2, so C3v's mirror 6 is the program's "mirror at 90 deg". The designer
needs only the results, worked out here from the generators (surface_operation.hpp):

| group | rotations (twelfths) | mirrors (twelfths) | gcd | multiple | operations |
|-------|----------------------|--------------------|-----|----------|------------|
| C2v   | 0, 6                 | 0, 6               | 6   | 2        | 4          |
| C4v   | 0, 3, 6, 9           | 0, 3, 6, 9         | 3   | 4        | 8          |
| C6v   | 0, 2, ..., 10        | 0, 2, ..., 10      | 2   | 6        | 12         |
| C3v   | 0, 4, 8              | 2, 6, 10           | 2   | 6        | 6          |
| C3    | 0, 4, 8              | -                  | 4   | 3        | 3          |
| none  | 0                    | -                  | 12  | 1        | 1          |
"""

PHI_POINTS_MULTIPLE = {"C2v": 2, "C4v": 4, "C6v": 6, "C3v": 6, "C3": 3, "none": 1}
GROUP_ORDER = {"C2v": 4, "C4v": 8, "C6v": 12, "C3v": 6, "C3": 3, "none": 1}

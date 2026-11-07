| Quantity | Fortran | C++ | Result |
|---|---|---|---|
| Z, r, x, y end points (`.pes` header, full precision) | r1 0.5, rn 8.3, x1 0, xn 4.43667, y1 0, yn 4.43667, z1 0, zn 14.88 | same | 8/8 bitwise identical |
| θ nodes (15, full precision) | `.pes` | helper grids.txt | 8/15 identical, max 3.3e-16 rad (6 ulp) |
| φ nodes (60) | `.pes` | helper grids.txt | 60/60 identical |
| nφ, nθ, J range, (x, y) points computed | 60, 15, 0..28 step 2, 19 | same | identical |
| Momentum grids dk_Z, dk_r, dk_x, dk_y | 6 decimals | %.16e | identical 8 (item 3) |

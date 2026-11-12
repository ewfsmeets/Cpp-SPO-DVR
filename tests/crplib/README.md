# Tests of the CRP library (`src/crplib/`)

The library's own tests: one program per `test_*.cpp`, linked only against `spodvr_crplib`, run
by CTest with the prefix `CrpLib.` (`ctest --test-dir build -R CrpLib`). No real data set is
needed; [docs/crp_library.md](../../docs/crp_library.md) "Testing without a real data set"
explains the approach.

| File | Tests |
|---|---|
| `test_list_directed_reader.cpp` | the subset of Fortran list-directed input the data files use, and the rejected forms |
| `test_units.cpp` | unit factors and the Fortran composition of conversions |
| `test_spline.cpp`, `test_periodic_spline_3d.cpp` | spline kernels: exact for cubics, derivatives, bitwise vs the Fortran |
| `test_fourier_tables.cpp` | the Fourier tables: bitwise vs the Fortran, inverse of their site bases |
| `test_dataset.cpp` | reading, writing, validation; the committed data sets equal the generator's output |
| `test_atom_surface.cpp` | atom-surface PES: exactness at the sites, symmetry, gradients, STOP paths, bitwise vs the Fortran |
| `test_site_interpolation.cpp` | the angular interpolation at each site type |
| `test_molecule_surface.cpp` | molecule-surface PES: impulse responses, the pairwise model, gradients, gas phase, switching, bitwise vs the Fortran |
| `test_qd_wrapper.cpp` | the quantum-dynamics wrapper: clamps, skewed coordinates, energy shift, range report |
| `test_cluster_pruning.cpp` | speed-up S2 (the cluster atoms within reach), bitwise |
| `test_real_data_set.cpp` | the real data set `potentials/H2Cu111SRP_C3v.dir/H2Cu111SRP_C3v` (H₂/Cu(111) SRP48) |

Shared headers: `crp_synthetic.hpp` (cut geometry, the pairwise model, synthetic data sets),
`reference_basis.hpp` (the Fourier bases and site positions, written down independently of the
library),
`fortran_reference.hpp` (reading the driver output, `EXPECT_BITS`), `temporary_directory.hpp`.
The program-side adapter is tested in `tests/test_crp_potential.cpp`.

## Reference data (`reference/`)

| Case | Data set |
|---|---|
| `kernels` | no data set: Fourier tables and spline kernels on generated inputs |
| `c6v_pairwise`, `c3v_pairwise`, `c4v_pairwise` | the pairwise model (bohr, Hartree; CutsInterpolation = F) |
| `c6v_pairwise_angstrom` | the same in Å and eV (unit conversion of data, slopes and gradients) |
| `c4v_negative_spacing` | C4v with a negative (relative) interlayer spacing |
| `<sym>_impulse_<k>` | molecular cut k equal to 1, everything else 0 (CutsInterpolation = T) |
| `c3v_missing_extra_cut` | C3v pairwise with two extra cuts: one on a shorter r grid (narrows Rmax), one without a file (skipped) |
| `cu111_srp48` | the real data set in `potentials/` (only the point lists and the results; no `dataset/` copy) |

Each case holds the data set (`dataset/`), the point lists (`*.in`, written by
`crplib_make_dataset`) and the Fortran driver's results (`*.out`, hexadecimal and decimal), plus
the compiler and glibc versions (`VERSIONS`). The `*_pairwise` data sets are also used by the gate
cases `tests/reference/crp_*_synthetic_small`.

Regenerate a case (data set, point lists and Fortran results) after a change to the generator:

```
cmake --build build
tools/fortran_reference/run_crp_reference.sh <case>      # or: kernels
```

The driver is linked against the objects of the Fortran reference install and checks that the
installed gfortran is the one that built them (`BUILDINFO.tsv`,
[tools/fortran_reference/README.md](../../tools/fortran_reference/README.md)); the results are only
bitwise comparable with that build.
`ReferenceDataSets.MatchTheGenerator` fails when a committed data set differs from what the
generator writes.

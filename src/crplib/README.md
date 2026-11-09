# `src/crplib/`: the CRP library

A self-contained C++ port of the Fortran CRP library of M. Wijzenbroek
(`SPO-DVR/procfiles.dir/SOURCES.dir/CRPlibrary.dir/`): CRP potential energy surfaces of a diatomic
molecule at a C3v, C4v or C6v surface, with analytic gradients. CMake target `spodvr_crplib`,
namespace `crplib`, standard library only, strict floating-point flags in every build. It keeps
the Fortran library's own constants and unit factors. SPO-DVR uses it only through
`src/potential/crp_potential.cpp`.

Documentation (API, data-set format, behaviour outside the data, cut tables):
[docs/crp_library.md](../../docs/crp_library.md). Deviations from the Fortran:
[docs/fortran_deviations.md](../../docs/fortran_deviations.md). Tests: `tests/crplib/`.

| Files | Content | Fortran |
|---|---|---|
| `dataset.*`, `list_directed_reader.*` | the data set as written in its files; reading (list-directed input), writing, validation | `CRP_Wrapper.f90` set-up, `Read_*_Data` |
| `units.*`, `constants.hpp`, `errors.hpp` | unit codes and factors, constants (incl. the values gfortran folded), errors | `UnitConversion.f90`, `CRP_Constants.f90` |
| `spline.*`, `periodic_spline_3d.*`, `fourier_tables.*` | spline kernels, the C4v periodic 3D spline, the Fourier tables | `Spline.f90`, `FOUR6D`, `FOUR3D` |
| `atom_surface*`, `hexagonal_atom_surface.cpp`, `square_atom_surface.cpp` | atom-surface PES | `CRP_C{3,4,6}v_AtomSurface.f90` |
| `site_interpolation.*`, `molecule_surface*`, `hexagonal_interaction.cpp`, `square_interaction.cpp` | molecule-surface PES: cuts, site functions, Fourier interpolation, gas phase and switching | `CRP_C{3,4,6}v_MoleculeSurface.f90`, `CRP_GasPhase.f90`, `CRP_General.f90` |
| `wrapper.*`, `qd_wrapper.*` | entry points with unit conversion; the quantum-dynamics wrapper | `CRP_Wrapper.f90`, `CRP_QD_Wrapper.f90` |

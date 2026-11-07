| # | Quantity | Fortran (label: value) | C++ (label: value) | Kind | Δ (C++ − F) | Tolerance | Verdict |
|---|---|---|---|---|---|---|---|
| 1 | lattice constant x | The surface lattice constant in X (bohr): `0.48400000E+01` | Lattice constants (x, y): `4.840000` | printed | 0.00e+00 | 5.55e-07 | printed |
| 1 | lattice constant y | The surface lattice constant in Y (bohr): `0.48400000E+01` | Lattice constants (x, y): `4.840000` | printed | 0.00e+00 | 5.55e-07 | printed |
| 1 | Z asymptotic | The asymptotic value in Z (bohr): `0.15000000E+02` | Z asymptotic: `15.000000` | printed | 0.00e+00 | 1.01e-06 | printed |
| 1 | skew angle | The unit cell skewing angle (degrees): `0.60000000E+02` | Skewing angle: `60.000000` | printed | 0.00e+00 | 1.01e-06 | printed |
| 1 | point group | The unit cell point group: `C6v` | Point group: `C6v` | text |  | exact | identical |
| 1 | r cut-off | The cutoff value in r (bohr): `0.10000000E+04` | -: `—` | n/a |  |  | n/a (UseRCutOff is not ported (forced off by the Fortran ProcessInput; fortran_deviations.md)) |
| 1 | projection scheme | The projection operator scheme used: `using two grids all the way` | ProjOp (input; only 2 is accepted): `using two grids all the way` | text |  | exact | identical |
| 1 | switching function | Using a switching function for the potential in Z: `F` | Switching to gas phase: `F` | text |  | exact | identical |
| 1 | only/also fine analysis | Only / Also performing a fine analysis: `FF` | Run mode: `FF` | text |  | exact | identical |
| 1 | analysed J, v, diffraction order | Maximum J/V/diffraction order analysed: `determined by the program (see info Analysis)` | (JAna, VAna, MDifAna absent from input.json): `determined by the program (see info Analysis)` | text |  | exact | identical |
| 1 | potential limits, mixing | Imposing a minimum / maximum value; mixing alfa: `FF0.00000000E+00` | (VMin, VMax absent; PESMixingAlfa rejected): `—` | n/a |  |  | n/a (C++ has no mixing (NotSupportedError); no limits in either input) |
| 1 | mJ-independent grid | A general MjIn independant grid will be used: `yes` | -: `—` | n/a |  |  | n/a (Fortran setupdvr.f90 message; the C++ grid is always MJIn independent) |
| 2 | number of phi points | Number of phi's: `60` | Phi points: `60` | int | 0 | exact | identical |
| 2 | number of theta points | Number of theta's: `15` | Theta points: `15` | int | 0 | exact | identical |
| 2 | J max in basis | JMax in basis set: `28` | J range (step): `28` | int | 0 | exact | identical |
| 2 | J min in basis | JMin in basis set: `0` | J range (step): `0` | int | 0 | exact | identical |
| 2 | J step in basis | Delta J in basis set: `2` | J range (step): `2` | int | 0 | exact | identical |
| 2 | phi step | Delta value in phi coordinate: `0.10471975511965977` | phi[1] - phi[0] (helper grids.txt): `0.10471975511965977` | full | 0 | 1e-13 rel. | identical (Fortran prints 17 significant digits) |
| 2 | (x, y) points with V computed | Number of points in (x,y) for which the potential is calculated: `19` | (x, y) points computed / total: `19` | int | 0 | exact | identical |
| 2 | (x, y) points in .pes vs C++ representatives | .pes n_xy_computed: `19` | C++ representatives: `19` | int | 0 | exact | identical (different fundamental regions: 5 of 19 points in common) |
| 2 | Z scattering, first | .pes header z1: `0.0` | helper grids.txt: `0` | full | 0 | 1e-13 rel. | identical |
| 2 | Z(nzpot), last point of the potential | .pes header zn: `14.879999999999999` | helper grids.txt: `14.879999999999999` | full | 0 | 1e-13 rel. | identical |
| 2 | r first | .pes header r1: `0.5` | helper grids.txt: `0.5` | full | 0 | 1e-13 rel. | identical |
| 2 | r last | .pes header rn: `8.3` | helper grids.txt: `8.3000000000000007` | full | 0 | 1e-13 rel. | identical |
| 2 | x first | .pes header x1: `0.0` | helper grids.txt: `0` | full | 0 | 1e-13 rel. | identical |
| 2 | x last | .pes header xn: `4.4366666666666665` | helper grids.txt: `4.4366666666666665` | full | 0 | 1e-13 rel. | identical |
| 2 | y first | .pes header y1: `0.0` | helper grids.txt: `0` | full | 0 | 1e-13 rel. | identical |
| 2 | y last | .pes header yn: `4.4366666666666665` | helper grids.txt: `4.4366666666666665` | full | 0 | 1e-13 rel. | identical |
| 2 | theta nodes (15) | .pes header thetas: `15 values` | helper grids.txt / log theta table: `—` | n/a |  |  | n/a (see full comparison: max \|d\| 3.33e-16 rad, 8/15 identical, max 6 ulp) |
| 3 | dk_Z scattering | grid interval size in kz dir. (scat. grid): `0.327249` | dkz (scattering / specular): `0.3272492347489368` | printed | 2.35e-07 | 5.05e-07 | printed |
| 3 | dk_Z specular | grid interval size in kz dir. (spec. grid): `0.204531` | dkz (scattering / specular): `0.20453077171808551` | printed | -2.28e-07 | 5.05e-07 | printed |
| 3 | initial k_Z | initial momentum in z direction: `13.205180` | Initial kz: `13.205179688916409` | printed | -3.11e-07 | 5.05e-07 | printed |
| 3 | dk_r | grid interval size in kr dir.: `0.785398` | dkr: `0.78539816339744828` | printed | 1.63e-07 | 5.05e-07 | printed |
| 3 | dk_x | grid interval size in kx dir.: `1.499008` | dkx, dky: `1.4990077390365437` | printed | -2.61e-07 | 5.05e-07 | printed |
| 3 | dk_y | grid interval size in ky dir.: `1.499008` | dkx, dky: `1.4990077390365437` | printed | -2.61e-07 | 5.05e-07 | printed |
| 3 | initial k_x | initial momentum in reciprocal x direction: `0.000000` | Parallel kinetic energy (hartree): `0` | printed | 0.00e+00 | 5.05e-07 | printed (C++ prints the parallel kinetic energy (0 for k_x = k_y = 0)) |
| 3 | initial k_y | initial momentum in reciprocal y direction: `0.000000` | Parallel kinetic energy (hartree): `0` | printed | 0.00e+00 | 5.05e-07 | printed |
| 3 | nmoder | For backward compatibility: nmoder: `4` | -: `—` | n/a |  |  | n/a (legacy diagnostic, no C++ counterpart) |
| 3 | analysis line Z | asymptotic analysis at z: `15.000000` | Analysis at Z (index): `15.000000` | printed | 0.00e+00 | 1.01e-06 | printed |
| 3 | flux analysis | No flux analysis performed.: `no` | (no FluxAnalysis key): `no` | text |  | exact | identical |
| 3 | nzpot | nzpot: `125` | Z points in interaction region: `125` | int | 0 | exact | identical |
| 3 | nzchsp | nzchsp: `126` | Optical potential Z (indices) (0-based): `125` | int | 0 | exact | identical |
| 3 | z(nzchsp) | z(nzchsp): `15.000000` | z_scattering[125] (helper grids.txt): `15` | printed | 0.00e+00 | 5.05e-07 | printed |
| 3 | nzchmx | nzchmx: `160` | Optical potential Z (indices) (0-based): `159` | int | 0 | exact | identical |
| 3 | z(nzchmx) | z(nzchmx): `19.080000` | z_scattering[159] (helper grids.txt): `19.079999999999998` | printed | 0.00e+00 | 5.05e-07 | printed |
| 3 | nschsp | nschsp: `185` | Optical potential Z specular (indices) (0-based): `184` | int | 0 | exact | identical |
| 3 | z(nschsp) | z(nschsp): `22.080000` | z_specular[184] (helper grids.txt): `22.079999999999998` | printed | 0.00e+00 | 5.05e-07 | printed |
| 3 | nschmx | nschmx: `256` | Optical potential Z specular (indices) (0-based): `255` | int | 0 | exact | identical |
| 3 | z(nschmx) | z(nschmx): `30.600000` | z_specular[255] (helper grids.txt): `30.599999999999998` | printed | -3.55e-15 | 5.05e-07 | printed |
| 3 | nrchsp | nrchsp: `19` | Optical potential r (indices) (0-based): `18` | int | 0 | exact | identical |
| 3 | r(nrchsp) | r(nrchsp): `4.100000` | r[18] (helper grids.txt): `4.1000000000000005` | printed | 8.88e-16 | 5.05e-07 | printed |
| 3 | nrchmx | nrchmx: `40` | Optical potential r (indices) (0-based): `39` | int | 0 | exact | identical |
| 3 | r(nrchmx) | r(nrchmx): `8.300000` | r[39] (helper grids.txt): `8.3000000000000007` | printed | 0.00e+00 | 5.05e-07 | printed |
| 3 | max kinetic energy in r | max kinetic energy in r (ev): `3.654734` | Max. kinetic energy in r: `3.654734` | printed | 0.00e+00 | 1.01e-06 | printed |
| 3 | max analysable total energy | the max total energy ... (eV): `4.274357` | Max. analysis energy: `4.274357` | printed | 0.00e+00 | 1.01e-06 | printed |
| 3 | ectmin | ectmin (ev): `0.400000` | 95% confidence interval: `0.400000` | printed | 0.00e+00 | 1.01e-06 | printed |
| 3 | ectmax | ectmax (ev): `0.950000` | 95% confidence interval: `0.950000` | printed | 0.00e+00 | 1.01e-06 | printed |
| 6 | initial state: energy initial v,j state | energy initial v,j state: `-4.477037` | energy: `-4.4770373027650265` | printed | -3.03e-07 | 5.05e-07 | printed |
| 6 | initial state: potential energy | potential energy: `-4.609423` | potential energy: `-4.6094230053894281` | printed | -5.39e-09 | 5.05e-07 | printed |
| 6 | initial state: rotational energy | rotational energy: `0.000000` | rotational energy: `0` | printed | 0.00e+00 | 5.05e-07 | printed |
| 6 | initial state: vibrational energy | vibrational energy: `0.132386` | vibrational energy: `0.1323857026244013` | printed | -2.97e-07 | 5.05e-07 | printed |
| 6 | initial state: rovibrational energy | rovibrational energy: `0.132386` | rovibrational kinetic energy: `0.1323857026244013` | printed | -2.97e-07 | 5.05e-07 | printed |
| 6 | initial state: expectation value r | expectation value r: `1.442768` | <r>: `1.4427680834424175` | printed | 8.34e-08 | 5.05e-07 | printed |
| 6 | initial rovibrational energy (hartree) | fort.21 record 2: eviroi: `-0.16452800307095536` | gas_phase_states.dat E(0,0): `-0.16452800307095652` | full | -1.17e-15 (42 ulp) | 1e-12 rel. | rounding (criterion 1e-12 relative (gate 2: LAPACK vs Jacobi)) |
| 7 | warning: basis too small | WARNING: Diffractional basis too small to describe all open channels!: `yes` | Warning: the diffraction basis is too small to describe all open channels.: `yes` | text |  | exact | identical |
| 7 | diffraction channels in basis | ndftow: `144` | Diffraction channels analysed (ndifa): '144 of N': `144` | int | 0 | exact | identical |
| 7 | diffraction channels analysed | ndifa: `144` | Diffraction channels analysed (ndifa): `144` | int | 0 | exact | identical |
| 7 | highest diffraction order | maxord: `11` | Highest diffraction order analysed (mdifa): `11` | int | 0 | exact | identical |
| 7 | largest number of channels of one order | maxnum: `30` | derived: max count per order of the 144 channels (bk_coefficients.bin labels): `30` | int | 0 | exact | identical |
| 8 | rotation channels in basis | nrotw: `435` | Rotational basis functions: `435` | int | 0 | exact | identical |
| 8 | rotational states analysed | nana: `66` | Rotational states analysed (nana): `66` | int | 0 | exact | identical |
| 8 | initial rotational state | ichn: `211` | Initial channel (index: J, mJ) (0-based): `210` | int | 0 | exact | identical |
| 8 | initial rotational state (fort.21) | fort.21 record 3: ichn: `211` | Initial channel (0-based): `210` | int | 0 | exact | identical |
| 9 | rotation-diffraction channels in basis | nstatw: `62640` | derived: rotational basis x diffraction channels: `62640` | int | 0 | exact | identical |
| 9 | rotation-diffraction states analysed | nstata: `9504` | Rotation-diffraction states (nstata): `9504` | int | 0 | exact | identical |
| 9 | initial rot-dif channel | inchna: `4321` | Initial channel (index) (0-based): `4320` | int | 0 | exact | identical |
| 10 | CAP coefficient: Z | Optical potential in Z absorbs with coefficient: `0.12709614E-01` | CAP strengths (log, %.8e): `1.27096140e-02` | printed | 0.00e+00 | 5.55e-10 | printed |
| 10 | CAP coefficient (full precision C++) | Optical potential in Z absorbs with coefficient: `0.12709614E-01` | bk_coefficients.bin: ZOptFactor (effective): `0.012709613996342364` | printed | -3.66e-12 | 5.05e-10 | printed |
| 10 | CAP coefficient: r | Optical potential in r absorbs with coefficient: `0.24559644E-01` | CAP strengths (log, %.8e): `2.45596439e-02` | printed | -1.00e-10 | 5.55e-10 | printed |
| 10 | CAP coefficient (full precision C++) | Optical potential in r absorbs with coefficient: `0.24559644E-01` | bk_coefficients.bin: ROptFactor (effective): `0.024559643875052824` | printed | -1.25e-10 | 5.05e-10 | printed |
| 10 | CAP coefficient: Z specular | Optical potential in Z specular absorbs with coefficient: `0.10193608E-01` | CAP strengths (log, %.8e): `1.01936078e-02` | printed | -2.00e-10 | 5.55e-10 | printed |
| 10 | CAP coefficient (full precision C++) | Optical potential in Z specular absorbs with coefficient: `0.10193608E-01` | bk_coefficients.bin: ZSpOptFactor (effective): `0.010193607830697777` | printed | -1.69e-10 | 5.05e-10 | printed |
| 11 | norm of the initial packet | Norm initial wave packet (Z,r): `0.10000000E+01` | Norm of the initial wave packet: `1.0000000000000002` | printed | 2.22e-16 | 5.05e-08 | printed |
| 11 | \|psi(k)\| at k = 0 | at k = 0: `0.14080155E-09` | \|psi(k)\| at k = 0: `1.4079984938068941e-10` | printed | -1.70e-15 | 5.05e-18 | explained (reproduced: the correctly rounded DFT of the samples with the Fortran's alpha and k_Z (11 and 3 ulp from the C++ values, quad-precision site) gives the Fortran value within its half unit (item11_psi_k.py)) |
| 11 | \|psi(k)\| at the grid edge | at the end/beginning of the grid: `0.30196069E-09` | \|psi(k)\| at the grid edge: `3.0196243966612511e-10` | printed | 1.75e-15 | 5.05e-18 | explained (reproduced: the correctly rounded DFT of the samples with the Fortran's alpha and k_Z (11 and 3 ulp from the C++ values, quad-precision site) gives the Fortran value within its half unit (item11_psi_k.py)) |
| 11 | \|psi(k)\| maximum | maximum value: `0.53214584E+00` | \|psi(k)\| maximum: `0.53214583713275454` | printed | -2.87e-09 | 5.05e-09 | printed |
| 11 | iLow | iLow: `126` | Allowed region (indices) (0-based): `125` | int | 0 | exact | identical |
| 11 | iHigh | iHigh: `184` | Allowed region (indices) (0-based): `183` | int | 0 | exact | identical |
| 11 | norm between iLow and iHigh | Norm of the initial wave packet between iLow and iHigh: `0.10000000E+01` | Norm inside the allowed region: `1.0000000000000002` | printed | 2.22e-16 | 5.05e-08 | printed |
| 11 | initial k_Z (fort.21, full) | fort.21 record 2: rkz: `13.205179688916404` | Initial kz (log %.16e): `13.205179688916409` | full | 5.33e-15 (3 ulp) | 1e-13 rel. | rounding |
| 11 | packet alpha (fort.21, full) | fort.21 record 2: alfwv: `1.9767175470965832` | Wave packet alpha (log %.16e): `1.9767175470965856` | full | 2.44e-15 (11 ulp) | 1e-13 rel. | rounding |
| 11 | translational energy (fort.21, full) | fort.21 record 3: enrgy: `0.6457207001484486` | Translational energy (log %.16e): `0.64572070014844907` | full | 4.44e-16 (4 ulp) | 1e-13 rel. | rounding |
| 11 | Z0 (fort.21) | fort.21 record 2: z0: `18.0` | input Z0: `18` | full | 0 | 1e-13 rel. | identical |
| 11 | analysis line (fort.21) | fort.21 record 2: ZAnalysis: `15.0` | Analysis at Z: `15` | full | 0 | 1e-13 rel. | identical |
| 11 | analysis time step (fort.21) | fort.21 record 2: deltat: `20.0` | dtAna (bk header): `20` | full | 0 | 1e-13 rel. | identical |
| 12 | tolerance for fixing the specular analysis line | Tolerance for fixing the dynamical analysis surface: `0.10000000E-06` | tolerance for fixing it: `9.9999999999999982e-08` | printed | -1.32e-23 | 5.05e-15 | printed |
| 12 | jana | Maximum j in analysis: jana: `10` | Highest j analysed (jana): `10` | int | 0 | exact | identical |
| 12 | nviba | Maximum v in analysis: nviba: `2` | Vibrational levels analysed (nviba): `2` | int | 0 | exact | identical |
| 12 | mdifa | Maximum diffraction order in analysis: mdifa: `11` | Highest diffraction order analysed (mdifa): `11` | int | 0 | exact | identical |
| 12 | number of analysis energies | per-energy blocks in stdout: `111` | Analysis energies: number: `111` | int | 0 | exact | identical |
| 12 | specular analysis line start | iHigh (Z = zsp[iHigh-1]): `184` | Analysis surface on the specular grid: starts at Z: `183` | int | 0 | exact | identical |

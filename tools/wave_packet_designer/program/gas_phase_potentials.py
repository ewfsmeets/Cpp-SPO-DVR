"""The gas-phase potential V(r) of the analytic potentials: what the program uses for the free molecule.

Shadows: src/potential/dai_light_cu111_potential.cpp (coulomb, exchange, DaiLightCu111Potential),
src/potential/hard_wall_potential.cpp (gas_phase, HardWallPotential),
src/potential/potential_grid.cpp (gas_phase_potential, impose_limits)

The program evaluates the potential at Z = Z_asymptotic for every point of the input's r grid (in eV, then
converted to hartree) and clips it at VMin/VMax when those are imposed. There:
- Dai-Light (H2Cu111DaiLight) is the H-H pair term of its LEPS form, U_M + A_M: with the Sato parameter it adds
  up to a Morse potential D (e^(-2 alpha (r - r0)) - 2 e^(-alpha (r - r0)));
- the hard wall (HardWallPES) is its Rydberg-type gas-phase function.
The expressions and the order of operations are the C++'s, so the values are bitwise the program's
(tests/test_gas_phase.py checks Dai-Light bitwise against the program's printed table, the hard wall against
the Fortran's printed table within 1e-14 of max(|V|, 1 eV)). Not copied: the Persson potential (its gas phase is the
full LEPS at Z = 8 bohr) and the CRP data sets; for them gas_phase_potential returns None.
"""

import math

from program.config import InputParameters
from program.real_compare import seq
from program.spatial_grids import UniformGrid
from program.unit_conversion import ev_to_hartree, hartree_to_ev

# The H2 pair term of the Dai-Light LEPS potential (hartree, bohr). Fortran: DAA, ALPAA, R0, DELAA
DAI_LIGHT_DEPTH = 0.1744
DAI_LIGHT_ALPHA = 1.028
DAI_LIGHT_R0 = 1.40
DAI_LIGHT_DELTA = -0.15


def coulomb(x: float, d: float, delta: float, alpha: float, x0: float) -> float:
    """Coulomb integral of a LEPS pair term with Sato parameter delta (hartree). Fortran: UPOT"""
    return (
        d
        / 4.0
        / (1.0 + delta)
        * ((3.0 + delta) * math.exp(-2.0 * alpha * (x - x0)) - (2.0 + 6.0 * delta) * math.exp(-alpha * (x - x0)))
    )


def exchange(x: float, d: float, delta: float, alpha: float, x0: float) -> float:
    """Exchange integral of a LEPS pair term with Sato parameter delta (hartree). Fortran: APOT"""
    return (
        d
        / 4.0
        / (1.0 + delta)
        * ((1.0 + 3.0 * delta) * math.exp(-2.0 * alpha * (x - x0)) - (6.0 + 2.0 * delta) * math.exp(-alpha * (x - x0)))
    )


def dai_light_gas_phase(r: float) -> float:
    """V(r) of the Dai-Light potential at Z >= Z_asymptotic (eV): its H-H pair term U_M + A_M."""
    u = coulomb(r, DAI_LIGHT_DEPTH, DAI_LIGHT_DELTA, DAI_LIGHT_ALPHA, DAI_LIGHT_R0)
    a = exchange(r, DAI_LIGHT_DEPTH, DAI_LIGHT_DELTA, DAI_LIGHT_ALPHA, DAI_LIGHT_R0)
    return hartree_to_ev(u + a)


def hard_wall_gas_phase(r_in: float) -> float:
    """V(r) of the hard-wall test potential (eV): an H2 Rydberg-type function. The real exponents are kept
    (r**2.0, r**3.0), as in the C++."""
    d, a1, a2, a3, a4, re = 4.8286, 2.282, 1.555, 0.7533, 2.23, 1.40
    r = r_in - re
    return -d * (1.0 + a1 * r + a2 * r**2.0 + a3 * r**3.0) * math.exp(-a4 * r) + d


GAS_PHASE = {"H2Cu111DaiLight": dai_light_gas_phase, "HardWallPES": hard_wall_gas_phase}


def has_gas_phase(system_id: str) -> bool:
    """True if the gas-phase potential of this SystemID is copied here."""
    return any(seq(system_id, known) for known in GAS_PHASE)


def gas_phase_potential(p: InputParameters, r_grid: UniformGrid) -> list[float] | None:
    """V(r) in hartree on the points of `r_grid`, clipped at VMin/VMax when imposed, as the program computes it
    (gas_phase_potential, impose_limits); None for a potential not copied here."""
    function = next((f for name, f in GAS_PHASE.items() if seq(p.system_id, name)), None)
    if function is None:
        return None
    values = [ev_to_hartree(function(r_grid.point(i))) for i in range(r_grid.size)]
    if p.impose_v_min:
        values = [max(v, ev_to_hartree(p.v_min)) for v in values]
    if p.impose_v_max:
        values = [min(v, ev_to_hartree(p.v_max)) for v in values]
    return values

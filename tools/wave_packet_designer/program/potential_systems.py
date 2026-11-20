"""The system parameters the designer needs from each potential: Z_asymptotic, the point group and, where
built in, the lattice constants.

Shadows: src/potential/potential_factory.cpp (potential_is_implemented),
src/potential/hard_wall_potential.cpp (z_asymptotic, HardWallPotential),
src/potential/dai_light_cu111_potential.cpp (z_asymptotic, DaiLightCu111Potential),
src/potential/persson_cu111_potential.cpp (z_asymptotic, PerssonCu111Potential)

The three analytic potentials have fixed values. A CRP data set (SystemID CRPlibrary, also under the SCM
correction CRPDistrtSurface) carries its own: the user gives Z_asymptotic and the point group (the program's
set-up output prints both); its lattice constants (with CRPDistrtSurface: the slab snapshot's) stay unknown
here, so the write2DPES upper bounds on x and y are not checked for it.
"""

from dataclasses import dataclass

from program.real_compare import seq
from program.unit_conversion import angstrom_to_bohr


@dataclass(frozen=True)
class SystemParameters:
    """What the designer uses of a potential's SystemParameters."""

    z_asymptotic: float  # bohr: beyond it V is the gas-phase V(r); the 6D potential is stored below it
    point_group: str  # C2v, C4v, C6v, C3v, C3 or none
    lattice_a: float | None = None  # bohr (alat), None when unknown here (a CRP data set, an SCM slab)
    lattice_b: float | None = None  # bohr (blat)


# The potentials with built-in parameters (hard_wall_id, dai_light_id, persson_id of potential_factory.cpp).
ANALYTIC_SYSTEMS = {
    "HardWallPES": SystemParameters(7.0, "C2v", 6.239660, 5.239660),
    "H2Cu111DaiLight": SystemParameters(15.0, "C6v", 4.84, 4.84),
    # Persson: cc = 2.5173 Å converted with the legacy cboan, as persson_cu111_potential.cpp does.
    "H2Cu111Persson": SystemParameters(8.0, "C6v", angstrom_to_bohr(2.5173), angstrom_to_bohr(2.5173)),
}
# Potentials whose parameters come from a CRP data set (crp_id, scm_id).
DATA_SET_SYSTEMS = ("CRPlibrary", "CRPDistrtSurface")


def is_known_system(system_id: str) -> bool:
    """True for every SystemID the program implements (potential_is_implemented)."""
    return any(seq(system_id, known) for known in (*ANALYTIC_SYSTEMS, *DATA_SET_SYSTEMS))


def is_data_set_system(system_id: str) -> bool:
    """True if Z_asymptotic and the point group come from a CRP data set."""
    return any(seq(system_id, known) for known in DATA_SET_SYSTEMS)


def system_parameters(
    system_id: str, z_asymptotic: float | None = None, point_group: str | None = None
) -> SystemParameters | None:
    """The parameters of a SystemID; the given values replace the built-in ones.

    Returns None when a value is unknown: an unknown SystemID, or a CRP data set without both values.
    """
    built_in = next((p for name, p in ANALYTIC_SYSTEMS.items() if seq(system_id, name)), None)
    if built_in is None and not is_data_set_system(system_id):
        return None
    z = z_asymptotic if z_asymptotic is not None else (built_in.z_asymptotic if built_in else None)
    group = point_group if point_group is not None else (built_in.point_group if built_in else None)
    if z is None or group is None:
        return None
    lattice = (built_in.lattice_a, built_in.lattice_b) if built_in else (None, None)
    return SystemParameters(z, group, *lattice)

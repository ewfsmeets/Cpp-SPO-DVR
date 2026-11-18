"""The system parameters the designer needs from each potential: Z_asymptotic and the point group.

Shadows: src/potential/potential_factory.cpp (potential_is_implemented),
src/potential/hard_wall_potential.cpp (z_asymptotic, HardWallPotential),
src/potential/dai_light_cu111_potential.cpp (z_asymptotic, DaiLightCu111Potential),
src/potential/persson_cu111_potential.cpp (z_asymptotic, PerssonCu111Potential)

The three analytic potentials have fixed values. A CRP data set (SystemID CRPlibrary, also under the SCM
correction CRPDistrtSurface) carries its own: the user gives them (the program's set-up output prints both).
"""

from dataclasses import dataclass

from program.real_compare import seq


@dataclass(frozen=True)
class SystemParameters:
    """What the designer uses of a potential's SystemParameters."""

    z_asymptotic: float  # bohr: beyond it V is the gas-phase V(r); the 6D potential is stored below it
    point_group: str  # C2v, C4v, C6v, C3v, C3 or none


# The potentials with built-in parameters (hard_wall_id, dai_light_id, persson_id of potential_factory.cpp).
ANALYTIC_SYSTEMS = {
    "HardWallPES": SystemParameters(7.0, "C2v"),
    "H2Cu111DaiLight": SystemParameters(15.0, "C6v"),
    "H2Cu111Persson": SystemParameters(8.0, "C6v"),
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
    return SystemParameters(z, group)

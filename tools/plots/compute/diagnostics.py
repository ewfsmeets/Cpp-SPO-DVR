"""The run verdict's findings: what the run's files say about its completeness and convergence. Each finding
has a severity (error, warning, information), its value and limit, what it means and the section of
inputs/README.md that explains it. Every threshold is a named constant here, marked "derived" (from the
physics or inputs/README.md) or "heuristic" (a judgement; change it if a case needs another); the texts are
formatted from the constants. NaN never passes a check: every test is written so that NaN raises the finding.
"""

from dataclasses import dataclass

import numpy as np

TURNING_POINT_BOHR = 2.0  # heuristic: Z_turn, "~2 bohr" in inputs/README.md section 8 (the designer's default)
UNANALYSED_LIMIT = 1e-3  # heuristic: norm still inside the analysis lines at tEnd
MOVING_WINDOW = 0.1  # heuristic: P_scattered is compared with its value this fraction of the run before tEnd ...
MOVING_RELATIVE = 1e-3  # heuristic: ... a change above this fraction of P_scattered is still moving ...
MOVING_FLOOR = 1e-4  # heuristic: ... and above this absolute change
UNITARITY_LIMIT = 1e-4  # heuristic: a P_reaction below -1e-4 is a unitarity error, not rounding
SLOW_CHANNELS_LIMIT = 1e-3  # heuristic: probability already seen in channels too slow to have all arrived

SECTION_TIME = "8. Total time and energy grid"
SECTION_STEPS = "4. dZ: the grid must hold the momenta"  # also the time step dtStep
SEVERITY_ORDER = {"error": 0, "warning": 1, "information": 2}


@dataclass
class Finding:
    """One finding of the verdict."""

    severity: str  # "error", "warning" or "information"
    name: str
    value: str
    limit: str  # "" when the finding has no limit
    meaning: str
    section: str  # the heading of inputs/README.md that explains it ("" for none)


def sorted_findings(findings: list[Finding | None]) -> list[Finding]:
    """The findings that apply, worst first."""
    return sorted((f for f in findings if f is not None), key=lambda f: SEVERITY_ORDER[f.severity])


def stopped_early(missing: list[str], last_bk_time_au: float | None, t_end_au: float) -> Finding | None:
    """Error when the files the program writes after the propagation are missing: the run stopped early, or is
    still running (the BK file, appended every analysis step, shows how far it got)."""
    if not missing:
        return None
    reached = "no BK record" if last_bk_time_au is None else f"BK records up to t = {last_bk_time_au:g} a.u."
    return Finding(
        "error",
        "run incomplete",
        f"{', '.join(missing)} missing; {reached}",
        f"tEnd = {t_end_au:g} a.u.",
        "the run stopped before tEnd or is still running: its results are not written yet",
        SECTION_TIME,
    )


def diverged(first_bad_time_au: float | None) -> Finding | None:
    """Error when the averages or the BK coefficients hold a NaN or inf."""
    if first_bad_time_au is None:
        return None
    return Finding(
        "error",
        "propagation diverged",
        f"NaN or inf from t = {first_bad_time_au:g} a.u.",
        "finite values",
        "every result after that time is meaningless (time step or grid too coarse for the potential?)",
        SECTION_STEPS,
    )


def not_checked(what: str, reason: str) -> Finding:
    """Information: a check that the run's files do not allow."""
    return Finding("information", f"not checked: {what}", reason, "", "", "")


ROWS_PER_CHUNK = 64  # (Z, r) points of the scattering wave squared at a time: memory stays small for any grid


def scattering_density_per_zr(scattering: np.ndarray) -> np.ndarray:
    """|psi|^2 summed over the (x, y) points and rotational states, per izr = ir + n_r iz (scattering:
    (n_zr, n_xy, n_rot), e.g. the memory-mapped wave of checkpoint.bin, read ROWS_PER_CHUNK rows at a time)."""
    density = np.empty(scattering.shape[0])
    for start in range(0, scattering.shape[0], ROWS_PER_CHUNK):
        rows = np.asarray(scattering[start : start + ROWS_PER_CHUNK])
        density[start : start + len(rows)] = np.sum(np.abs(rows) ** 2, axis=(1, 2))
    return density


def norm_inside_analysis_lines(
    scattering_density: np.ndarray,
    specular_density: np.ndarray,
    iz_analysis: int,
    surface_index: int,
    n_r: int,
    *,
    dz_bohr: float,
    dr_bohr: float,
    xy_volume_element_bohr2: float,
) -> float:
    """The norm at tEnd not yet past the analysis lines (an upper bound on the flux not yet analysed): the
    scattering wave at Z < Zana plus the specular wave (the initial channel, analysed at its own surface) at
    every Z up to and including that surface. The program's rest norm (src/propagation/averages.cpp) is the
    part of the latter between Zana and the surface. Both grids share dZ.

    scattering_density: |psi|^2 per (izr = ir + n_r iz) summed over (x, y) points and rotational states;
    specular_density: |chi|^2 per (iz, ir). Returns: a probability (dimensionless).
    """
    scattering = float(np.sum(scattering_density[: iz_analysis * n_r])) * dz_bohr * dr_bohr * xy_volume_element_bohr2
    specular = float(np.sum(specular_density[: surface_index + 1])) * dz_bohr * dr_bohr
    return scattering + specular


def unanalysed_flux(norm_inside: float) -> Finding | None:
    """Warning when more than UNANALYSED_LIMIT of the norm has not crossed the analysis lines at tEnd."""
    if norm_inside <= UNANALYSED_LIMIT:
        return None
    return Finding(
        "warning",
        "unanalysed flux, upper bound",
        f"{norm_inside:.2e}",
        f"<= {UNANALYSED_LIMIT:g}",
        "the norm still inside the analysis lines at tEnd: an upper bound on the flux not yet analysed, summed "
        "over all energies; at a single energy the error of P can be much larger (P_scattered still moving is "
        "the check per energy)",
        SECTION_TIME,
    )


def lowest_converged_normal_energy_ev(
    elapsed_time_au: float,
    *,
    z0_bohr: float,
    z_analysis_bohr: float,
    slowest_normal_energy_ev: float,
    total_mass_me: float,
    hartree_to_ev: float,
    z_turn_bohr: float = TURNING_POINT_BOHR,
) -> float:
    """E_n,min (eV), the lowest outgoing normal energy whose flux has crossed Zana by the elapsed time T:
    T = (Z0 - Z_turn) M / k_min + (Zana - Z_turn) M / k_n,min solved for k_n,min (inputs/README.md section 8),
    with k_min the momentum of the slowest incident normal energy. inf when that part of the packet has not
    reached Z_turn by T."""
    k_min = np.sqrt(2.0 * total_mass_me * slowest_normal_energy_ev / hartree_to_ev)
    time_left_au = elapsed_time_au - (z0_bohr - z_turn_bohr) * total_mass_me / k_min
    if time_left_au <= 0.0:
        return float("inf")
    k_n = (z_analysis_bohr - z_turn_bohr) * total_mass_me / time_left_au
    return float(k_n * k_n / (2.0 * total_mass_me) * hartree_to_ev)


def slow_probability(
    energy_index: np.ndarray,
    outgoing_normal_ev: np.ndarray,
    probability: np.ndarray,
    n_energies: int,
    e_n_min_ev: float,
) -> np.ndarray:
    """Per analysis energy (n_energies,): the probability in the channels that leave with a normal energy below
    E_n,min. The arguments are per row of s_matrix.tsv: the index of its analysis energy, its channel's
    outgoing normal energy E_total - E_channel (eV) and its P."""
    slow = outgoing_normal_ev < e_n_min_ev
    return np.bincount(energy_index[slow], weights=probability[slow], minlength=n_energies)


def slowest_converged(e_n_min_ev: float, slow_p: np.ndarray, collision_energy_ev: np.ndarray) -> Finding:
    """E_n,min and the largest probability the channels slower than it already carry: a warning above
    SLOW_CHANNELS_LIMIT (or NaN), else information."""
    worst = int(np.argmax(np.where(np.isnan(slow_p), np.inf, slow_p)))  # a NaN is the worst
    largest = float(slow_p[worst])
    return Finding(
        "information" if largest <= SLOW_CHANNELS_LIMIT else "warning",
        "slowest converged energy",
        f"E_n,min {e_n_min_ev:.3g} eV; the channels leaving slower carry up to {largest:.2e} "
        f"(E_coll {collision_energy_ev[worst]:.3g} eV)",
        f"<= {SLOW_CHANNELS_LIMIT:g} in those channels",
        "flux leaving with a normal energy below E_n,min has not all crossed the analysis line by tEnd: those "
        "channels' probabilities are too low and P_reaction too high (the initial channel, analysed at the "
        "specular line further out, is slower still)",
        SECTION_TIME,
    )


def norm_left(norm: float) -> Finding:
    """Information: the total norm at tEnd."""
    return Finding(
        "information",
        "norm left",
        f"{norm:.2e}",
        "",
        "it includes flux already analysed but not yet absorbed, and bounds only the error of 1 - sum R averaged "
        "over the packet's spectrum: at a single energy the error can be much larger",
        SECTION_TIME,
    )


def moving_window_start(times_au: np.ndarray) -> int:
    """The index of the record closest to (1 - MOVING_WINDOW) of the last time."""
    return int(np.argmin(np.abs(times_au - (1.0 - MOVING_WINDOW) * times_au[-1])))


def still_moving(p_end: np.ndarray, p_earlier: np.ndarray, energies_ev: np.ndarray) -> Finding | None:
    """Warning when P_scattered changed over the last MOVING_WINDOW of the run by more than
    max(MOVING_RELATIVE P_s, MOVING_FLOOR) at some analysis energy (every analysis energy lies within the
    packet's energy range, where its weight is at least exp(-2))."""
    change = np.abs(p_end - p_earlier)
    moving = ~(change <= np.maximum(MOVING_RELATIVE * p_end, MOVING_FLOOR))  # NaN counts as moving
    if not moving.any():
        return None
    return Finding(
        "warning",
        "P_scattered still moving",
        f"{int(moving.sum())} energies, up to {np.max(change[moving]):.2e}",
        f"change <= max({MOVING_RELATIVE:g} P_s, {MOVING_FLOOR:g}) over the last {MOVING_WINDOW:.0%} of the run",
        f"may not be converged at E_coll {energies_ev[moving].min():.3g}-{energies_ev[moving].max():.3g} eV: "
        "flux still arriving (a longer run), or flux reflected by an absorber that is too short (inputs/README.md "
        "section 7); the energy resolution is about 2 pi / tEnd",
        SECTION_TIME,
    )


def unitarity(error: float) -> Finding | None:
    """Warning when P_reaction falls below -UNITARITY_LIMIT."""
    if error <= UNITARITY_LIMIT:
        return None
    return Finding(
        "warning",
        "unitarity error",
        f">= {error:.2e}",
        f"<= {UNITARITY_LIMIT:g}",
        "P_reaction < 0: at the worst energy sum |S|^2 exceeds 1 by this much (time sums cut off before "
        "convergence, flux reflected by an absorber, or a numerical error)",
        SECTION_TIME,
    )

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
EDGE_LIMIT = 1e-3  # heuristic: a basis edge (top j', outermost diffraction channels) holding more of P_scattered
R_EDGE_LIMIT = 1e-3  # heuristic: |psi_v,j| at an end of the r grid above this fraction of its maximum
FRAME_EDGE_LIMIT = 1e-6  # heuristic: probability at a basis edge in a frame, relative to the initial norm
BUDGET_LIMIT = 1e-4  # heuristic: probability that may leave the grids beyond the absorbers' losses and their error
# indicator, relative to the initial norm (the scale of the unitarity error's limit): it covers the time stepping's
# own error of the absorbers (second order in dtStep: 4.1e-7 by t = 6000 a.u. on the hard-wall fixture at dtStep 2,
# 1.05e-7 at dtStep 1, ROpt eating the bound state's tail on the specular grid)

SECTION_TIME = "8. Total time and energy grid"
SECTION_STEPS = "4. dZ: the grid must hold the momenta"  # also the time step dtStep
SECTION_VERIFY = "10. Verifying a new input"
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
        "over all energies; at a single energy the error of P is first order in the amplitude, about "
        "2 sqrt(norm) (P_scattered still moving is the check per energy)",
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
        "it includes flux already analysed but not yet absorbed; the norm not yet analysed bounds only the loss "
        "averaged over the packet's spectrum: at a single energy the error of P is about 2 sqrt(N), and "
        "absorber reflection (about 2 sqrt(R) P) shows in no norm",
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
    ranking = np.where(moving, change, -1.0)  # only moving energies; a NaN change ranks first
    worst = int(np.argmax(np.where(np.isnan(ranking), np.inf, ranking)))
    return Finding(
        "warning",
        "P_scattered still moving",
        f"{int(moving.sum())} of {len(change)} energies, up to {change[worst]:.2e} at E_coll {energies_ev[worst]:.3g} eV",
        f"change <= max({MOVING_RELATIVE:g} P_s, {MOVING_FLOOR:g}) over the last {MOVING_WINDOW:.0%} of the run",
        f"may not be converged at E_coll {energies_ev[moving].min():.3g}-{energies_ev[moving].max():.3g} eV: "
        "flux still arriving (a longer run), or flux reflected by an absorber that is too short (inputs/README.md "
        "section 7); the energy resolution is about 2 pi / tEnd",
        SECTION_TIME,
    )


def unitarity(error: float, worst_energy_ev: float) -> Finding | None:
    """Warning when P_reaction falls below -UNITARITY_LIMIT (worst_energy_ev: where it is lowest)."""
    if error <= UNITARITY_LIMIT:
        return None
    return Finding(
        "warning",
        "unitarity error",
        f">= {error:.2e} (P_reaction lowest at E_coll {worst_energy_ev:.3g} eV)",
        f"<= {UNITARITY_LIMIT:g}",
        "P_reaction < 0: at the worst energy sum |S|^2 exceeds 1 by this much (time sums cut off before "
        "convergence, flux reflected by an absorber, or a numerical error)",
        SECTION_TIME,
    )


def basis_edge(what: str, largest_fraction: float, energy_ev: float, meaning: str) -> Finding | None:
    """Warning when a basis edge holds more than EDGE_LIMIT of P_scattered at some energy (or NaN)."""
    if largest_fraction <= EDGE_LIMIT:
        return None
    return Finding(
        "warning",
        f"flux at the edge of the basis: {what}",
        f"up to {largest_fraction:.2e} of P_scattered (E_coll {energy_ev:.3g} eV)",
        f"<= {EDGE_LIMIT:g}",
        meaning,
        SECTION_VERIFY,
    )


def frame_edge(
    severity: str, what: str, shares: np.ndarray, times_au: np.ndarray, meaning: str, section: str = SECTION_VERIFY
) -> Finding | None:
    """A finding of the given severity when, in some frame, a basis edge holds more than FRAME_EDGE_LIMIT of the
    initial norm (shares: per frame; NaN counts as above), citing the given section of inputs/README.md."""
    worst = int(np.argmax(np.nan_to_num(shares, nan=np.inf)))
    if shares[worst] <= FRAME_EDGE_LIMIT:
        return None
    value = f"up to {shares[worst]:.2e} of the initial norm (t = {times_au[worst]:g} a.u.)"
    return Finding(
        severity, f"packet at the edge of the basis: {what}", value, f"<= {FRAME_EDGE_LIMIT:g}", meaning, section
    )


def rotational_edge(what: str, shares: np.ndarray, times_au: np.ndarray) -> Finding | None:
    """Information when an edge of the rotational basis (the top J or |mJ|) holds more than FRAME_EDGE_LIMIT of the
    initial norm in some frame: while the molecule is at the surface J is a basis label and such a passing
    population is normal; the S matrix's edge check (EDGE_LIMIT) decides whether the basis is large enough."""
    meaning = (
        f"populated while the packet passes: not an error by itself; the check of j' in the S matrix (above "
        f"{EDGE_LIMIT:g} of P_scattered) tells whether the basis is large enough"
    )
    return frame_edge("information", what, shares, times_au, meaning)


def grid_edge(what: str, rises: np.ndarray, times_au: np.ndarray) -> Finding | None:
    """Warning when a grid edge (the first or last r points, the first Z points) rises by more than FRAME_EDGE_LIMIT
    of the initial norm above its first frame's share (rises: compute.frames.grid_edge_rises): the packet reaches the
    edge, where the FFT wraps it around to the other end."""
    meaning = "the packet reaches the edge of the grid, where the FFT wraps it around to the other end"
    return frame_edge("warning", f"{what} (above the first frame's)", rises, times_au, meaning)


def momentum_edge(what: str, rises: np.ndarray, times_au: np.ndarray) -> Finding | None:
    """Warning when, in some frame, the probability at the edge of a momentum grid rises by more than
    FRAME_EDGE_LIMIT of the initial norm above the first frame's (rises: compute.frame_extras.momentum_edges; |k|
    near the grid's largest momentum, or the outer diffraction orders): the FFT cannot hold larger momenta and folds
    them over. Cites inputs/README.md section 4 (the grid must hold the momenta)."""
    meaning = "momenta reach the edge of the grid: the FFT folds larger ones over to the other side (aliasing)"
    return frame_edge(
        "warning", f"{what} (momentum grid, above the first frame's)", rises, times_au, meaning, SECTION_STEPS
    )


def top_rotational_level(fraction: np.ndarray, energies_ev: np.ndarray, j_top: int) -> Finding | None:
    """P(j' = the top j' of the rotational basis) / P_scattered per energy: warning above EDGE_LIMIT."""
    worst = int(np.argmax(np.nan_to_num(fraction, nan=np.inf)))
    meaning = f"j' = {j_top} is the highest j of the rotational basis (JMax): a larger JMax may change the results"
    return basis_edge(f"j' = {j_top}", float(fraction[worst]), energies_ev[worst], meaning)


def outer_diffraction_channels(fraction: np.ndarray, energies_ev: np.ndarray) -> Finding | None:
    """P in the outermost analysed diffraction channels (n or m at an end of its range, (0, 0) excepted) /
    P_scattered per energy: warning above EDGE_LIMIT."""
    worst = int(np.argmax(np.nan_to_num(fraction, nan=np.inf)))
    meaning = "the outermost analysed diffraction channels carry flux: a larger nX, nY may change the results"
    return basis_edge("outermost diffraction channels", float(fraction[worst]), energies_ev[worst], meaning)


def r_grid_ends(largest_ratio: float, level: tuple[int, int], inner: bool | None) -> Finding | None:
    """Warning when |psi_v,j(r)| of an analysed level at an end of the r grid exceeds R_EDGE_LIMIT of its
    maximum: the grid cuts the gas-phase states. inner: the larger end is the first r point (RStart too large;
    on the periodic grid that tail also reappears at the last point), False the last (the grid too short), None
    unknown (a NaN state: both advices)."""
    if largest_ratio <= R_EDGE_LIMIT:
        return None
    ends = {
        True: ("the first r point", "a smaller RStart"),
        False: ("the last r point", "a longer r grid"),
        None: ("an end of the r grid", "a smaller RStart or a longer r grid"),
    }
    end, advice = ends[inner]
    return Finding(
        "warning",
        "gas-phase states reach the ends of the r grid",
        f"|psi| at {end} {largest_ratio:.2e} of its maximum (v = {level[0]}, j = {level[1]})",
        f"<= {R_EDGE_LIMIT:g}",
        "the r grid cuts the gas-phase states: their energies shift, and the r absorber removes part of the "
        f"initial state from the start (counted as reaction); use {advice}",
        SECTION_VERIFY,
    )


def norm_by_region(time_au: float, regions: dict[str, float], specular_inside_line: float) -> Finding:
    """Information: where the norm left is at the last densities/ record (the regions of
    compute.densities.region_norms: r >= r_c, and the Z bands for r < r_c), and the specular wave between Zana
    and its analysis line (the rest norm, compute.densities.specular_inside_line)."""
    listed = "; ".join(f"{label}: {value:.1e}" for label, value in regions.items())
    return Finding(
        "information",
        "norm by region at the end",
        f"t = {time_au:g} a.u.: {listed}; of the bands beyond Zana, the specular wave inside its analysis line: "
        f"{specular_inside_line:.1e}",
        "",
        "the bands below Zana are not yet analysed; the bands beyond Zana are analysed, except the specular wave "
        "inside its own analysis line (listed last: the rest norm of averages.dat), and wait to be absorbed; "
        "r >= r_c is dissociating, counted as reaction",
        SECTION_TIME,
    )


def budget_open(*, excess: float, residual: float, error_estimate: float, t_au: float) -> Finding | None:
    """Warning when the norm budget does not close: at some window end the residual (the norm that left the grids
    minus the absorbers' summed losses) exceeds its error estimate by more than BUDGET_LIMIT of the initial norm
    (excess > 0, a NaN counts; compute.norm_budget.budget_check, at its worst window)."""
    if excess <= 0.0:
        return None
    return Finding(
        "warning",
        "norm budget does not close",
        f"|residual| {abs(residual):.2e} at t = {t_au:g} a.u., error estimate {error_estimate:.2e}",
        f"|residual| <= error estimate + {BUDGET_LIMIT:g} of the initial norm",
        "the absorbers' losses do not account for the norm that left the grids: probability left another way, most "
        "likely through the angular transform while the packet is at the surface (it keeps only the rotational "
        "basis's part at every time step; less with a smaller dtStep), or the run diverged; it counts as reaction in "
        "P_reaction = 1 - P_scattered",
        SECTION_VERIFY,
    )

"""Quantities of the FrameExtras datasets of frames/ (readers/frames.py): the alignment of the molecular axis per
steradian (angular), the probability per parallel-momentum order (diffraction) and per gas-phase level (vibration)
against Z, <P2(cos theta)> and the rotational kinetic energy against Z (rotation_moments), the momentum distributions
in r and in Z (momentum), and how much of the wave packet reaches the edges of the momentum grids.

The functions on frames take records before any NaN or inf (compute.frames.usable_frames); probabilities are per grid
cell or per momentum point, as the program writes them.
"""

from dataclasses import dataclass

import numpy as np

from compute.frames import EdgeSeries
from compute.gas_phase import Levels
from readers.frames import Frames

# Heuristic: the outer 20 % of a momentum grid (|k| above this fraction of its largest, pi / d) counts as its edge.
# A packet designed by inputs/README.md section 4 keeps its momenta well inside; the FFT folds momenta beyond the
# grid over to the other side, and what approaches the edge is the warning sign. The same band of orders in (x, y).
NYQUIST_FRACTION = 0.8
MAX_ORDERS = 6  # the diffraction panels draw the six largest orders, the rest summed: six line styles stay legible
ISOTROPIC_PER_SR = 1.0 / (4.0 * np.pi)  # an isotropic distribution of the axis, per steradian


def probability_z(frames: Frames, record: int) -> np.ndarray:
    """(nz,): the probability at each Z point of the scattering grid in one record (pop_jmj_z summed over the
    rotational states): what the extras' Z profiles sum to."""
    return np.asarray(frames.series.arrays["pop_jmj_z"][record]).sum(axis=1)


def densest_z(frames: Frames, record: int, stop: int) -> int:
    """The Z index below `stop` with the most probability in one record."""
    return int(np.argmax(probability_z(frames, record)[:stop]))


@dataclass(frozen=True, kw_only=True)
class StericDistribution:
    """The distribution of the molecular axis at one Z point of one record: probability per steradian, normalised to
    1 over the sphere (isotropic: ISOTROPIC_PER_SR), at the DVR nodes (phi from the first lattice vector). For a
    homonuclear molecule the nodes cover theta <= 90 degrees (an alignment, not an orientation): the other
    hemisphere is the mirror (pi - theta, phi + pi)."""

    theta_deg: np.ndarray  # (ntheta,)
    phi_deg: np.ndarray  # (nphi,)
    per_sr: np.ndarray  # (ntheta, nphi)
    z_bohr: float
    t_au: float
    homonuclear: bool


def steric_distribution_at(frames: Frames, record: int, iz: int) -> StericDistribution:
    """The distribution at Z index iz: prob_angle_z divided by each node's solid angle w_t 2 pi / nphi and by the
    probability at that Z."""
    nodes = frames.extras.angular
    probability = np.asarray(frames.series.arrays["prob_angle_z"][record, iz])  # (ntheta, nphi)
    solid_angle_sr = nodes.theta_weights[:, None] * (2.0 * np.pi / len(nodes.phi_rad))
    return StericDistribution(
        theta_deg=np.degrees(nodes.theta_rad),
        phi_deg=np.degrees(nodes.phi_rad),
        per_sr=probability / solid_angle_sr / probability.sum(),
        z_bohr=float(frames.z_bohr[iz]),
        t_au=float(frames.series.times_au[record]),
        homonuclear=nodes.homonuclear,
    )


def fft_orders(n: int) -> np.ndarray:
    """(n,): the order at each FFT index as the program labels it (src/analysis/analysis_basis.cpp): k for k <= n / 2,
    else k - n (the Nyquist index n / 2 is the order +n / 2)."""
    k = np.arange(n)
    return np.where(k <= n // 2, k, k - n)


@dataclass(frozen=True, kw_only=True)
class ZProfiles:
    """Probability against Z of one record split into labelled parts that sum to the probability at each Z."""

    labels: list[str]
    values: np.ndarray  # (n_parts, nz)
    t_au: float


def diffraction_profiles(frames: Frames, records: list[int]) -> list[ZProfiles]:
    """P(n, m; Z) of each record: the MAX_ORDERS orders (n, m) with the most probability in all the records together,
    largest first (the same orders, in the same order, for every record), and the other orders summed."""
    pop = [np.asarray(frames.series.arrays["pop_nm_z"][k]).reshape(frames.nz_scattering, -1).T for k in records]
    ranking = np.argsort(-sum(p.sum(axis=1) for p in pop), kind="stable")
    top, rest = ranking[:MAX_ORDERS], ranking[MAX_ORDERS:]
    n, m = np.meshgrid(fft_orders(frames.nx), fft_orders(frames.ny), indexing="ij")
    labels = [f"(n, m) = ({n.flat[k]}, {m.flat[k]})" for k in top] + (["other orders"] if len(rest) else [])
    profiles = []
    for k, per_order in zip(records, pop, strict=True):
        values = [per_order[i] for i in top] + ([per_order[rest].sum(axis=0)] if len(rest) else [])
        profiles.append(ZProfiles(labels=labels, values=np.array(values), t_au=float(frames.series.times_au[k])))
    return profiles


def bound_levels(levels: Levels, j_values: np.ndarray, nr: int) -> np.ndarray:
    """(nr, nJ) bool in pop_vj_z's order (v, then the J of j_values): the level lies below V_asy at the last r point
    (compute.gas_phase.levels; the r absorber and the centrifugal term there are ignored). The FGH gives nr levels
    per J; ValueError when gas_phase_states.dat lacks one."""
    bound = np.zeros((nr, len(j_values)), dtype=bool)
    found = np.zeros((nr, len(j_values)), dtype=bool)
    for v, j, is_bound in zip(levels.v, levels.j, levels.bound, strict=True):
        if v < nr and j in j_values:
            j_index = int(np.searchsorted(j_values, j))
            bound[v, j_index], found[v, j_index] = is_bound, True
    if not found.all():
        raise ValueError("gas_phase_states.dat lacks a level of the frames' vibration extra")
    return bound


def vibration_profiles(frames: Frames, record: int, bound: np.ndarray, n_analysed: int) -> ZProfiles:
    """P(v; Z) of the vibrational levels v < n_analysed (summed over J), then the other levels below V_asy(r_end) and
    those above it (bound: bound_levels; a level v < n_analysed above V_asy(r_end) stays in its v)."""
    pop = np.asarray(frames.series.arrays["pop_vj_z"][record])  # (nz, nr, nJ)
    analysed = np.arange(pop.shape[1]) < n_analysed
    labels = [f"v = {v}" for v in range(n_analysed)]
    values = [pop[:, v, :].sum(axis=1) for v in range(n_analysed)]
    labels += ["other levels below V_asy(r_end)", "levels above V_asy(r_end)"]
    values.append((pop * (bound & ~analysed[:, None])).sum(axis=(1, 2)))
    values.append((pop * (~bound & ~analysed[:, None])).sum(axis=(1, 2)))
    return ZProfiles(labels=labels, values=np.array(values), t_au=float(frames.series.times_au[record]))


@dataclass(frozen=True, kw_only=True)
class RotationMoments:
    """<P2(cos theta)> and the rotational kinetic energy <J(J+1) / (2 mu r^2)> at each Z of one record (NaN where
    nothing is), with the probability at each Z (the panels mask where it is small)."""

    p2: np.ndarray  # (nz,)
    rotational_energy_ev: np.ndarray  # (nz,)
    probability: np.ndarray  # (nz,)
    t_au: float


def ratio_or_nan(numerator: np.ndarray, probability: np.ndarray) -> np.ndarray:
    """numerator / probability, NaN where the probability is zero."""
    return np.divide(numerator, probability, out=np.full_like(numerator, np.nan), where=probability > 0)


def rotation_moments(frames: Frames, record: int, hartree_to_ev: float) -> RotationMoments:
    """The moments: the program's numerators divided by the probability at each Z."""
    probability = probability_z(frames, record)
    p2_numerator = np.asarray(frames.series.arrays["p2_numerator_z"][record])
    rot_numerator_ev = np.asarray(frames.series.arrays["rot_numerator_z"][record]) * hartree_to_ev
    return RotationMoments(
        p2=ratio_or_nan(p2_numerator, probability),
        rotational_energy_ev=ratio_or_nan(rot_numerator_ev, probability),
        probability=probability,
        t_au=float(frames.series.times_au[record]),
    )


def momenta_per_bohr(n: int, spacing_bohr: float) -> np.ndarray:
    """(n,): the momentum at each FFT index of a grid of n points spaced spacing_bohr apart, bohr^-1."""
    return 2.0 * np.pi * np.fft.fftfreq(n, spacing_bohr)


def edge_per_bohr(momenta: np.ndarray) -> float:
    """Where a momentum grid's edge begins: NYQUIST_FRACTION of its largest momentum (pi / d for the even sizes the
    program uses)."""
    return NYQUIST_FRACTION * float(np.max(np.abs(momenta)))


def probability_near_edge(probability: np.ndarray, momenta: np.ndarray) -> np.ndarray:
    """The probability at |k| > edge_per_bohr (summed over the last axis of probability: the momentum points)."""
    return probability[..., np.abs(momenta) > edge_per_bohr(momenta)].sum(axis=-1)


@dataclass(frozen=True, kw_only=True)
class MomentumDistribution:
    """A probability per momentum point of one record, sorted by momentum, and where the grid's edge begins."""

    k_per_bohr: np.ndarray  # (n,), ascending
    probability: np.ndarray  # (n,)
    edge_per_bohr: float
    t_au: float


def momentum_distribution(frames: Frames, record: int, dataset: str) -> MomentumDistribution:
    """prob_kr (on the r grid) or prob_kz_spec (on the specular Z grid) of one record, sorted by momentum."""
    grid = frames.r_bohr if dataset == "prob_kr" else frames.z_bohr
    k = momenta_per_bohr(len(grid), grid[1] - grid[0])
    order = np.argsort(k, kind="stable")
    return MomentumDistribution(
        k_per_bohr=k[order],
        probability=np.asarray(frames.series.arrays[dataset][record])[order],
        edge_per_bohr=edge_per_bohr(k),
        t_au=float(frames.series.times_au[record]),
    )


def momentum_edges(frames: Frames, reference: float) -> list[EdgeSeries]:
    """Per record, the probability that the run brings to the edges of its momentum grids, relative to the initial
    norm (reference): the rise above the first record's share (as for the grid edges; an initial state that is not
    smooth on its grid, e.g. a gas-phase state reaching the r grid's ends, has a share there of its own) at |k| above
    edge_per_bohr in r and in Z of chi (momentum), and in the orders with |n| or |m| above NYQUIST_FRACTION of the
    largest, nX / 2 or nY / 2 (diffraction)."""
    arrays = frames.series.arrays
    shares = []
    if "momentum" in frames.extras.names:
        k_r = momenta_per_bohr(len(frames.r_bohr), frames.r_bohr[1] - frames.r_bohr[0])
        k_z = momenta_per_bohr(len(frames.z_bohr), frames.z_bohr[1] - frames.z_bohr[0])
        shares.append(("k_r", probability_near_edge(np.asarray(arrays["prob_kr"]), k_r)))
        shares.append(("k_Z of the specular wave", probability_near_edge(np.asarray(arrays["prob_kz_spec"]), k_z)))
    if "diffraction" in frames.extras.names:
        pop = np.asarray(arrays["pop_nm_z"]).sum(axis=1)  # (records, nx, ny)
        n, m = np.meshgrid(fft_orders(frames.nx), fft_orders(frames.ny), indexing="ij")
        edge = (np.abs(n) >= NYQUIST_FRACTION * (frames.nx // 2)) | (np.abs(m) >= NYQUIST_FRACTION * (frames.ny // 2))
        shares.append(("outer diffraction orders", pop[:, edge].sum(axis=1)))
    return [EdgeSeries(label=label, values=(share - share[0]) / reference) for label, share in shares]

"""compute/frame_extras.py on synthetic frames/: the distribution per steradian, the order and level profiles over Z,
the rotational moments and the momentum edges, each against a hand count, with a slip that a test must catch."""

from pathlib import Path

import numpy as np
import pytest

from compute.frame_extras import (
    MAX_ORDERS,
    NYQUIST_FRACTION,
    bound_levels,
    densest_z,
    diffraction_profiles,
    fft_orders,
    momentum_distribution,
    momentum_edges,
    probability_near_edge,
    rotation_moments,
    steric_distribution_at,
    vibration_profiles,
)
from compute.gas_phase import Levels
from readers.arrays import Series
from readers.frames import AngularNodes, FrameExtrasInfo, Frames

NZ, NR, NX, NY = 3, 4, 4, 2
NZ_SPECULAR = 6  # an even size, as every FFT size of the program


def synthetic_frames(arrays: dict[str, np.ndarray], names: tuple[str, ...]) -> Frames:
    """Two records on NZ x NR x NX x NY points and the states (J, mJ) = (0, 0), (2, 0); the given extras' datasets;
    the angular extra on 2 theta nodes (weights 0.5, 1.5) x 4 phi nodes; vibration with J = 0, 2. The probability per
    Z is 1, 3, 0 in record 0 and 1, 1, 4 in record 1."""
    pop = np.array([[[1.0, 0.0], [2.0, 1.0], [0.0, 0.0]], [[0.5, 0.5], [1.0, 0.0], [3.0, 1.0]]])  # (2, NZ, 2)
    series = Series(
        directory=Path("frames"),
        manifest={},
        times_au=np.array([0.0, 100.0]),
        arrays={"pop_jmj_z": pop, **arrays},
        first_non_finite=None,
    )
    nodes = AngularNodes(
        theta_rad=np.radians([30.0, 90.0]),
        theta_weights=np.array([0.5, 1.5]),
        phi_rad=np.radians([0.0, 90.0, 180.0, 270.0]),
        homonuclear=True,
    )
    return Frames(
        series=series,
        z_bohr=0.5 * np.arange(NZ_SPECULAR),
        r_bohr=0.5 + 0.25 * np.arange(NR),
        nz_scattering=NZ,
        nx=NX,
        ny=NY,
        states=np.array([[0, 0], [2, 0]]),
        initial_index=0,
        sites=[],
        starts_au=[0.0],
        extras=FrameExtrasInfo(names=names, angular=nodes, vibration_j_values=np.array([0, 2])),
    )


def test_densest_z_and_the_fft_orders() -> None:
    """The densest Z of record 1 is index 2, below index 2 it is index 0 (the first of equals); of record 0 index 1.
    The FFT orders as the program labels them: 0, 1, 2, -1 for 4 points (the Nyquist index is +2), 0, 1 for 2."""
    frames = synthetic_frames({}, ())
    assert densest_z(frames, 1, NZ) == 2 and densest_z(frames, 1, 2) == 0 and densest_z(frames, 0, NZ) == 1
    assert list(fft_orders(4)) == [0, 1, 2, -1] and list(fft_orders(2)) == [0, 1]


def test_steric_distribution_is_normalised_over_the_sphere() -> None:
    """The probabilities 1 everywhere but 3 at (theta 30, phi 90), total 10: per steradian prob / (w_t pi / 2) / 10, so
    that node's value is 3 / (0.5 pi / 2) / 10 and the sum of per_sr x the solid angle is 1. Dividing by w_t alone
    (no dphi, a slip) breaks the normalisation."""
    angle = np.ones((2, NZ, 2, 4))
    angle[0, 1, 0, 1] = 3.0
    frames = synthetic_frames({"prob_angle_z": angle}, ("angular",))
    steric = steric_distribution_at(frames, 0, 1)
    solid_angle = np.array([0.5, 1.5])[:, None] * (np.pi / 2)
    assert abs(steric.per_sr[0, 1] - 3.0 / (0.5 * np.pi / 2) / 10.0) < 1e-15
    assert abs(np.sum(steric.per_sr * solid_angle) - 1.0) < 1e-15
    assert list(steric.theta_deg.round(12)) == [30.0, 90.0] and steric.z_bohr == 0.5 and steric.homonuclear
    slip = angle[0, 1] / np.array([0.5, 1.5])[:, None] / 10.0
    assert abs(np.sum(slip * solid_angle) - 1.0) > 0.1


def test_diffraction_profiles_rank_the_orders_once() -> None:
    """Record 0 holds k + 1 at the flat (ix, iy) index k (twice that at its last Z point: the ranking sums over Z),
    record 1 three times 8 - k: summed over Z and both records index k holds 4 (k + 1) + 9 (8 - k) = 76 - 5 k, so
    index 0, (n, m) = (0, 0), comes first in both records' panels: the same orders in the same order. MAX_ORDERS of
    them, the rest summed; each Z sums to the record's total. Ranking record 0 alone (a slip) puts index 7, (n, m) =
    (-1, 1), first."""
    orders = np.arange(1.0, NX * NY + 1).reshape(NX, NY)
    pop = np.stack([np.tile(orders, (NZ, 1, 1)), 3.0 * np.tile(orders[::-1, ::-1], (NZ, 1, 1))])
    pop[0, 2] *= 2.0
    profiles = diffraction_profiles(synthetic_frames({"pop_nm_z": pop}, ("diffraction",)), [0, 1])
    assert profiles[0].labels == profiles[1].labels and profiles[0].labels[0] == "(n, m) = (0, 0)"
    assert len(profiles[0].labels) == MAX_ORDERS + 1 and profiles[0].labels[-1] == "other orders"
    for record, profile in enumerate(profiles):
        assert np.allclose(profile.values.sum(axis=0), pop[record].sum(axis=(1, 2)), rtol=0, atol=1e-13)
    alone = diffraction_profiles(synthetic_frames({"pop_nm_z": pop}, ("diffraction",)), [0])[0]
    assert alone.labels[0] == "(n, m) = (-1, 1)"


def test_vibration_profiles_split_below_and_above() -> None:
    """Levels v = 0..3 of J = 0 and 2; below V_asy(r_end): v < 2 for J = 0, none for J = 2 (so v = 0 of J = 2, an
    analysed level, is above it and stays in v = 0). With 1 analysed level: v = 0 (1 + 10), the other levels below
    (J = 0 v = 1: 2) and the rest (3 + 4 + 20 + 30 + 40). Each Z sums to the total (110). Without the analysed mask on
    the levels above (a slip), v = 0 of J = 2 would also count there. A missing level is a ValueError."""
    pop = np.zeros((2, NZ, NR, 2))
    pop[0, :, :, 0] = [1.0, 2.0, 3.0, 4.0]
    pop[0, :, :, 1] = [10.0, 20.0, 30.0, 40.0]
    v = np.repeat(np.arange(NR), 2)
    j = np.tile([0, 2], NR)
    below = (j == 0) & (v < 2)
    levels = Levels(v=v, j=j, energy_ev=np.arange(2.0 * NR), bound=below)
    bound = bound_levels(levels, np.array([0, 2]), NR)
    frames = synthetic_frames({"pop_vj_z": pop}, ("vibration",))
    profiles = vibration_profiles(frames, 0, bound, 1)
    assert profiles.labels == ["v = 0", "other levels below V_asy(r_end)", "levels above V_asy(r_end)"]
    assert np.allclose(profiles.values[:, 0], [11.0, 2.0, 97.0], rtol=0, atol=1e-13)
    assert np.allclose(profiles.values.sum(axis=0), 110.0, rtol=0, atol=1e-13)
    assert (pop[0, 0] * ~bound).sum() != 97.0  # the slip
    partial = Levels(v=v[:-1], j=j[:-1], energy_ev=np.arange(2.0 * NR - 1), bound=below[:-1])
    with pytest.raises(ValueError, match="lacks a level"):
        bound_levels(partial, np.array([0, 2]), NR)


def test_rotation_moments_are_ratios() -> None:
    """Record 0's probability per Z is 1, 3, 0: <P2> = numerator / probability, NaN where nothing is; the rotational
    energy in eV with the given factor."""
    numerators = {
        "p2_numerator_z": np.array([[0.5, 1.5, 0.0]] * 2),
        "rot_numerator_z": np.array([[2.0, 3.0, 0.0]] * 2),
    }
    moments = rotation_moments(synthetic_frames(numerators, ("rotation_moments",)), 0, 10.0)
    assert np.allclose(moments.p2[:2], [0.5, 0.5], rtol=0, atol=1e-15) and np.isnan(moments.p2[2])
    assert np.allclose(moments.rotational_energy_ev[:2], [20.0, 10.0], rtol=0, atol=1e-13)


def test_momentum_edges_by_hand() -> None:
    """probability_near_edge on 8 momenta (largest pi): only |k| > NYQUIST_FRACTION pi counts (3 pi / 4 does not, the
    Nyquist point does). The edges, the rise above the first record over the reference 2: k_r (4 r points: only k =
    -pi / dr), k_Z of chi (6 Z points: k = +-2 pi / 3 dz are inside 0.8 pi / dz, -pi / dz is not), and the orders
    with |n| >= 0.8 x 2 or |m| >= 0.8 x 1, i.e. n = 2 or m = 1. Without the first record subtracted (a slip) k_r would
    be 0.5, not 0. The momentum distribution comes sorted, with its edge."""
    k = 2 * np.pi * np.fft.fftfreq(8, 1.0)
    assert 0.75 < NYQUIST_FRACTION < 1.0
    assert probability_near_edge(np.arange(8.0)[None, :], k)[0] == 4.0  # index 4: k = -pi
    pop = np.zeros((2, NZ, NX, NY))
    pop[1, 0, 2, 0] = 1.0  # n = 2: an outer order
    pop[1, 1, 1, 1] = 3.0  # m = 1: an outer order
    pop[1, 2, 1, 0] = 5.0  # (1, 0): inside
    k_z_spec = np.array([np.zeros(NZ_SPECULAR), np.arange(1.0, NZ_SPECULAR + 1)])
    arrays = {"prob_kr": np.ones((2, NR)), "prob_kz_spec": k_z_spec, "pop_nm_z": pop}
    frames = synthetic_frames(arrays, ("diffraction", "momentum"))
    edges = momentum_edges(frames, 2.0)
    assert [e.label for e in edges] == ["k_r", "k_Z of the specular wave", "outer diffraction orders"]
    assert list(edges[0].values) == [0.0, 0.0]
    assert list(edges[1].values) == [0.0, 4.0 / 2.0]  # index 3 of 6: k = -pi / dz
    assert list(edges[2].values) == [0.0, 2.0]
    distribution = momentum_distribution(frames, 1, "prob_kz_spec")
    assert list(distribution.probability) == [4.0, 5.0, 6.0, 1.0, 2.0, 3.0]  # FFT order 0, 1, 2, -3, -2, -1 sorted
    assert abs(distribution.edge_per_bohr - NYQUIST_FRACTION * np.pi / 0.5) < 1e-12

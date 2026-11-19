"""The final-state quantities on hand-made values: the alignment A0(2) and the sums of s_matrix.tsv rows."""

import numpy as np

from compute.states import alignment_a0, helicopter_limit, level_map, mj_fractions, state_to_state
from readers.tables import Table


def test_alignment_of_hand_computed_distributions() -> None:
    """j = 1: only mJ = 0 gives A0 = (0 - 2) / 2 = -1 ("cartwheel"); only |mJ| = 1 gives (3 - 2) / 2 = 1/2
    ("helicopter"); equal populations give 0. j = 2, P(mJ) = (1, 0, 2, 0, 1): (2*(12-6) + 2*(0-6)) / (6*4) = 0."""
    assert alignment_a0(np.array([0.0, 1.0, 0.0]), 1) == -1.0
    assert alignment_a0(np.array([0.5, 0.0, 0.5]), 1) == 0.5
    assert abs(alignment_a0(np.array([1.0, 1.0, 1.0]), 1)) < 1e-15
    assert abs(alignment_a0(np.array([1.0, 0.0, 2.0, 0.0, 1.0]), 2)) < 1e-15
    assert abs(alignment_a0(np.array([0.0, 1.0, 0.0, 1.0, 0.0]), 2) - (3 - 6) / 6) < 1e-15


def test_j_squared_instead_of_j_j_plus_1_fails() -> None:
    """Negative test: j^2 for j(j+1) (a realistic slip) gives (3 - 1) / 1 = 2, not 1/2, for |mJ| = 1 of j = 1."""
    m_j, p, j = np.array([-1, 0, 1]), np.array([0.5, 0.0, 0.5]), 1
    wrong = np.sum(p * (3 * m_j * m_j - j * j)) / (j * j * np.sum(p))
    assert wrong != alignment_a0(p, j)


def table(rows: list[tuple]) -> Table:
    """An s_matrix.tsv table of (E_coll, v, j, mJ, n, m, P) rows."""
    columns = np.array(rows, dtype=float).T
    return Table(dict(zip(("E_coll", "v", "j", "mJ", "n", "m", "P"), columns, strict=True)), None)


def test_sums_over_channels_and_truncated_basis() -> None:
    """P(v', j') sums mJ' and the diffraction channels; with MJMax = 1 the mJ' = +-2 of j' = 2 are not in the
    basis (NaN) and its alignment is left out; with MJMax = 2 it is defined."""
    rows = [(0.3, 0, 0, 0, 0, 0, 0.5), (0.3, 0, 2, 1, 0, 0, 0.1), (0.3, 0, 2, -1, 1, 0, 0.2), (0.4, 1, 0, 0, 0, 0, 0.7)]
    truncated = state_to_state(table(rows), np.array([0.3, 0.4]), mj_max=1, j_top=2)
    assert truncated.levels == [(0, 0), (0, 2), (1, 0)]
    assert np.allclose(truncated.p_vj, [[0.5, 0.3, 0.0], [0.0, 0.0, 0.7]], rtol=0, atol=1e-15)
    assert np.array_equal(truncated.p_mj[2][0], [np.nan, 0.2, 0.0, 0.1, np.nan], equal_nan=True)
    assert np.isnan(truncated.alignment[2]).all()
    full = state_to_state(table(rows), np.array([0.3, 0.4]), mj_max=2, j_top=2)
    assert list(full.p_mj[2][0]) == [0.0, 0.2, 0.0, 0.1, 0.0]
    assert abs(full.alignment[2][0] - alignment_a0(full.p_mj[2][0], 2)) == 0 and np.isnan(full.alignment[2][1])


def test_rounding_noise_is_left_out() -> None:
    """A level, and an alignment, with probability below LEVEL_FLOOR of P_scattered are rounding noise: the
    level is not listed and the alignment of that energy is NaN."""
    rows = [(0.3, 0, 0, 0, 0, 0, 1.0), (0.3, 0, 2, 1, 0, 0, 1e-26)]
    result = state_to_state(table(rows), np.array([0.3]), mj_max=2, j_top=2)
    assert result.levels == [(0, 0)] and np.isnan(result.alignment[2][0])
    assert mj_fractions(result, 0) == {}


def test_level_map_and_mj_fractions() -> None:
    """The map holds P(v', j') at its (v', j') cell, NaN elsewhere; each j''s mJ' fractions sum to 1; the
    helicopter limit of j' = 2 is 1."""
    rows = [(0.3, 0, 0, 0, 0, 0, 0.5), (0.3, 0, 2, 1, 0, 0, 0.1), (0.3, 0, 2, -2, 0, 0, 0.3), (0.3, 1, 0, 0, 0, 0, 0.1)]
    result = state_to_state(table(rows), np.array([0.3]), mj_max=2, j_top=2)
    v_values, j_values, grid = level_map(result, 0)
    assert (v_values, j_values) == ([0, 1], [0, 2])
    assert np.allclose(grid, [[0.5, 0.4], [0.1, np.nan]], equal_nan=True, rtol=0, atol=1e-15)
    fractions = mj_fractions(result, 0)
    assert list(fractions) == [2] and abs(np.sum(fractions[2]) - 1.0) < 1e-15
    assert helicopter_limit(2) == 1.0

"""The absorbers' factors per time step at the grid points, from run_info.json alone: the independent input of the
numpy checks that need them (test_runs_norm_budget, test_runs_frame_extras). Not a test module."""

import numpy as np

from readers.run_info import RunInfo


def absorbers_on(info: RunInfo, coordinate: str, grid: str) -> list[dict]:
    """The absorbers of run_info.json on `coordinate` ("Z" or "r") that act on `grid` (ROpt's grid is "both")."""
    return [a for a in info.data["absorbers"]["list"] if a["coordinate"] == coordinate and a["grid"] in (grid, "both")]


def absorber_factors(info: RunInfo, absorbers: list[dict], points: np.ndarray) -> np.ndarray:
    """exp(-V_opt dt) per time step at the grid points (1 outside), from each absorber's start, end, strength and the
    power p: V_opt = strength (p + 1) / 2 x^p, x rising from 0 at the start of a high-end absorber to 1 at its end
    (falling for a low-end one, *Low; a later absorber overwrites an earlier one where they meet). The formula of
    src/propagation/propagation_tables.cpp set_optical_factors (Fortran: Absorption.f90 SetupAbsorpPot)."""
    power, dt = info.data["absorbers"]["power"], info.data["times"]["dt_step"]
    factors = np.ones(len(points))
    for a in absorbers:
        inside = (points >= a["start"]) & (points <= a["end"])
        length = a["end"] - a["start"]
        x = (a["end"] - points) / length if a["name"].endswith("Low") else (points - a["start"]) / length
        factors = np.where(inside, np.exp(-(a["strength"] * 0.5 * (power + 1) * x**power) * dt), factors)
    return factors

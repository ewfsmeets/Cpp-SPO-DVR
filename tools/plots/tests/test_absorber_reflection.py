"""R and T of a Z absorber (compute/absorber_reflection.py): a free wave passes, each energy is computed on its
own, the wave-packet designer's scalar loop agrees (run in a process of its own: its packages are not the plot
tool's), and the numbers of inputs/README.md section 7 for the DZ absorbers and of the hard-wall fixture."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from compute import absorber_reflection as a

DESIGNER = Path(__file__).resolve().parents[2] / "wave_packet_designer"
EV = 27.2114  # the program's legacy hartree -> eV
H2_TOTAL_MASS_ME = 3674.217084
H2_REDUCED_MASS_ME = 918.554271
# (name, length_bohr, strength_hartree, mass_me): the DZ case's absorbers and the hard-wall fixture's ZSpOpt
ABSORBERS = [
    ("DZ ZOpt", 4.08, 1.27096140e-02, H2_TOTAL_MASS_ME),
    ("DZ ZSpOpt", 8.52, 1.01936078e-02, H2_TOTAL_MASS_ME),
    ("DZ ROpt", 4.2, 2.45596439e-02, H2_REDUCED_MASS_ME),
    ("hardwall_tiny ZSpOpt", 3.52, 0.01829896432306322, H2_TOTAL_MASS_ME),
]
DESIGNER_SCRIPT = """
import json, sys
sys.path.insert(0, sys.argv[1])
from design.transmission import reflection_transmission
cases = json.load(sys.stdin)
print(json.dumps([[reflection_transmission(e, c["length"], c["strength"], c["mass"]) for e in c["energies"]]
                  for c in cases]))
"""


def response(
    energies_ev: list[float] | np.ndarray, length_bohr: float, strength_hartree: float, mass_me: float
) -> tuple[np.ndarray, np.ndarray]:
    """R and T at energies in eV."""
    energies_hartree = np.asarray(energies_ev, dtype=float) / EV
    return a.reflection_transmission(
        energies_hartree, length_bohr=length_bohr, strength_hartree=strength_hartree, mass_me=mass_me
    )


def test_free_wave_passes_unchanged() -> None:
    """With strength 0 the transfer matrix gives R = 0 and T = 1 (to the rounding of its 2000 slices)."""
    reflection, transmission = response([0.3], 4.0, 0.0, 3674.2)
    assert reflection[0] < 1e-24 and abs(transmission[0] - 1.0) < 1e-12


def test_each_energy_on_its_own() -> None:
    """The vectorised result of each energy equals, bitwise, that energy computed alone."""
    energies_ev = np.geomspace(0.001, 3.0, 25)
    for _, length_bohr, strength_hartree, mass_me in ABSORBERS:
        together = response(energies_ev, length_bohr, strength_hartree, mass_me)
        for i, energy_ev in enumerate(energies_ev):
            alone = response([energy_ev], length_bohr, strength_hartree, mass_me)
            assert together[0][i] == alone[0][0] and together[1][i] == alone[1][0]


def designer_responses(cases: list[dict]) -> list[list[list[float]]]:
    """The designer's reflection_transmission (design/transmission.py) per case and energy, as [R, T]: run in a
    Python process of its own with the designer's directory first on sys.path."""
    completed = subprocess.run(
        [sys.executable, "-I", "-B", "-c", DESIGNER_SCRIPT, str(DESIGNER)],
        input=json.dumps(cases),
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, f"the designer's process exited with {completed.returncode}: {completed.stderr}"
    return json.loads(completed.stdout)


def test_agrees_with_the_designer() -> None:
    """R and T agree with the designer's scalar loop within 1e-8 of R + T (numpy's complex arithmetic rounds
    differently from cmath's) for the DZ absorbers and the hard-wall fixture's ZSpOpt, 1 meV to 3 eV."""
    energies_hartree = np.geomspace(0.001, 3.0, 25) / EV
    cases = [
        {"length": length, "strength": strength, "mass": mass, "energies": energies_hartree.tolist()}
        for _, length, strength, mass in ABSORBERS
    ]
    responses = designer_responses(cases)
    for (name, length_bohr, strength_hartree, mass_me), expected in zip(ABSORBERS, responses, strict=True):
        reflection, transmission = a.reflection_transmission(
            energies_hartree, length_bohr=length_bohr, strength_hartree=strength_hartree, mass_me=mass_me
        )
        designer = np.array(expected)
        scale = designer[:, 0] + designer[:, 1]
        assert np.all(np.abs(reflection - designer[:, 0]) <= 1e-8 * scale), name
        assert np.all(np.abs(transmission - designer[:, 1]) <= 1e-8 * scale), name


def test_dz_absorbers_as_the_readme_says() -> None:
    """Section 7, for the DZ absorbers (as the designer's tests/test_design.py): scattering R + T <= 4e-5 for
    0.05-1.0 eV and 8e-5 at 1.5 eV; specular R + T <= 4e-9 for 0.2-1.0 eV (one printed digit)."""
    _, z_length_bohr, z_strength_hartree, mass_me = ABSORBERS[0]
    _, zsp_length_bohr, zsp_strength_hartree, _ = ABSORBERS[1]
    scattering = sum(response([0.05, 0.1, 0.2, 0.4, 0.7, 1.0], z_length_bohr, z_strength_hartree, mass_me))
    at_1_5_ev = sum(response([1.5], z_length_bohr, z_strength_hartree, mass_me))
    specular = sum(response([0.2, 0.4, 0.7, 1.0], zsp_length_bohr, zsp_strength_hartree, mass_me))
    assert (f"{scattering.max():.0e}", f"{at_1_5_ev[0]:.0e}", f"{specular.max():.0e}") == ("4e-05", "8e-05", "4e-09")


def test_hard_wall_fixture_specular_absorber() -> None:
    """The hard-wall fixture's 3.52 bohr ZSpOpt reflects R = 4.86e-7 at 0.2 eV, above the section 7 target
    (R + T <= 4e-9): the verdict's "absorber reflection" warning on that fixture."""
    _, length_bohr, strength_hartree, mass_me = ABSORBERS[3]
    reflection, _ = response([0.2], length_bohr, strength_hartree, mass_me)
    assert f"{reflection[0]:.2e}" == "4.86e-07"


def test_rescaling_keeps_r_and_t(monkeypatch: pytest.MonkeyPatch) -> None:
    """psi and psi' rescaled at every slice (a low RESCALE_ABOVE) give R and T within 1e-10 of R + T of the
    unscaled loop; a long strong absorber, where the unscaled wave would overflow, gives finite R and T."""
    energies_ev = np.geomspace(0.01, 2.0, 7)
    reference = response(energies_ev, 4.08, 1.27096140e-02, H2_TOTAL_MASS_ME)
    monkeypatch.setattr(a, "RESCALE_ABOVE", 1.0)
    rescaled = response(energies_ev, 4.08, 1.27096140e-02, H2_TOTAL_MASS_ME)
    scale = reference[0] + reference[1]
    for before, after in zip(reference, rescaled, strict=True):
        assert np.all(np.abs(after - before) <= 1e-10 * scale)  # measured: 9.3e-12 (2040 divisions)
    monkeypatch.undo()
    reflection, transmission = response([0.2], 200.0, 1.0, H2_TOTAL_MASS_ME)
    assert np.isfinite(reflection[0]) and 0.0 <= reflection[0] <= 1.0 and 0.0 <= transmission[0] < 1e-100


def test_returned_amplitude_per_row() -> None:
    """a = sqrt(R) + T per row, computed once per distinct energy and mapped back to the rows."""
    energies_hartree = np.array([0.3, 0.1, 0.3, 0.2]) / EV
    reflection, amplitude = a.returned_amplitude(
        energies_hartree, length_bohr=3.52, strength_hartree=0.01829896432306322, mass_me=H2_TOTAL_MASS_ME
    )
    r, t = a.reflection_transmission(
        energies_hartree, length_bohr=3.52, strength_hartree=0.01829896432306322, mass_me=H2_TOTAL_MASS_ME
    )
    assert np.array_equal(reflection, r) and np.array_equal(amplitude, np.sqrt(r) + t)
    assert reflection[0] == reflection[2] and amplitude[0] == amplitude[2]

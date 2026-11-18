"""The gas-phase potential V(r) of the DZ input, as the program printed it: a stand-in and test data.

The designer uses it for a potential whose gas phase it does not copy (program/gas_phase_potentials.py: the
Persson potential, the CRP data sets), as an indicative stand-in; and test case `gas-phase` checks the copied
Dai-Light formula against it bitwise.

Data layer (README, "Architecture"). Source: the block GAS-PHASE POTENTIAL (lines 142-188) of the committed
DZ run output docs/validation_thesis_fig5_3_dz/data/logs/cleaned_cuda_strict.txt (2026-11-08), that is, the Dai-Light H2/Cu(111) potential (SystemID H2Cu111DaiLight) at Z = Z_asymptotic = 15 bohr and
theta = pi/2, on that input's r grid (RStart 0.5, dR 0.2, nR 40), in eV as printed (17 digits). It is the
Born-Oppenheimer curve of the free H2 molecule, so it serves D2 and HD too. For another potential the states
it gives are only indicative; the report says so.
"""

from program.gas_phase import GasPhaseCurve
from program.unit_conversion import ev_to_hartree

SYSTEM_ID = "H2Cu111DaiLight"
SOURCE = "the gas-phase curve of the DZ input (Dai-Light H2/Cu(111) potential)"
R_START = 0.5  # bohr
DR = 0.2  # bohr
POTENTIAL_EV = (  # V(r) at r = R_START + i * DR, i = 0 .. 39
    6.2529805913068461e00,  # r = 0.5
    5.2248575994617608e-01,  # r = 0.7
    -2.6028191166022210e00,  # r = 0.9
    -4.1263668968755214e00,  # r = 1.1
    -4.6900378476700180e00,  # r = 1.3
    -4.7003762586442397e00,  # r = 1.5
    -4.4114507038166062e00,  # r = 1.7
    -3.9791245809655948e00,  # r = 1.9
    -3.4964973223126630e00,  # r = 2.1
    -3.0169634172775024e00,  # r = 2.3
    -2.5691443802100493e00,  # r = 2.5
    -2.1664992302563184e00,  # r = 2.7
    -1.8134609808937108e00,  # r = 2.9
    -1.5093138405066528e00,  # r = 3.1
    -1.2506082173052870e00,  # r = 3.3
    -1.0326353095663074e00,  # r = 3.5
    -8.5030179225637448e-01,  # r = 3.7
    -6.9862595048924414e-01,  # r = 3.9
    -5.7299843418226382e-01,  # r = 4.1
    -4.6929964842292554e-01,  # r = 4.3
    -3.8393241946409906e-01,  # r = 4.5
    -3.1380689092758207e-01,  # r = 4.7
    -2.5630058648760901e-01,  # r = 4.9
    -2.0920757297925505e-01,  # r = 5.1
    -1.7068492656390874e-01,  # r = 5.3
    -1.3920109763353219e-01,  # r = 5.5
    -1.1348853563229702e-01,  # r = 5.7
    -9.2501581141926553e-02,  # r = 5.9
    -7.5379838840396174e-02,  # r = 6.1
    -6.1416803205046548e-02,  # r = 6.3
    -5.0033284732241709e-02,  # r = 6.5
    -4.0755091841170814e-02,  # r = 6.7
    -3.3194407790419676e-02,  # r = 6.9
    -2.7034328517344171e-02,  # r = 7.1
    -2.2016075119155026e-02,  # r = 7.3
    -1.7928450963799115e-02,  # r = 7.5
    -1.4599170728089798e-02,  # r = 7.7
    -1.1887742940500930e-02,  # r = 7.9
    -9.6796368501893132e-03,  # r = 8.1
    -7.8815078848688738e-03,  # r = 8.3
)


def dz_gas_phase_curve() -> GasPhaseCurve:
    """The curve in atomic units, converted with the program's factor (as the program prints it: in eV)."""
    return GasPhaseCurve(R_START, DR, tuple(ev_to_hartree(v) for v in POTENTIAL_EV), SYSTEM_ID, SOURCE)

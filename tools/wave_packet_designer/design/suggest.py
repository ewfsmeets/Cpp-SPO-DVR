"""Suggested Z grids, packet position and Z absorbers: the recipe of inputs/README.md, section 2, steps 3-11.

The energy range, the masses, ZStart and Z_asymptotic come from the input; the rest follows:
- dZ: the input's, if it holds the momenta ((k_grid - k0)/sigma_k >= 6, section 4); else the largest step
  that does, rounded down to DZ_STEP. Section 4 has two more limits that the designer does not check: the
  local momenta in wells or the exit channel, and the time step;
- Zana: the first grid point >= Z_asymptotic;
- Z0: Zana + 8 sigma_Z, rounded up to Z0_STEP;
- ZSpOptStart: the first grid point >= Z0 + 5.2 sigma_Z + 2 bohr;
- nZ: the smallest valid FFT size whose last point lies SCATTERING_ABSORBER_LENGTH beyond Zana; the
  scattering absorber runs from Zana to that point;
- nZSpecular: the largest valid FFT size whose specular absorber (to its last point) stays within the
  strength table at EMinWP, L/lambda <= 15 ("as long as the table allows"; more points cost nothing);
- the target energies 0.2 eV (scattering) and EMinWP (specular), where the input tunes them by energy;
  absorbers tuned by a factor keep it;
- tEnd: long enough for outgoing normal energies down to E_n,min (0.03 eV: section 8, step 11), in whole averaging
  intervals dtAve, so nStep, dtAna and dtAve stay consistent.
Lengths are rounded to 10 decimals for the eye where the rounded value equals the exact one within the
program's tolerance (1e-12); otherwise the exact value is kept.
"""

import math

from design.rules import LOWEST_OUTGOING_ENERGY, TURNING_POINT, t_end_for
from program.config import DerivedParameters, InputParameters
from program.fft_sizes import is_valid_fft_size, nearest_fft_size
from program.real_compare import req, rge
from program.unit_conversion import PI, ev_to_hartree
from program.validate import steps_per_analysis

DZ_STEP = 0.01  # bohr: a suggested dZ is rounded down to a multiple of this
Z0_STEP = 0.1  # bohr: Z0 is rounded up to a multiple of this
SCATTERING_ABSORBER_LENGTH = 4.0  # bohr, section 7
TAIL_GAP = 2.0  # bohr between the packet's outer tail and ZSpOptStart, section 5
LARGEST_REDUCED_LENGTH = 15.0  # the last L/lambda of the strength table
SCATTERING_TARGET = 0.2  # eV, section 7
# The InputParameters fields a suggestion may set (the energies only where the input tunes by energy).
SUGGESTED_FIELDS = (
    "dz",
    "nz",
    "nz_specular",
    "z_analysis",
    "z0",
    "z_opt_start",
    "z_opt_end",
    "zsp_opt_start",
    "zsp_opt_end",
    "z_opt_energy",
    "zsp_opt_energy",
    "t_end",
    "n_steps",
)


def suggested_keys(
    p: InputParameters,
    d: DerivedParameters,
    z_asymptotic: float,
    dz: float | None = None,
    e_n_min: float = LOWEST_OUTGOING_ENERGY,
    z_turn: float = TURNING_POINT,
) -> dict[str, float | int]:
    """The suggested values of the Z keys and tEnd, as InputParameters field names; dZ = `dz` if given."""
    sigma_k = math.sqrt(d.wave_packet_alpha)
    sigma_z = 1.0 / (2.0 * sigma_k)
    largest_dz = PI / (d.kz + 6.0 * sigma_k)  # (k_grid - k0)/sigma_k = 6
    if dz is None:
        dz = p.dz if rge(largest_dz, p.dz) else round(math.floor(largest_dz / DZ_STEP) * DZ_STEP, 10)

    def point(i: int) -> float:
        return grid_point(p.z_start, dz, i)

    def first_index_at_or_above(z: float) -> int:
        return first_index_at_or_beyond(p.z_start, dz, z)

    i_analysis = analysis_index(p.z_start, dz, z_asymptotic)
    z_analysis = point(i_analysis)
    z0 = math.ceil((z_analysis + 8.0 * sigma_z) / Z0_STEP) * Z0_STEP
    i_specular_start = first_index_at_or_above(z0 + 5.2 * sigma_z + TAIL_GAP)

    nz = nearest_fft_size(first_index_at_or_above(z_analysis + SCATTERING_ABSORBER_LENGTH) + 1)
    wavelength = 2.0 * PI / math.sqrt(2.0 * d.total_mass * ev_to_hartree(p.e_min_wave_packet))
    longest = LARGEST_REDUCED_LENGTH * wavelength
    nz_specular = largest_fft_size(math.floor((point(i_specular_start) + longest - p.z_start) / dz) + 1)
    nz_specular = max(nz_specular, nearest_fft_size(max(nz, i_specular_start + 2)))

    changes: dict[str, float | int] = {
        "dz": dz,
        "nz": nz,
        "nz_specular": nz_specular,
        "z_analysis": tidy(z_analysis),
        "z0": tidy(z0),
        "z_opt_start": tidy(z_analysis),
        "z_opt_end": tidy(point(nz - 1)),
        "zsp_opt_start": tidy(point(i_specular_start)),
        "zsp_opt_end": tidy(point(nz_specular - 1)),
    }
    if p.z_opt_energy is not None:
        changes["z_opt_energy"] = SCATTERING_TARGET
    if p.zsp_opt_energy is not None:
        changes["zsp_opt_energy"] = p.e_min_wave_packet
    k_min = d.kz - 2.0 * sigma_k  # the slow end of the packet
    needed = t_end_for(changes["z0"], changes["z_analysis"], d.total_mass, k_min, e_n_min, z_turn)
    changes.update(fit_time(p, needed))
    return changes


def fit_time(p: InputParameters, t_end: float) -> dict[str, float | int]:
    """tEnd and nStep for a propagation of at least `t_end` (a.u.): whole averaging intervals dtAve after
    tStart, so tEnd = tStart + nStep dtStep and the propagation time is a whole multiple of dtAna and dtAve
    (the program's checks), in whole multiples of the program's truncated time steps per analysis too. For a valid
    dtAna that is dtAna / dtStep itself; a decimal ratio such as 1.2 / 0.1 truncates to 11, and the program refuses
    such a dtAna (docs/fortran_deviations.md, analysis D5), which the fitted nStep leaves as the only time error."""
    steps_per_average = math.lcm(round(p.dt_averages / p.time_step), steps_per_analysis(p))
    intervals = max(1, math.ceil((t_end - p.t_start) / (steps_per_average * p.time_step) - 1e-9))
    n_steps = intervals * steps_per_average
    return {"t_end": p.t_start + n_steps * p.time_step, "n_steps": n_steps}


def fit_z_grids(
    p: InputParameters, z_asymptotic: float, scattering_length: float, specular_length: float
) -> dict[str, float | int]:
    """The keys that follow from dZ, ZSpOptStart and the two absorber lengths (bohr), as the recipe places
    them (the window's sliders): Zana on the first grid point >= Z_asymptotic, the scattering absorber from
    Zana and the specular one from ZSpOptStart, each ending on the last point of the smallest valid FFT grid
    that holds it."""
    dz = p.dz
    i_analysis = analysis_index(p.z_start, dz, z_asymptotic)
    z_analysis = grid_point(p.z_start, dz, i_analysis)
    nz = nearest_fft_size(first_index_at_or_beyond(p.z_start, dz, z_analysis + scattering_length) + 1)
    specular_end = first_index_at_or_beyond(p.z_start, dz, p.zsp_opt_start + specular_length)
    nz_specular = nearest_fft_size(max(nz, specular_end + 1))
    return {
        "nz": nz,
        "nz_specular": nz_specular,
        "z_analysis": tidy(z_analysis),
        "z_opt_start": tidy(z_analysis),
        "z_opt_end": tidy(grid_point(p.z_start, dz, nz - 1)),
        "zsp_opt_end": tidy(grid_point(p.z_start, dz, nz_specular - 1)),
    }


def fit_r_grid(p: InputParameters) -> dict[str, float]:
    """ROptEnd on the last point of the r grid (the window's dR and nR sliders): the r absorber ends there."""
    return {"r_opt_end": tidy(p.r_start + (p.nr - 1) * p.dr)}


def grid_point(z_start: float, dz: float, i: int) -> float:
    """Point i of a Z grid, as the program computes it."""
    return z_start + i * dz


def first_index_at_or_beyond(z_start: float, dz: float, z: float) -> int:
    """The first grid index whose point lies at or beyond z (the 1e-9 keeps a point that lies on z)."""
    return math.ceil((z - z_start) / dz - 1e-9)


def analysis_index(z_start: float, dz: float, z_asymptotic: float) -> int:
    """The index of the first grid point >= Z_asymptotic by the program's tolerant comparison (Zana)."""
    i = first_index_at_or_beyond(z_start, dz, z_asymptotic)
    while not rge(grid_point(z_start, dz, i), z_asymptotic):
        i += 1
    return i


def tidy(x: float) -> float:
    """x rounded to 10 decimals if that equals x within the program's tolerance (else x itself)."""
    rounded = round(x, 10)
    return rounded if req(rounded, x) else x


def largest_fft_size(n: int) -> int:
    """The largest valid FFT size <= n (at least 2, the smallest valid size)."""
    while n > 2 and not is_valid_fft_size(n):
        n -= 1
    return max(n, 2)

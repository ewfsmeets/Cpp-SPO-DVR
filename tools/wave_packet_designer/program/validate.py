"""The program's consistency checks of the input, in the order and with the comparisons of the C++.

Shadows: src/input/validate.cpp (validate_input)

Each `checks.error_if(violated, key, message)` line copies one `check(...)` line of validate_input, under
the same section headings and with the same key and message (tests/test_validate.py compares them). Not
copied: the flux-analysis block (not ported; a propagating run with doFluxAnalysis stops anyway, see
evaluation.py). One difference in a check: the write2DPES bounds x2DPES <= alat and y2DPES <= blat are skipped
when the lattice is unknown here (a CRP data set; with CRPDistrtSurface, the slab snapshot's lattice).
"""

import math

from program.config import InputParameters
from program.fft_sizes import is_valid_fft_size
from program.findings import Finding, FindingList
from program.potential_systems import SystemParameters, is_known_system
from program.real_compare import nint, req, rge, rgt, rle, rlt, rne, seq
from program.unit_conversion import PI, ev_to_hartree


def validate_input(p: InputParameters, system: SystemParameters, z_asymptotic: float) -> list[Finding]:
    """Every consistency error of the input (the program stops at the first one)."""
    checks = FindingList("validate_input")
    check_system_and_grids(p, system, checks)
    check_rotational_basis(p, checks)
    check_propagation_and_analysis(p, z_asymptotic, checks)
    check_wave_packet_and_state(p, z_asymptotic, checks)
    check_potential_and_absorbers(p, z_asymptotic, checks)
    check_pes_cut(p, system, checks)
    return checks.findings


def check_system_and_grids(p: InputParameters, system: SystemParameters, checks: FindingList) -> None:
    # --- system ID and degrees of freedom ---
    checks.error_if(not is_known_system(p.system_id), "SystemID", "Unknown system ID")
    checks.error_if(not p.z_dof, "ZDoF", "Currently only 6D calculations are supported.")
    checks.error_if(not p.r_dof, "RDoF", "Currently only 6D calculations are supported.")
    checks.error_if(not p.x_dof, "XDoF", "Currently only 6D calculations are supported.")
    checks.error_if(not p.y_dof, "YDoF", "Currently only 6D calculations are supported.")
    checks.error_if(not p.theta_dof, "ThetaDoF", "Currently only 6D calculations are supported.")
    checks.error_if(not p.phi_dof, "PhiDoF", "Currently only 6D calculations are supported.")

    # --- masses ---
    checks.error_if(rle(p.mass1, 0.0), "Mass1", "Mass of atom 1 <= 0.")
    checks.error_if(rle(p.mass2, 0.0), "Mass2", "Mass of atom 2 <= 0.")
    checks.error_if(
        p.homonuclear and rne(p.mass1, p.mass2), "HomoNuc", "masses are different, but molecule is homo-nuclear."
    )
    checks.error_if(
        not p.homonuclear and req(p.mass1, p.mass2), "HomoNuc", "masses are the same, but molecule is hetero-nuclear."
    )

    # --- scattering and specular grid in Z ---
    checks.error_if(p.nz < 1, "nZ", "nZ < 1.")
    checks.error_if(not is_valid_fft_size(p.nz), "nZ", "nZ is not a valid FFT size.")
    checks.error_if(rle(p.dz, 0.0), "dZ", "dZ <= 0.")
    checks.error_if(p.nz_specular < 1, "nZSpecular", "nZSpecular < 1.")
    checks.error_if(not is_valid_fft_size(p.nz_specular), "nZSpecular", "nZSpecular is not a valid FFT size.")

    # --- grid in r ---
    checks.error_if(rle(p.r_start, 0.0), "RStart", "RStart <= 0.")
    checks.error_if(p.nr < 1, "nR", "nR < 1.")
    checks.error_if(not is_valid_fft_size(p.nr), "nR", "nR is not a valid FFT size.")
    checks.error_if(rle(p.dr, 0.0), "dR", "dR <= 0.")

    # --- grid in x and y ---
    checks.error_if(p.nx < 1, "nX", "nX < 1.")
    checks.error_if(p.ny < 1, "nY", "nY < 1.")
    checks.error_if(not is_valid_fft_size(p.nx), "nX", "nX is not a valid FFT size.")
    checks.error_if(not is_valid_fft_size(p.ny), "nY", "nY is not a valid FFT size.")
    if any(seq(system.point_group, group) for group in ("C4v", "C6v", "C3v")):
        checks.error_if(p.x_dof and p.y_dof and p.nx != p.ny, "nY", "nX /= nY.")


def check_rotational_basis(p: InputParameters, checks: FindingList) -> None:
    # --- rotational basis ---
    checks.error_if(p.j_in < 0, "JIn", "JIn < 0.")
    checks.error_if(p.j_max < 0, "JMax", "JMax < 0.")
    checks.error_if(abs(p.mj_max) > p.j_max, "MJMax", "|MJMax| > JMax.")
    checks.error_if(abs(p.mj_in) > p.j_in, "MJIn", "|MJIn| > JIn.")
    checks.error_if(p.j_in >= p.j_max, "JIn", "JIn >= JMax.")  # upstream r1726 (was >)
    checks.error_if(abs(p.mj_in) > abs(p.mj_max), "MJIn", "|MJIn| > |MJMax|.")
    if p.homonuclear:
        checks.error_if(
            p.j_max % 2 != p.j_in % 2,
            "JMax",
            "the molecule is homo-nuclear but JMax is even and JIn is odd, or vice versa.",
        )
        checks.error_if(
            p.j_max % 2 != abs(p.mj_max) % 2,
            "MJMax",
            "the molecule is homo-nuclear, but JMax is even and |MJMax| is odd, or vice versa.",
        )

    # --- projection operator ---
    # Kept as in the C++, although the reader (json_input.py) already rejects 0 and 1: ProjOp = 1, merging the
    # specular wave into the scattering grid, is a legacy option of the Fortran and not ported.
    checks.error_if(p.projection_operator not in (0, 1, 2), "ProjOp", "ProjOp must be 0, 1 or 2.")


def steps_per_analysis(p: InputParameters) -> int:
    """Time steps per analysis, truncated as the Fortran INT (nanspo). Shadows: src/input/process.cpp."""
    return int(p.dt_analysis / p.time_step)


def check_propagation_and_analysis(p: InputParameters, z_asymptotic: float, checks: FindingList) -> None:
    # --- split operator propagation ---
    n_steps = float(p.n_steps)
    checks.error_if(rlt(p.t_start, 0.0), "tStart", "tStart < 0.")
    checks.error_if(p.new_calculation and rne(p.t_start, 0.0), "tStart", "This is a new calculation but tStart /= 0.")
    checks.error_if(not p.new_calculation and req(p.t_start, 0.0), "tStart", "This is a restart but tStart = 0.")
    checks.error_if(
        not p.new_calculation and p.restart_dir == "",
        "RestartDir",
        "This is a restart (NewCalc = false) but RestartDir is not given.",
    )
    checks.error_if(
        p.only_analysis and p.restart_dir == "",
        "RestartDir",
        "An analysis-only run (onlyAnalysis) needs RestartDir: the directory with the BK coefficient file.",
    )
    checks.error_if(
        p.restart_format not in ("native", "fortran"), "RestartFormat", 'RestartFormat must be "native" or "fortran".'
    )
    checks.error_if(
        p.checkpoint_format not in ("native", "fortran", "both"),
        "CheckpointFormat",
        'CheckpointFormat must be "native", "fortran" or "both".',
    )
    checks.error_if(rlt(p.t_end, 0.0), "tEnd", "tEnd < 0.")
    checks.error_if(rgt(p.t_start, p.t_end), "tStart", "tStart > tEnd.")
    checks.error_if(rne(p.t_end, p.t_start + n_steps * p.time_step), "tEnd", "tEnd /= tStart + nStep*dtStep.")
    checks.error_if(rle(p.time_step, 0.0), "dtStep", "dtStep <= 0.")
    checks.error_if(p.n_steps < 0, "nStep", "nStep < 0.")

    # --- analysis ---
    checks.error_if(rle(p.dt_analysis, 0.0), "dtAna", "dtAna <= 0.")
    checks.error_if(rlt(p.dt_analysis, p.time_step), "dtAna", "dtAna < dtStep.")
    # analysis D4 (as validate.cpp, with the ratios themselves), only where the C++ gets that far: it stops at the
    # first error, so after "dtAna < dtStep." it never reports the band below dtStep.
    analysis_ratio = (
        p.dt_analysis / p.time_step
        if not rle(p.time_step, 0.0) and not rle(p.dt_analysis, 0.0) and not rlt(p.dt_analysis, p.time_step)
        else None
    )
    analysis_steps = math.floor(analysis_ratio) if analysis_ratio is not None else None
    checks.error_if(
        analysis_ratio is not None and analysis_ratio < 1.0,
        "dtAna",
        f"dtAna = {p.dt_analysis} lies within the 1e-12 tolerance below dtStep = {p.time_step}, so dtAna / dtStep "
        "truncates to 0 time steps per analysis; use dtAna >= dtStep.",
    )
    checks.error_if(rle(p.dt_averages, 0.0), "dtAve", "dtAve <= 0.")
    checks.error_if(rlt(p.dt_averages, p.time_step), "dtAve", "dtAve < dtStep.")
    if p.j_analysis is not None:
        checks.error_if(p.j_analysis < 0, "JAna", "JAna < 0.")
    if p.v_analysis is not None:
        checks.error_if(p.v_analysis < 0, "VAna", "VAna < 0.")
    if p.diffraction_order_analysis is not None:
        checks.error_if(p.diffraction_order_analysis < 0, "DifOrdAna", "DifOrdAna < 0.")
    # NrOfEs < 2, as the Fortran since r1719 (docs/fortran_deviations.md, analysis D1).
    checks.error_if(
        p.n_energies < 2,
        "NrOfEs",
        f"NrOfEs = {p.n_energies}: the analysis needs at least 2 energies, (EMaxWP-EMinWP)/(NrOfEs-1) apart; use "
        "NrOfEs >= 2.",
    )
    if rle(p.dt_analysis, 0.0) or rle(p.dz, 0.0):
        return  # the divisions below would fail; the program stopped at the checks above already
    propagation_time = n_steps * p.time_step
    checks.error_if(
        rne(nint(propagation_time / p.dt_analysis) * p.dt_analysis, propagation_time),
        "dtAna",
        "Propagation time must be an integer times dtAna.",
    )
    checks.error_if(
        analysis_ratio is not None
        and analysis_ratio >= 1.0
        and (p.n_steps != 0 if math.floor(analysis_ratio) > p.n_steps else p.n_steps % math.floor(analysis_ratio) != 0),
        "dtAna",
        f"nStep = {p.n_steps} is not a multiple of the {analysis_steps} time steps per analysis that dtAna / dtStep = "
        f"{analysis_ratio} truncates to; the run would stop before tEnd. Use a dtStep and dtAna whose ratio is exact in "
        "binary, e.g. dtStep 0.125 or 0.25.",
    )
    checks.error_if(rlt(p.dt_averages, p.dt_analysis), "dtAve", "dtAve < dtAna.")
    checks.error_if(
        rne(nint(p.dt_averages / p.dt_analysis) * p.dt_analysis, p.dt_averages),
        "dtAve",
        "dtAve must be an integer times dtAna.",
    )
    average_ratio = p.dt_averages / p.dt_analysis if not rlt(p.dt_averages, p.dt_analysis) else None
    checks.error_if(
        average_ratio is not None and average_ratio < 1.0,
        "dtAve",
        f"dtAve = {p.dt_averages} lies within the 1e-12 tolerance below dtAna = {p.dt_analysis}, so dtAve / dtAna "
        "truncates to 0 analyses per average; use dtAve >= dtAna.",
    )

    checks.error_if(rlt(p.z_analysis, z_asymptotic), "Zana", "ZAna < ZAsymptotic.")
    checks.error_if(
        rne(nint((p.z_analysis - p.z_start) / p.dz) * p.dz, p.z_analysis - p.z_start),
        "Zana",
        "ZAna not on a grid point.",
    )
    checks.error_if(rlt(p.z_analysis, p.z_start), "Zana", "ZAna < ZStart.")
    checks.error_if(
        rgt(p.z_analysis, p.z_start + (p.nz - 1) * p.dz), "Zana", "ZAnalysis is too big for the scattering grid in Z."
    )
    checks.error_if(rgt(p.z_analysis, p.z_opt_start), "Zana", "ZAnalysis > ZOptStart.")

    if p.j_analysis is not None:
        checks.error_if(p.j_analysis > p.j_max, "JAna", "JAna > JMax.")
    if p.v_analysis is not None:
        checks.error_if(p.v_analysis > p.nr - 1, "VAna", "VAna > nR-1.")
    if p.diffraction_order_analysis is not None:
        checks.error_if(
            p.diffraction_order_analysis > p.nx + p.ny - 2,
            "DifOrdAna",
            "DifOrdAna > nX+nY-2 (the maximum diffraction order possible).",
        )
    checks.error_if(
        rlt(PI / p.dt_analysis, ev_to_hartree(p.e_max_wave_packet)), "dtAna", "dtAna is too big for fine analysis."
    )


def check_wave_packet_and_state(p: InputParameters, z_asymptotic: float, checks: FindingList) -> None:
    # --- initial wave packet ---
    checks.error_if(rle(p.e_min_wave_packet, 0.0), "EMinWP", "EMinWP <= 0.")
    checks.error_if(rle(p.e_max_wave_packet, 0.0), "EMaxWP", "EMaxWP <= 0.")
    checks.error_if(rge(p.e_min_wave_packet, p.e_max_wave_packet), "EMinWP", "EMinWP >= EMaxWP.")
    checks.error_if(rlt(p.z0, z_asymptotic), "Z0", "Z0 < ZAsymptotic.")
    checks.error_if(rgt(p.z0, p.zsp_opt_start), "Z0", "Z0 > ZSpOptStart.")

    # --- initial quantum state ---
    checks.error_if(p.v_in < 0, "VIn", "VIn < 0.")
    checks.error_if(p.j_in < 0, "JIn", "JIn < 0.")
    checks.error_if(abs(p.mj_in) > p.j_in, "MJIn", "|MJIn| > JIn.")
    checks.error_if(p.v_in > p.nr - 1, "VIn", "VIn > nR-1.")
    checks.error_if(p.j_in > p.j_max, "JIn", "JIn > JMax.")  # unreachable in the program since r1726, as there
    checks.error_if(abs(p.mj_in) > p.mj_max, "MJIn", "MJIn > MJMax.")
    # Upstream r1726: VAna counts levels (v = 0 to VAna-1) and JAna is checked the same way.
    if p.v_analysis is not None:
        checks.error_if(p.v_in >= p.v_analysis, "VIn", "VIn >= VAna.")
    if p.j_analysis is not None:
        checks.error_if(p.j_in >= p.j_analysis, "JIn", "JIn >= JAna.")
    # analysis D3 (C++ only): the analysis takes at most nR-1 vibrational levels (docs/fortran_deviations.md).
    checks.error_if(
        p.v_in == p.nr - 1,
        "VIn",
        f"VIn = nR-1 = {p.v_in}: the analysis takes at most nR-1 vibrational levels (v = 0 to nR-2), so the "
        "initial level would not be analysed; use VIn < nR-1.",
    )


def check_potential_and_absorbers(p: InputParameters, z_asymptotic: float, checks: FindingList) -> None:
    # --- potential ---
    if p.impose_v_min and p.impose_v_max:
        checks.error_if(rgt(p.v_min, p.v_max), "VMin", "VMin > VMax.")
    checks.error_if(p.symmetry_tolerance <= 0.0, "SymmetryTolerance", "SymmetryTolerance <= 0.")

    # --- optical potential in Z, scattering grid ---
    z_end = p.z_start + (p.nz - 1) * p.dz
    checks.error_if(rgt(p.z_opt_start, p.z_opt_end), "ZOptStart", "ZOptStart > ZOptEnd.")
    checks.error_if(rgt(p.z_opt_end, z_end), "ZOptEnd", "ZOptEnd is too big for the scattering grid in Z.")
    checks.error_if(rlt(p.z_opt_start, p.z_analysis), "ZOptStart", "ZOptStart < ZAna.")
    checks.error_if(rgt(p.z_analysis, p.z_opt_start), "Zana", "ZAna > ZOptStart.")

    # --- optical potential in r ---
    checks.error_if(rgt(p.r_opt_start, p.r_opt_end), "ROptStart", "ROptStart > ROptEnd.")
    checks.error_if(rgt(p.r_opt_end, p.r_start + (p.nr - 1) * p.dr), "ROptEnd", "ROptEnd is too big for the grid in r.")

    # --- optical potential in Z, specular grid (ProjOp = 2) ---
    if p.projection_operator == 2:
        checks.error_if(rgt(p.zsp_opt_start, p.zsp_opt_end), "ZSpOptStart", "ZSpOptStart > ZSpOptEnd.")
        checks.error_if(
            rgt(p.zsp_opt_end, p.z_start + (p.nz_specular - 1) * p.dz),
            "ZSpOptEnd",
            "ZSpOptEnd is too big for the specular grid in Z.",
        )
        checks.error_if(rlt(p.zsp_opt_start, z_asymptotic), "ZSpOptStart", "ZSpOptStart < ZAsymptotic.")

    # --- parallel motion ---
    checks.error_if(rlt(p.e_parallel, 0.0), "EParallel", "EParallel < 0.")
    checks.error_if(rlt(p.phi_parallel, 0.0), "PhiParallel", "PhiParallel < 0.")
    checks.error_if(rgt(p.phi_parallel, 360.0), "PhiParallel", "PhiParallel > 360 degrees.")

    # --- switching function ---
    if p.use_switching:
        checks.error_if(rge(p.z_min_switch, p.z_max_switch), "ZminSwitch", "ZminSwitch >= ZmaxSwitch.")


def check_pes_cut(p: InputParameters, system: SystemParameters, checks: FindingList) -> None:
    # --- 2D cut of the PES (write2DPES) ---
    # alat and blat are unknown here for a CRP data set (with CRPDistrtSurface: the slab snapshot's lattice;
    # potential_systems.py): then the upper bounds of x and y are not checked.
    if p.write_2d_pes:
        checks.error_if(rlt(p.x_2d_pes, 0.0), "x2DPES", "x2DPES < 0.")
        checks.error_if(system.lattice_a is not None and rgt(p.x_2d_pes, system.lattice_a), "x2DPES", "x2DPES > alat.")
        checks.error_if(rlt(p.y_2d_pes, 0.0), "y2DPES", "y2DPES < 0.")
        checks.error_if(system.lattice_b is not None and rgt(p.y_2d_pes, system.lattice_b), "y2DPES", "y2DPES > blat.")
        checks.error_if(rlt(p.theta_2d_pes, 0.0), "theta2DPES", "theta2DPES < 0.")
        checks.error_if(rgt(p.theta_2d_pes, 360.0), "theta2DPES", "theta2DPES > 360.")
        checks.error_if(rlt(p.phi_2d_pes, 0.0), "phi2DPES", "phi2DPES < 0.")
        checks.error_if(rgt(p.phi_2d_pes, 360.0), "phi2DPES", "phi2DPES > 360.")

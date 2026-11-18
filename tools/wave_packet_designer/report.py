"""Text report of an Evaluation, for the command line.

Shadows: src/report/format.cpp (format_bytes, Table::print), src/report/sections.cpp (print_memory_estimate)

Presentation layer (README, "Architecture"): formats what evaluate() computed and computes nothing itself
(apart from unit conversions for display). One function per section; each returns its lines.
"""

import textwrap

from evaluation import SETTINGS_STAGE, SYMMETRY_GIVEN, SYMMETRY_NOT_GIVEN, Evaluation, Suggestion
from program.findings import Finding, Severity
from program.json_input import key_of_field
from program.memory_estimate import MemoryEstimate
from program.unit_conversion import hartree_to_ev

LABEL_WIDTH = 46  # label plus dotted leader, as in the program's set-up output
TEXT_WIDTH = 116  # long messages are wrapped to this width


def format_check_report(e: Evaluation, source: str, reading_findings: list[Finding]) -> str:
    """The full report of the `check` command."""
    lines = input_section(e, source)
    lines += findings_section(reading_findings + list(e.findings), e.stopped_at)
    if e.derived is not None and e.grids is not None:
        lines += wave_packet_section(e)
        lines += z_grids_section(e)
        lines += absorbers_section(e)
        lines += r_grid_section(e)
        lines += gas_phase_section(e)
        lines += design_section(e)
        lines += absorber_response_section(e)
    return "\n".join(lines) + "\n"


def format_suggestion_report(s: Suggestion, original: Evaluation, source: str) -> str:
    """The report of the `suggest` command: the changed keys, then the checks of the changed input."""
    lines = section(f"SUGGESTED Z GRIDS AND ABSORBERS for {source} (inputs/README.md, section 2)")
    p, new = original.parameters, s.evaluation.parameters
    lines.append(f"  {'key':14} {'input':>14} {'suggested':>14}")
    for field in s.changes:
        old, value = getattr(p, field), getattr(new, field)
        mark = "" if old == value else "  *"
        lines.append(f"  {key_of_field(field):14} {old!s:>14} {value!s:>14}{mark}")
    lines.append("  (* changed)")
    for key, factor in (("ZOptFactor", p.z_opt_factor), ("ZSpOptFactor", p.zsp_opt_factor)):
        if factor is not None and key.replace("Factor", "Energy") not in [key_of_field(f) for f in s.changes]:
            lines.append(f"  {key} {factor:g} kept: the absorber is tuned by its factor, not by an energy")
    lines.append("  Not checked (section 4): the local momenta in wells or the exit channel, and the time step.")
    lines += findings_section(list(s.evaluation.findings), s.evaluation.stopped_at)
    lines += design_section(s.evaluation)
    lines += absorber_response_section(s.evaluation)
    return "\n".join(lines) + "\n"


def format_memory_report(e: Evaluation, source: str, reading_findings: list[Finding]) -> str:
    """The report of the `memory` command: the sizes that set the memory, and the CPU and GPU estimates."""
    lines = input_section(e, source)
    lines += findings_section(reading_findings + list(e.findings), e.stopped_at)
    if e.resources is not None:
        lines += sizes_section(e)
        if e.resources.memory is None:
            lines.append("  (no estimate: an analysis-only run needs the fine analysis's memory, not copied)")
        else:
            lines += memory_section("MEMORY ESTIMATE", e.resources.memory)
        if e.resources.gpu_memory is not None:
            lines += memory_section("GPU MEMORY ESTIMATE (with --device cuda)", e.resources.gpu_memory)
    return "\n".join(lines) + "\n"


def format_reading_failure(source: str, reading_findings: list[Finding]) -> str:
    """The report when the program would stop while reading the file."""
    lines = section(f"INPUT  {source}")
    lines += findings_section(reading_findings, "parse_input")
    return "\n".join(lines) + "\n"


def section(title: str) -> list[str]:
    return ["", title, "-" * len(title)]


def row(label: str, value: str, unit: str = "") -> str:
    """'  label ........ value unit', aligned like the program's key-value lists."""
    dots = "." * max(2, LABEL_WIDTH - len(label) - 1)
    return f"  {label} {dots} {value}{' ' + unit if unit else ''}"


def input_section(e: Evaluation, source: str) -> list[str]:
    p = e.parameters
    lines = section(f"INPUT  {source}")
    lines.append(row("SystemID", p.system_id))
    if e.system is not None and e.z_asymptotic is not None:
        if p.use_switching:
            origin = "ZmaxSwitch, as UseSwitching is on"
        elif e.settings.z_asymptotic is not None:
            origin = "given"
        else:
            origin = "built in"
        lines.append(row("Z_asymptotic", f"{e.z_asymptotic:g} ({origin})", "bohr"))
        lines.append(row("Point group", e.system.point_group))
    return lines


def findings_section(findings: list[Finding], stopped_at: str | None) -> list[str]:
    lines = section("PROGRAM CHECKS (what spodvr would do with this input)")
    if stopped_at is None:
        lines.append("  spodvr accepts the input: every set-up stage the designer copies passes.")
    elif stopped_at == SETTINGS_STAGE:
        lines.append("  The designer needs more settings to evaluate this input (this is not a program error):")
    else:
        lines.append(f"  spodvr stops in {stopped_at} (it reports the first error; all of this stage are listed):")
    for f in findings:
        mark = "ERROR  " if f.severity is Severity.ERROR else "WARNING"
        key = f"{f.key}: " if f.key else ""
        text = f"{mark} [{f.origin}] {key}{f.message}"
        lines.append(textwrap.fill(text, TEXT_WIDTH, initial_indent="  ", subsequent_indent="            "))
    return lines


def wave_packet_section(e: Evaluation) -> list[str]:
    p, d, s = e.parameters, e.derived, e.packet_summary
    lines = section("INITIAL WAVE PACKET (Z)")
    lines.append(row("Energy range (EMinWP .. EMaxWP)", f"{p.e_min_wave_packet:g} .. {p.e_max_wave_packet:g}", "eV"))
    lines.append(row("Mean momentum kz (towards the surface)", f"{d.kz:.4f}", "1/bohr"))
    lines.append(row("Width in k, sigma_k = sqrt(alpha)", f"{s.sigma_k:.4f}", "1/bohr"))
    lines.append(row("Width in Z, sigma_Z = 1/(2 sigma_k)", f"{s.sigma_z:.4f}", "bohr"))
    lines.append(row("Centre Z0", f"{p.z0:g}", "bohr"))
    if e.packet is None:
        lines.append("  (no wave packet: none in this run mode, or the program stops before building it)")
        return lines
    info, z = e.packet.info, e.grids.z_specular
    lines.append(row("Allowed region (specular indices)", f"{info.i_low} .. {info.i_high}"))
    lines.append(row("Allowed region (Z)", f"{z.point(info.i_low):.4f} .. {z.point(info.i_high):.4f}", "bohr"))
    lines.append(row("Norm on the specular grid", f"{info.norm:.16f}"))
    lines.append(row("Norm outside the region, norm - norm inside", f"{s.norm_outside:.3e}"))
    lines.append(row("1 - norm inside (must be <= TolSp)", f"{s.tolsp_measure:.3e} (TolSp {p.tol_specular:g})"))
    lines.append(row("|psi(k)| at the grid edge / maximum", f"{s.edge_to_maximum:.3e}"))
    return lines


def z_grids_section(e: Evaluation) -> list[str]:
    g = e.grids
    z, zsp = g.z_scattering, g.z_specular
    lines = section("Z GRIDS")
    lines.append(
        row("Scattering grid (nZ)", f"{z.size} points, {z.start:.4f} .. {z.last:.4f}, step {z.spacing:g}", "bohr")
    )
    lines.append(row("Specular grid (nZSpecular)", f"{zsp.size} points, {zsp.start:.4f} .. {zsp.last:.4f}", "bohr"))
    lines.append(row("Largest momentum k_grid = pi/dZ", f"{e.packet_summary.grid_momentum_limit:.4f}", "1/bohr"))
    lines.append(row("Kinetic energy at k_grid", f"{hartree_to_ev(g.kz_scattering.max_kinetic_energy):.4f}", "eV"))
    lines.append(row("Interaction region nzpot, Z < Z_asymptotic", f"{g.nz_potential} points"))
    lines.append(row("Analysis line Zana (index)", f"{g.z_analysis:.4f} ({g.iz_analysis})", "bohr"))
    return lines


def absorbers_section(e: Evaluation) -> list[str]:
    lines = section("ABSORBERS")
    lines.append(
        f"  {'name':10} {'grid':13} {'start .. end (bohr)':>20} {'indices':>11} {'target (eV)':>12} "
        f"{'L/lambda':>9} {'strength (Ha)':>14}"
    )
    for a in e.absorbers:
        target = f"{a.target_energy:g}" if a.target_energy is not None else "factor"
        reduced = f"{a.reduced_length:.3f}" if a.reduced_length is not None else "-"
        lines.append(
            f"  {a.name:10} {a.grid:13} {a.start:9.4f} .. {a.end:7.4f} "
            f"{a.indices.first:5} .. {a.indices.last:<4} {target:>8} {reduced:>9} {a.strength:14.6e}"
        )
    return lines


def r_grid_section(e: Evaluation) -> list[str]:
    r = e.grids.r
    lines = section("r GRID")
    lines.append(row("r grid (nR)", f"{r.size} points, {r.start:.4f} .. {r.last:.4f}, step {r.spacing:g}", "bohr"))
    lines.append(
        row("Kinetic energy at pi/dR (reduced mass)", f"{hartree_to_ev(e.grids.kr.max_kinetic_energy):.4f}", "eV")
    )
    return lines


def gas_phase_section(e: Evaluation) -> list[str]:
    lines = section("GAS PHASE (r)")
    gas = e.gas_phase
    if gas is None:
        lines.append("  (not computed: VIn lies beyond the levels of the gas-phase curve)")
        return lines
    lines.append(f"  States of {gas.curve.source},")
    lines.append(f"  for the input's masses, on the curve's r grid (r = {gas.r[0]:g} .. {gas.r[-1]:g} bohr):")
    if e.gas_phase_indicative:
        lines.append(
            "  (the program solves its own potential's curve on the input's r grid: these are only indicative)"
        )
    lines.append(row("Open vibrational levels", f"{gas.v_max_open}"))
    lines.append(row("Highest open j", f"{gas.j_max_open}"))
    lines.append(row("Vibrational levels analysed (VAna)", f"{gas.n_v_analysis}"))
    lines.append(row("Highest j analysed (JAna)", f"{gas.j_max_analysis}"))
    for warning in gas.warnings:
        lines.append(f"  WARNING {warning}")
    return lines


def sizes_section(e: Evaluation) -> list[str]:
    p, r = e.parameters, e.resources
    lines = section("SIZES")
    lines.append(row("Theta points", f"{r.angular.n_theta}"))
    lines.append(row("Phi points", f"{r.angular.n_phi}"))
    lines.append(row("J range (step)", f"{r.angular.j_min} .. {r.angular.j_max} ({r.angular.delta_j})"))
    lines.append(row("MJMax (basis)", f"{r.angular.mj_max}"))
    lines.append(row("Rotational basis functions", f"{r.rotational_basis_size}"))
    lines.append(row("Grid points (x, y, Z, r)", f"{p.nx} x {p.ny} x {p.nz} x {p.nr}"))
    lines.append(row("Interaction region nzpot", f"{e.grids.nz_potential} Z points"))
    n_xy = p.nx * p.ny
    if r.symmetry_source == SYMMETRY_GIVEN:
        symmetry = f"{r.symmetry.n_representatives} of {n_xy} (x, y) points, {r.symmetry.n_operations} operations"
    elif r.symmetry_source == SYMMETRY_NOT_GIVEN:
        symmetry = f"not given: all {n_xy} (x, y) points (an upper bound; see --representatives)"
    else:
        symmetry = "none (UseSymmetry false)"
    lines.append(row("Surface symmetry of the potential", symmetry))
    lines.append(row("OpenMP threads (work arrays)", f"{r.threads}"))
    return lines


def memory_section(title: str, estimate: MemoryEstimate) -> list[str]:
    """The program's table of the large arrays and the totals (print_memory_estimate)."""
    lines = section(title)
    header = ("array", "values", "size", "note")
    rows = [(item.name, str(item.values), format_bytes(item.bytes), item.note) for item in estimate.items]
    widths = [max(len(cells[column]) for cells in [header, *rows]) for column in range(4)]
    for cells in [header, *rows]:
        lines.append(("  " + "  ".join(cell.rjust(width) for cell, width in zip(cells, widths, strict=True))).rstrip())
        if cells is header:
            lines.append("-" * sum(width + 2 for width in widths))
    lines.append("")
    if estimate.propagation_bytes == 0:  # potential-only run
        lines.append(row("Total", format_bytes(estimate.setup_bytes)))
    else:
        lines.append(row("Setting up the propagation", format_bytes(estimate.setup_bytes)))
        lines.append(row("During the propagation", format_bytes(estimate.propagation_bytes)))
    lines.append("  Small tables and FFT plans are not included.")
    return lines


def design_section(e: Evaluation) -> list[str]:
    lines = section("DESIGN (the recommendations of inputs/README.md)")
    for c in e.design:
        mark = {True: "ok     ", False: "ADVICE ", None: "info   "}[c.passed]
        source = "the designer's limit" if c.designer_limit else "recommended"
        limit = "" if c.relation is None else f" ({source} {c.relation} {c.limit:g})"
        unit = f" {c.unit}" if c.unit else ""
        lines.append(f"  {mark} {c.name}: {c.value:.4g}{unit}{limit}  [section {c.section}]")
        if c.passed is not True:
            lines.append(f"            {c.advice}")
    return lines


def absorber_response_section(e: Evaluation) -> list[str]:
    lines = section("ABSORBERS: REFLECTION + TRANSMISSION (1D, inputs/README.md section 7)")
    lines.append("  Energies: specular absorbers EMinWP .. EMaxWP; the others 0.05 eV .. the available energy")
    lines.append("  (for ROpt the designer's choice: the README sets no range).")
    for r in e.absorber_responses:
        lines.append(f"  {r.name:10} E (eV)  " + " ".join(f"{x:8.3f}" for x in r.energies))
        lines.append(f"  {'':10} R + T   " + " ".join(f"{x:8.1e}" for x in r.total))
    return lines


def format_bytes(n_bytes: int) -> str:
    """A size in B, KiB, MiB, GiB or TiB with one decimal (format_bytes)."""
    if n_bytes < 1024:
        return f"{n_bytes} B"
    units = ("KiB", "MiB", "GiB", "TiB")
    size = n_bytes / 1024.0
    unit = 0
    while size >= 1024.0 and unit < 3:
        size = size / 1024.0
        unit = unit + 1
    return f"{size:.1f} {units[unit]}"

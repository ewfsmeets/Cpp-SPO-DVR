"""Text report of an Evaluation, for the command line.

Presentation layer (README, "Architecture"): formats what evaluate() computed and computes nothing itself
(apart from unit conversions for display). One function per section; each returns its lines.
"""

import textwrap

from evaluation import SETTINGS_STAGE, Evaluation
from program.findings import Finding, Severity
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

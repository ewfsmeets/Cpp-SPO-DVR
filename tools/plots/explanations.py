"""What each figure shows, how to read it and what good looks like (layer 3), keyed by its PLOTS name: a short
title (the report's heading) and the explanation the report prints under the figure. tests/test_runs_figures.py
checks that every PLOTS entry has both and that every text has an entry."""

TITLES = {
    "verdict": "Verdict",
    "results": "Results against energy",
    "averages": "Norms and mean positions in time",
    "gas_phase": "Gas phase",
    "bk": "Balint-Kurti coefficients and convergence in time",
}

FIGURES = {
    "verdict": (
        "The findings about the run's completeness and convergence, worst first, each with its value, limit and "
        "the section of inputs/README.md that explains it. Good: no errors; warnings that you understand."
    ),
    "results": (
        "The scattering and reaction probabilities against the collision energy. P_reaction is 1 - P_scattered, so "
        "flux not yet analysed counts as reaction; a negative P_reaction shows the unitarity error. When more than "
        "one vibrational level is scattered into, P(v') in a panel of its own. The mean final energies of the "
        "scattered molecules (the diffraction energy is the final parallel kinetic energy), left out where "
        "P_scattered < 1e-6."
    ),
    "averages": (
        "Left: the norm, the absorbed part 1 - norm and the rest norm (the specular wave between the two analysis "
        "lines) against time. Right: <Z> and <r>. Good: the norm drops mainly once the packet reaches the "
        "absorbers (an early slow drop is the initial state's tail in the r absorber: see the gas phase), and the "
        "rest norm goes to zero by tEnd. Dotted lines: where restart segments start."
    ),
    "gas_phase": (
        "Left: the gas-phase potential V_asy(r) with the levels of the initial j (dashed: not bound on the r grid) "
        "and the r absorber. Right: |psi(r)| of the lowest levels on a log axis. Good: the analysed levels are bound "
        "and |psi| has fallen by many decades at the r absorber; what is left there is absorbed from the start, "
        "also from the initial state (counted as reaction)."
    ),
    "bk": (
        "Left: |c(t)|^2 at the analysis line of the strongest channels: when each channel's flux passes (the "
        "initial channel is analysed at the specular line, which the incoming packet's slow tail also crosses, "
        "inward, before the scattered flux arrives). Right: "
        "P_scattered(E; t) / P_scattered(E; tEnd), the scattering probability if the time sums stopped at t "
        "(white = 1, the final value; purple above, brown below). Good: every |c|^2 has fallen to its floor and "
        "every energy reaches 1 well before tEnd."
    ),
}

"""What each figure shows, how to read it and what good looks like (layer 3), keyed by its PLOTS name: a short
title (the report's heading) and the explanation the report prints under the figure. tests/test_runs_figures.py
checks that every PLOTS entry has both and that every text has an entry."""

from compute.potential import CONTOUR_LINE_EVERY, CONTOUR_STEP_EV, ENERGY_CAP_EV, MAX_ELBOW_SITES

TITLES = {
    "verdict": "Verdict",
    "results": "Results against energy",
    "state_to_state": "Final rovibrational states",
    "diffraction": "Diffraction",
    "scattering_angles": "Scattering angles",
    "time_delay": "Time delays",
    "fine_analysis": "On-the-fly and fine analysis",
    "averages": "Norms and mean positions in time",
    "gas_phase": "Gas phase",
    "densities": "The packet in Z and r against time",
    "density_movie": "The packet over (Z, r)",
    "potential_maps": "The potential over the unit cell",
    "elbows": "Elbow plots and the minimax barrier",
    "orientation": "The potential against the orientation",
    "bk": "Balint-Kurti coefficients and convergence in time",
}

FIGURES = {
    "verdict": (
        "The findings about the run's completeness and convergence, worst first, each with its value, limit and "
        "the section of inputs/README.md that explains it. Good: no errors; warnings that you understand."
    ),
    "results": (
        "Top: the scattering probability, and the reaction probability P_reaction = 1 - P_scattered on an axis of "
        "its own. Flux not yet analysed counts as reaction; the grey band is the unitarity error (-min P_reaction): "
        "a P_reaction inside it cannot be told from 0. Bottom: the mean final energies of the scattered molecules "
        "(left out where P_scattered < 1e-6; the diffraction energy is the final parallel kinetic energy, the line "
        "the initial one), and P(v') when more than one vibrational level is scattered into."
    ),
    "state_to_state": (
        "Top left: P(v', j') of the most probable final levels against the collision energy. Top right: P(v', j') "
        "at one energy (the middle one, or --energy); the highest j' of the basis is marked JMax. Bottom left: "
        "P(mJ') within each j', with the surface normal as quantisation axis. Bottom right: the alignment A0(2): "
        "-1 for a molecule rotating in a plane that contains the normal (cartwheel), up to (2j' - 1)/(j' + 1) "
        "(dotted; 1 for j' = 2) for one rotating about the normal (helicopter), 0 for no preference. Left out: "
        "energies where P(j') is rounding noise, and a j' above the run's MJMax (its mJ' basis is truncated)."
    ),
    "diffraction": (
        "Left: the analysed diffraction channels at their parallel momentum K (open circles), the probability at "
        "one energy as filled circles of area proportional to P (the strongest labelled (n,m)), summed over v', j' "
        "and mJ'; the arrow is the incident parallel momentum, the dotted line the plane of incidence (at normal "
        "incidence a chosen azimuth, --azimuth). Right: P of each diffraction order, and the fraction of "
        "P_scattered out of the plane, against energy. At an azimuth off the mirror lines only the specular "
        "channel (0,0) lies in the plane. Good: the outermost analysed channels carry no flux (the verdict checks)."
    ),
    "scattering_angles": (
        "Every open final state (v', j', n, m) at one energy, summed over mJ', at its azimuth phi_f (from the first "
        "lattice vector, the x axis of the spot map) and its polar angle theta_f from the surface normal (the "
        "radius), marker area proportional to P. The open circle is the specular direction: the angle of "
        "incidence at this energy (the parallel momentum is fixed, so it changes with the energy)."
    ),
    "time_delay": (
        "Left: the time delay d arg S / dE of the most probable channels, each S referred to the line where the "
        "channel is analysed (the initial state: the specular line), relative to a reflection at Z = 0 (negative: "
        "the molecule turned before reaching Z = 0); shown where the channel carries at least 1e-3 of P_scattered. "
        "A spike where its probability dips (right) is a near-zero of S, not a time. The initial state's delay is "
        "only as good as its specular line: fixed close to Z0, it also counts the incoming packet's slow tail. "
        "Right: their probabilities."
    ),
    "fine_analysis": (
        "Left: P_scattered of the on-the-fly analysis, of the fine analysis on its spectral grid (dots) and "
        "interpolated to the analysis energies (squares). Right: their difference. Good: the two agree to well "
        "below the effect you study; a difference that changes sign from point to point is the interpolation "
        "between the fine grid's points. The analysis energies outside the fine grid are listed below."
    ),
    "bk": (
        "Left: |c(t)|^2 at the analysis line of the strongest channels: when each channel's flux passes. Right: "
        "how far P_scattered(E) would be from its final value if the time sums stopped at t, "
        "|P(E; t) / P(E; tEnd) - 1|, on a log colour scale. Good: every |c|^2 has fallen to its floor and every "
        "energy is dark (below 1e-3) well before tEnd."
    ),
    "averages": (
        "Left: the norm, 1 - norm (the norm removed since t = 0: by the absorbers and, where the potential couples "
        "rotational states, by the angular transform, which keeps only the rotational basis up to JMax and MJMax) "
        "and the rest norm (the specular wave between the two analysis lines) against time. Right: <Z> and <r>, "
        "averaged over the norm left at that time. Good: the norm drops mainly once the packet reaches the "
        "absorbers (an early slow drop is the initial state's tail in the r absorber: see the gas phase), and the "
        "rest norm goes to zero by tEnd. A thin vertical line marks where a restart segment starts."
    ),
    "gas_phase": (
        "Left: the gas-phase potential V_asy(r) with the levels of the initial j (dashed: not bound on the r grid) "
        "and the r absorber. Right: |psi(r)| of the lowest levels on a log axis. Good: the analysed levels are bound "
        "and |psi| has fallen by many decades at both ends of the r grid; what is left at the outer end is absorbed "
        "from the start, also from the initial state (counted as reaction). The verdict checks both ends."
    ),
    "densities": (
        "Top: the probability per bohr of Z (summed over r) and of r (summed over Z) against time, from densities/ "
        "(writeDensities: the probability per (Z, r) cell, the sums behind averages.dat; the scattering and the "
        "specular wave summed): the packet comes in, turns and leaves. The marks name the wave each absorber acts "
        "on: ZOpt only the scattering wave, ZSpOpt only the specular wave (the initial channel), ROpt both. On the Z "
        "map the specular analysis line starts before ZSpOpt and moves inward behind the incoming packet until it "
        "is fixed; the specular wave is analysed there only once it is fixed. A faint band at the outer end of r "
        "from t = 0 is the initial state's tail that wraps around the periodic r grid into the r absorber, not "
        "dissociation. Bottom: the norm per region (r beyond r_c, where the r absorber starts: dissociating, counted "
        "as reaction; for smaller r the Z bands below zasymp, up to Zana, up to the scattering grid's end and "
        "beyond it), 1 - norm (removed: by the absorbers and by the rotational basis cut at JMax, MJMax, see the "
        "averages) and, in the absorbers' colours, the probability inside each absorber. Good: the norm leaves "
        "through the absorbers, and little reaches r >= r_c before the packet reaches the surface."
    ),
    "density_movie": (
        "rho(Z, r) of the scattering and the specular wave summed (bohr^-2), on one log colour scale for all "
        "frames, with zasymp, Zana, the scattering grid's end, where each absorber begins to absorb (with the wave "
        "it acts on) and the specular analysis line of each frame marked. The PNG beside the GIF is the last frame."
    ),
    "potential_maps": (
        "Over the unit cell (Cartesian coordinates, the first lattice vector along X; each (x, y) grid point's cell "
        "drawn around it, periodically): left, the minimax barrier on V_min(Z, r) of each (x, y) point (see the "
        "elbows; hatched where it has no saddle); middle, V_min at the (Z, r) of the lowest of these barriers (or at "
        "--z, --r); right, V_min at the "
        "same Z and the gas-phase bond length. V_min is the lower envelope of V over the run's DVR nodes at fixed "
        "(x, y): each point has its own orientation, contours kink where the node changes, and no rotational or "
        "zero-point energy is in it. V is the propagation's (switching and the VMin/VMax limits included), relative "
        "to min V_asy, the lowest V on the r grid of asymptotic_potential.dat (Z = zasymp, theta = 90 deg); "
        f"hatched: more than {ENERGY_CAP_EV:g} eV above it (not drawn). The sites of PlotSites and the lowest "
        "barrier's point are marked."
    ),
    "elbows": (
        f"Top: V_min(Z, r) at each site of PlotSites (the first {MAX_ELBOW_SITES}; or, without sites, at the (x, y) "
        "point with the lowest barrier), with the minimax barrier: the lowest possible highest point of a path from "
        "the gas-phase side (the lowest point just below zasymp) to the dissociated side (r >= r_c, where the r "
        "absorber starts), marked x, and one path that reaches it (dashed; below the barrier any route will do: it "
        "is not a minimum-energy path). It is the barrier at the run's nodes: the PES's own can be lower when the "
        "best orientation lies between nodes. When the path's highest point is its end at r_c (no saddle before the "
        "dissociated side) or its start (V falls from the gas-phase side: no barrier), the title says so and no x is "
        "drawn. Bottom: V at the orientation of the top barrier's point (its DVR node, the first of equal ones), with "
        "its own pass: a fixed orientation, usually higher than the top one, along which V_min may change "
        "orientation. theta is measured from the surface "
        "normal (90 deg: parallel), phi from the first lattice vector. A barrier equal to VMax is the clamp. Filled "
        f"contours every {CONTOUR_STEP_EV:g} eV, lines every {CONTOUR_LINE_EVERY * CONTOUR_STEP_EV:g} eV; hatched: "
        f"more than {ENERGY_CAP_EV:g} eV above the gas-phase minimum."
    ),
    "orientation": (
        "V at every DVR node (theta from the surface normal, phi from the first lattice vector) at the barrier "
        "point of the site with the lowest barrier: the orientation dependence the propagation sees there; the dots "
        "are the nodes, the only orientations the run has, each drawn as the cell around it. A homonuclear molecule "
        "has theta <= 90 deg only, and at theta = 90 deg phi and phi + 180 deg are the same orientation. Hatched: "
        f"more than {ENERGY_CAP_EV:g} eV above the gas-phase minimum."
    ),
}

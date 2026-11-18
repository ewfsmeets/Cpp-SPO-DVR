"""The explanations of the window: what each panel shows and how to read it, and what each slider is for.

Presentation data (README, "Architecture"): plain text, one entry per panel (keyed by its draw function in
plots.py) and per slider or button (keyed by its field in gui_mpl.py). The physics stays in inputs/README.md;
each text says where to read more. Edit the texts here; the window only shows them.
"""

# What a panel shows, how to read it, and what a good design looks like there.
PANELS = {
    "draw_z_layout": (
        "Z: the two Z grids as lanes at the top. The scattering grid (nZ points) holds the full 6D wave; the "
        "6D potential is stored only in its interaction region, Z < Z_asymptotic. The specular grid reaches "
        "further out and holds only the incoming molecule, so the packet starts there, inside the scattering "
        "grid's length but on the other grid. Each grid ends in an absorber (ZOpt, ZSpOpt). Below: the initial "
        "packet |psi(Z)|^2 on a log scale; Zana is the analysis line, where outgoing flux is analysed. Good: "
        "the packet lies inside its allowed region (arrow), Z0 >= Zana + 8 sigma_Z, and its outer tail "
        ">= 2 bohr before ZSpOpt. inputs/README.md sections 1, 2, 5."
    ),
    "draw_momentum": (
        "Z momentum: |psi(k)|^2 of the initial packet. It moves towards the surface, so its momenta lie around "
        "-k0, between the momenta of EMaxWP and EMinWP (ticks at the bottom). The grid holds only |k| <= "
        "k_grid = pi/dZ: whatever lies beyond -k_grid re-appears near +k_grid, moving outwards (shaded), and "
        "runs into the specular absorber. Good: the curve has fallen far below 1e-15 before -k_grid, and the "
        "outward part is many orders below TolSp. inputs/README.md sections 3 and 4."
    ),
    "draw_absorber_response": (
        "Absorbers: each absorber removes waves that reach the end of its grid, but never perfectly: a "
        "fraction R is reflected back into the grid and a fraction T passes through to the grid's end and "
        "re-enters at the other side. Both end up as errors in the probabilities. Solid: R; dotted: T; per "
        "absorber over the energies it must handle (specular: the packet's; scattering: 0.05 eV up to the "
        "available energy; r: up to the r kinetic energy the r grid holds). Lower is better; dashed: the "
        "levels of the README's H2 examples. R rising at low energy: the absorber is too short for slow waves "
        "(longer, or a lower target); T rising at high energy: too weak for fast ones (a higher target). "
        "inputs/README.md section 7."
    ),
    "draw_r_layout": (
        "r (the H-H distance): a bound molecule vibrates around r = 1.4 bohr; its wave function falls off "
        "further out. The r grid ends in an absorber (orange) that removes dissociating molecules, and it must "
        "not reach into the bound ones. Each curve: for one analysed vibrational level v, the norm of its most "
        "extended analysed state beyond r, i.e. what an absorber starting at r would remove. Black: the states "
        "as the program computes them on the input's r grid; grey: the same states converged on a finer grid. "
        "Good: the black curves are below the 1e-6 limit at ROptStart. Where black stays flat above grey, the "
        "r grid is too coarse (a finer dR), not ROptStart too small. For a potential whose gas phase is not "
        "copied the DZ input's curve stands in (indicative, no grey). inputs/README.md section 7."
    ),
    "draw_memory": (
        "Memory: the peak memory of the run on the CPU and with --device cuda on the GPU (bars), against what "
        "this machine has (orange marks: the RAM and the free GPU memory; spodvr refuses a GPU run that does "
        "not fit). It grows with nZ x nR x nX x nY x the rotational basis, the potential with the interaction "
        "region. With UseSymmetry the potential is stored at the representatives only: without "
        "--representatives the bars are an upper bound. inputs/README.md section 9; `memory` lists every array."
    ),
    "draw_time": (
        "Time: the classical paths Z(t) of the packet's fast, central and slow parts (the momenta of EMaxWP, "
        "k0 and EMinWP): in from Z0, reflected at the turning point Z_turn, and, if the reflection is elastic, "
        "out at the same speed. The slow one is back across Zana at 'elastic': a lower bound for tEnd. "
        "Rotational excitation and diffraction take normal energy, so many outgoing waves are slower; the "
        "orange path leaves with the lowest normal energy E_n,min you want converged and is back at 'E_n'. "
        "Good: tEnd (thick line) beyond both. Not included: trapping in wells, the spreading of the packet. "
        "inputs/README.md section 8."
    ),
    "draw_verdict": (
        "Verdict: first what spodvr would do with this input (its errors and set-up warnings), then the "
        "recommendations of inputs/README.md that the input misses, worst first: the run fails, the results "
        "suffer, it costs more than needed. Each with its value, the limit (marked where it is the designer's "
        "own choice) and what to change. Then the values the sliders led to, and what nothing here checks."
    ),
}

# What you are trying to achieve with a slider or button, which way to move it, and what it costs.
CONTROLS = {
    "e_min_wave_packet": (
        "EMinWP: the lowest normal energy you want results for (with a margin: the ends of the range are "
        "less accurate). Lower: a wider packet in k, slower, so a longer tEnd and a longer specular absorber "
        "(L/lambda at EMinWP). Keep EMaxWP <= 4 EMinWP. inputs/README.md section 3."
    ),
    "e_max_wave_packet": (
        "EMaxWP: the highest normal energy. Higher: faster momenta, which dZ must still hold (momentum "
        "panel), and more open channels (r panel). Keep EMaxWP <= 4 EMinWP, else the packet itself moves "
        "partly outwards. inputs/README.md section 3."
    ),
    "dz": (
        "dZ: the step of both Z grids. It must hold the packet's momenta, (pi/dZ - k0)/sigma_k >= 6 (momentum "
        "panel). Smaller: more Z points, so memory and time grow as 1/dZ; larger: wrap-around. Zana, nZ, "
        "nZSpecular and the absorber ends follow the new grid. inputs/README.md section 4."
    ),
    "z0": (
        "Z0: the centre of the initial packet. At least Zana + 8 sigma_Z, so no norm starts below the "
        "analysis line (spodvr's TolSp check), and its outer tail >= 2 bohr before the specular absorber. "
        "Further out only costs propagation time. inputs/README.md section 5."
    ),
    "zsp_opt_start": (
        "ZSpOptStart: where the specular absorber starts. Beyond the packet's outer tail plus ~2 bohr (Z panel); "
        "the specular grid grows to keep the absorber's length. inputs/README.md sections 5 and 7."
    ),
    "scattering absorber length": (
        "Scattering absorber length (from Zana to the end of nZ): >= ~4 bohr, so slow outgoing channels are "
        "absorbed too (R + T panel, low energies). Longer: a larger nZ and more memory. inputs/README.md "
        "section 7."
    ),
    "specular absorber length": (
        "Specular absorber length: as long as the strength table allows, L/lambda <= 15 at its target "
        "energy (the report warns beyond). It sees only the incoming channel; extra specular points are "
        "cheap. inputs/README.md section 7."
    ),
    "dr": (
        "dR: the step of the r grid. It must resolve the bound states up to the analysed v and j: on a too "
        "coarse grid their tails stay flat beyond the well (r panel: the program's curves above the "
        "converged grey ones), and the r absorber removes that part. Smaller: more r points (memory and time "
        "grow with nR). pi/dR must also hold the dissociating fragments' momenta. The r absorber keeps "
        "ending on the last r point. inputs/README.md sections 4 and 7."
    ),
    "nr": (
        "nR: the number of r points (a valid FFT size). With dR it sets the r range, RStart to RStart + "
        "(nR - 1) dR, which must hold the bound states and leave room for the r absorber. Memory and time grow "
        "with nR. inputs/README.md sections 7 and 9."
    ),
    "r_opt_start": (
        "ROptStart: where the r absorber starts. Beyond the extent of the analysed bound states (r panel), "
        "so bound molecules are not removed; the r grid must leave room to absorb dissociating fragments "
        "after it. inputs/README.md section 7."
    ),
    "z_opt_energy": (
        "ZOptEnergy: the target energy of the scattering absorber; it sets the strength. Lower absorbs slow "
        "waves but lets fast ones through, higher reflects slow ones; ~0.2 eV for H2. Watch its R + T "
        "curve. inputs/README.md section 7."
    ),
    "zsp_opt_energy": (
        "ZSpOptEnergy: the target energy of the specular absorber, about EMinWP: it absorbs the incoming "
        "channel reflected elastically at the packet's energies. Watch its R + T curve. inputs/README.md "
        "section 7."
    ),
    "r_opt_energy": (
        "ROptEnergy: the target energy of the r absorber, ~0.3 eV. It only keeps dissociating fragments from "
        "wrapping around in r. Watch its R + T curve. inputs/README.md section 7."
    ),
    "t_end": (
        "tEnd: how long the propagation runs. Outgoing flux that has not crossed Zana by then counts as "
        "reaction, so tEnd sets the lowest outgoing normal energy that is converged (the verdict says which; "
        "0.03 eV in the README's example). Longer: more cost, linear in tEnd. It moves in whole averaging "
        "intervals dtAve. inputs/README.md section 8."
    ),
    "lowest_outgoing_energy": (
        "E_n,min: the lowest outgoing normal energy tEnd must converge (a setting of the time budget, not an "
        "input key). Outgoing waves slower than this have not crossed Zana by tEnd and count as reaction. "
        "Lower: a longer tEnd (time panel). Near each channel threshold the outgoing waves are arbitrarily "
        "slow, so some cut is needed; 0.03 eV is the README's example. inputs/README.md section 8."
    ),
    "sliders": (
        "Sliders: each starts at the input's value, which its black tick keeps marking. Moving one changes the "
        "input in memory; the Z grids then follow as the recipe places them (Zana on the first grid point at "
        "or beyond Z_asymptotic, nZ and nZSpecular the smallest valid sizes that hold the absorbers). The "
        "verdict shows the values that result."
    ),
    "Suggest": "Suggest: applies the recipe of inputs/README.md section 2 to the Z grids, Z0, the Z absorbers and tEnd.",
    "Reset": "Reset: back to the input as read.",
    "Save": "Save: writes the input with every changed key to the file named in the box (never over one).",
    "Explain": "Explain: shows or hides these explanations on the panels.",
}

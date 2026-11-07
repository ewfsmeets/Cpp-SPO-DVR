#!/usr/bin/env python3
"""Item 11: why |psi(k)| of the initial packet at k = 0 and at the grid edge differ between the
codes by ~1.7e-15 (about 340 half units of the Fortran's 8 printed digits).

Both are tails of the momentum distribution, ~1e-10 of its peak, obtained by an FFT of the
Gaussian on the specular grid (C++ initial_wave_packet.cpp, Fortran CoherentStateFunc): a sum
of 256 terms of size up to ~1 that cancel to 1e-10. At k = 0 the DFT factor is 1 and at the
edge (k = pi/dZ, Z_j = j dZ) it is (-1)^j, so the correctly rounded value of the DFT of the
double-precision samples is available with math.fsum. The script computes it, an independent
FFT (numpy/pocketfft), and compares with the C++ (oneMKL) and Fortran (FFTW) values. Writes
data/items/item11_psi_k.json."""
from __future__ import annotations

import cmath
import math

import numpy as np

from common import RUN_DIRS, WORK, FortranBk, clean_log_lines, parse_fortran_real, printed_half_unit, save_item
from parse_logs import cpp_dotted, cpp_section, section


def samples(z, z0, alpha, kz):
    amp = (2.0 * alpha / math.pi) ** 0.25
    return [amp * cmath.exp(complex(-(alpha * ((zj - z0) * (zj - z0))), -(kz * zj))) for zj in z]


def exact_sum(g, sign_alternating):
    re = math.fsum((-1) ** j * x.real if sign_alternating else x.real for j, x in enumerate(g))
    im = math.fsum((-1) ** j * x.imag if sign_alternating else x.imag for j, x in enumerate(g))
    return complex(re, im)


def main():
    gtxt = (WORK / "potential" / "build-fast" / "grids.txt").read_text().split("\n")
    i = gtxt.index("z_specular 256")
    z = [float(x) for x in gtxt[i + 1:i + 257]]
    dz = 0.12
    scale = dz / math.sqrt(2.0 * math.pi)
    CL = clean_log_lines("cpu")
    cin = cpp_dotted(cpp_section(CL, "INPUT"))
    rbk = cpp_dotted(cpp_section(CL, "ROTATIONAL BASIS AND INITIAL WAVE PACKET"))
    alpha_c = float(cin["Wave packet alpha"].split()[0])
    kz_c = float(cin["Initial kz"].split()[0])
    fb = FortranBk(RUN_DIRS["fortran"] / "fort.21")
    wf = section(clean_log_lines("fortran"), "INFO FROM SETTING UP THE WAVE FUNCTION")
    f_text = {"zero": next(l.split()[-1] for l in wf if "at k = 0" in l),
              "edge": next(l.split()[-1] for l in wf if "end/beginning of the grid" in l)}
    c_val = {"zero": float(rbk["|psi(k)| at k = 0"]), "edge": float(rbk["|psi(k)| at the grid edge"])}
    out = {}
    for params, (alpha, kz) in (("cpp", (alpha_c, kz_c)), ("fortran", (fb.alpha, fb.kz))):
        g = samples(z, 18.0, alpha, kz)
        spec = np.fft.fft(np.array(g)) * scale
        out[params] = {"exact_zero": abs(exact_sum(g, False)) * scale, "exact_edge": abs(exact_sum(g, True)) * scale,
                       "numpy_zero": float(abs(spec[0])), "numpy_edge": float(abs(spec[128])),
                       "max": float(np.max(np.abs(spec))), "sum_abs_terms": math.fsum(abs(x) for x in g) * scale}
    ref = out["cpp"]
    res = {"inputs": {"alpha_cpp": alpha_c, "kz_cpp": kz_c, "alpha_fortran": fb.alpha, "kz_fortran": fb.kz},
           "by_parameters": out}
    for k in ("zero", "edge"):
        e = ref[f"exact_{k}"]
        f = parse_fortran_real(f_text[k])
        res[k] = {"exact_fsum": e, "cpp_onemkl": c_val[k], "fortran_fftw": f, "fortran_text": f_text[k],
                  "numpy_pocketfft": ref[f"numpy_{k}"],
                  "cpp_minus_exact": c_val[k] - e, "fortran_minus_exact": f - e, "numpy_minus_exact": ref[f"numpy_{k}"] - e,
                  "fortran_half_unit": printed_half_unit(f_text[k]),
                  "effect_of_fortran_alpha_kz": out["fortran"][f"exact_{k}"] - e,
                  "fortran_minus_exact_with_fortran_alpha_kz": f - out["fortran"][f"exact_{k}"],
                  "within_fortran_half_unit_with_fortran_alpha_kz": abs(f - out["fortran"][f"exact_{k}"]) <= printed_half_unit(f_text[k]),
                  "cpp_minus_fortran": c_val[k] - f}
    # ulps of the Fortran alpha and k_Z relative to the C++ values (quad-precision site)
    res["alpha_ulps"] = abs(fb.alpha - alpha_c) / np.spacing(alpha_c)
    res["kz_ulps"] = abs(fb.kz - kz_c) / np.spacing(kz_c)
    # rounding of the phase kz*Z near the packet centre (Z0 = 18): ulp(kz * Z0)
    res["phase_ulp_at_z0"] = float(np.spacing(kz_c * 18.0))
    # rounding scale of an FFT output element: eps * log2(N) * sum |terms| (standard bound)
    res["fft_bound"] = np.finfo(float).eps * math.log2(256) * ref["sum_abs_terms"]
    res["max_cpp"] = float(rbk["|psi(k)| maximum"])
    save_item("item11_psi_k", res)
    for k in ("zero", "edge"):
        r = res[k]
        print(k, "exact %.10e  C++-exact %.2e  F-exact %.2e  numpy-exact %.2e  (F alpha/kz effect %.1e, F half unit %.1e)" % (
            r["exact_fsum"], r["cpp_minus_exact"], r["fortran_minus_exact"], r["numpy_minus_exact"],
            r["effect_of_fortran_alpha_kz"], r["fortran_half_unit"]))
        print("   F - exact(with Fortran alpha, kz) = %.2e  within half unit: %s" % (
            r["fortran_minus_exact_with_fortran_alpha_kz"], r["within_fortran_half_unit_with_fortran_alpha_kz"]))
    print("alpha differs by %.0f ulp, kz by %.0f ulp; ulp of the phase kz*Z0: %.1e rad" % (res["alpha_ulps"], res["kz_ulps"], res["phase_ulp_at_z0"]))
    print("FFT rounding bound eps*log2(N)*sum|terms| = %.2e; sum|terms| = %.3f, max|psi(k)| = %.3f" % (
        res["fft_bound"], ref["sum_abs_terms"], ref["max"]))


if __name__ == "__main__":
    main()

"""
examples/mkdv.py

Modified KdV equation:  u_t + 6 u^2 u_x + u_xxx = 0.

Traveling-wave reduction (U(Lambda), Lambda = x - c t), integrated once
with zero integration constant:
    U'' = c U - 2 U^3   =>  a3 = -2, a2 = 0, a1 = c, a0 = 0.

Two wave-speed regimes, same equation:
    c > 0:  a3 = -2 < 0, f'(R0) = a1 = c > 0  -> APPLICABLE (hyperbolic branch)
    c < 0:  a3 = -2 < 0, f'(R0) = a1 = c < 0  -> NOT APPLICABLE (forbidden quadrant)

This pair demonstrates that applicability is governed purely by the
sign of the wave speed relative to a3, not by any special property of
mKdV itself -- the same equation succeeds or fails depending only on
which sign of c is chosen.

Run:  python3 examples/mkdv.py
"""

import sympy as sp
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gedf_applicability_checker import check_applicability, report

Lambda = sp.symbols("Lambda")

if __name__ == "__main__":
    print("Modified KdV: u_t + 6 u^2 u_x + u_xxx = 0\n")

    print("--- c = 1 (c > 0) -- APPLICABLE ---")
    result = check_applicability(a0=0, a1=1, a2=0, a3=-2)
    report(result)
    if result.applicable:
        H1, k = result.H1, result.k
        U = H1 / sp.cosh(k * Lambda)
        residual = sp.simplify(sp.diff(U, Lambda, 2) - (1 * U - 2 * U**3))
        print(f"\nCandidate solution: u(x,t) = {H1} * sech({k}(x - t))")
        print(f"Symbolic check U'' - (U - 2U^3) = {residual}  "
              f"({'verified' if residual == 0 else 'MISMATCH'})\n")

    print("--- c = -1 (c < 0) -- NOT APPLICABLE (forbidden quadrant) ---")
    result = check_applicability(a0=0, a1=-1, a2=0, a3=-2)
    report(result)
    print("\nSame equation, opposite wave-speed sign: a3 < 0 and f'(R0) < 0")
    print("simultaneously. No auxiliary function or wavenumber k applies.")
    print("This confirms applicability is a sign phenomenon, not a property")
    print("tied to mKdV specifically.")

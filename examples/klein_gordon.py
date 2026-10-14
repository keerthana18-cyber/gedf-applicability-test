"""
examples/klein_gordon.py

Nonlinear Klein-Gordon (phi^4-type) equation:  u_tt - u_xx - u - u^3 = 0.

Traveling-wave reduction (U(Lambda), Lambda = x + c t):
    (c^2 - 1) U'' = U + U^3
    =>  a3 = a1 = 1/(c^2-1),  a2 = a0 = 0.

No integration constant arises (already second order), so Theorem 7.1's
recovery mechanism is structurally unavailable here -- this makes the
equation a clean test of the genuine, k-independent sign-quadrant
obstruction (Proposition 8.1 of the paper).

Two regimes:
    |c| > 1  (superluminal): a3 > 0, f'(R0) = a3 > 0  -> sinh branch applies.
    |c| < 1  (subluminal):   a3 < 0, f'(R0) = a3 < 0  -> FORBIDDEN QUADRANT,
                              no auxiliary function applies, for any k.

Run:  python3 examples/klein_gordon.py
"""

import sympy as sp
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gedf_applicability_checker import check_applicability, report

Lambda = sp.symbols("Lambda")

if __name__ == "__main__":
    print("Nonlinear Klein-Gordon (phi^4-type): u_tt - u_xx - u - u^3 = 0\n")

    print("--- Case A: c^2 = 5 (superluminal, |c| > 1) ---")
    c2 = 5
    a3 = sp.Rational(1, c2 - 1)
    a1 = a3
    result = check_applicability(a0=0, a1=a1, a2=0, a3=a3)
    report(result)
    if result.applicable:
        H1 = result.H1
        k = result.k
        U = H1 / sp.sinh(k * Lambda)
        residual = sp.simplify((c2 - 1) * sp.diff(U, Lambda, 2) - (U + U**3))
        print(f"\nCandidate solution: u(x,t) = {H1} / sinh({k}(x + ct)), c = sqrt(5)")
        print(f"Symbolic check (c^2-1)U'' - (U+U^3) = {residual}  "
              f"({'verified' if residual == 0 else 'MISMATCH'})\n")

    print("--- Case B: c^2 = 1/2 (subluminal, |c| < 1) -- genuine failure ---")
    c2 = sp.Rational(1, 2)
    a3 = sp.Rational(1, c2 - 1)   # = -2
    a1 = a3
    result = check_applicability(a0=0, a1=a1, a2=0, a3=a3)
    report(result)
    print("\nAs predicted by Proposition 8.1: a3 < 0 and f'(R0) < 0 simultaneously.")
    print("No auxiliary function (cos, -sin, sinh, cosh), and no wavenumber k,")
    print("can produce an n=1 GEDF solution of this equation for |c| < 1.")
    print("This is a structural obstruction that holds for every wavenumber k,")
    print("and for every value of the integration constant (none arises here).")

"""
examples/gardner.py

Gardner equation (combined KdV-mKdV):  u_t + 6 u u_x + 6 u^2 u_x + u_xxx = 0.

Traveling-wave reduction (U(Lambda), Lambda = x - c t), integrated once:
    U'' = c U - 3 U^2 - 2 U^3 - K   =>  a3 = -2, a2 = -3, a1 = c, a0 = -K.

Uses c = -1/2, K = -1/4 (hyperbolic branch, a3 < 0, f'(R0) > 0), and
reports the resulting exact soliton solution together with the
wavenumber k recovered by the checker.

Run:  python3 examples/gardner.py
"""

import sympy as sp
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gedf_applicability_checker import check_applicability, report

Lambda = sp.symbols("Lambda")

if __name__ == "__main__":
    print("Gardner equation: u_t + 6 u u_x + 6 u^2 u_x + u_xxx = 0")
    print("Reduced master ODE (c=-1/2, K=-1/4): U'' = -U/2 - 3U^2 - 2U^3 + 1/4\n")

    a0 = sp.Rational(1, 4)
    a1 = sp.Rational(-1, 2)
    result = check_applicability(a0=a0, a1=a1, a2=-3, a3=-2)
    report(result)

    if result.applicable:
        R0, H1, k = result.R0, result.H1, result.k
        U = R0 + H1 / sp.cosh(k * Lambda)
        residual = sp.simplify(sp.diff(U, Lambda, 2) - (a1 * U - 3 * U**2 - 2 * U**3 + a0))
        print(f"\nCandidate solution: u(x,t) = {R0} + ({H1}) * sech({k}(x + t/2))")
        print(f"Symbolic check of the master ODE = {residual}  "
              f"({'verified' if residual == 0 else 'MISMATCH'})")

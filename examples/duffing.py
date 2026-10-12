"""
examples/duffing.py

Undamped, unforced Duffing oscillator:  u'' + alpha u + beta u^3 = 0.

Already in master-ODE form directly (no traveling-wave reduction
needed): a3 = -beta, a2 = 0, a1 = -alpha, a0 = 0.

Checks three physically distinct regimes:
    (a) beta < 0, alpha > 0  -- softening spring:        APPLICABLE (trig branch)
    (b) beta > 0, alpha < 0  -- hardening, inverted:      APPLICABLE (hyperbolic branch)
    (c) beta > 0, alpha > 0  -- hardening spring:         NOT APPLICABLE
                                 (forbidden sign quadrant, Proposition 8.1)

Case (c) is a genuine, physically meaningful failure: the classical
hardening Duffing oscillator has purely bounded, oscillatory orbits
around a single equilibrium and admits no homoclinic sech/csch-type
traveling wave, which the theory correctly predicts.

Run:  python3 examples/duffing.py
"""

import sympy as sp
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from gedf_applicability_checker import check_applicability, report

t = sp.symbols("t")

if __name__ == "__main__":
    print("Undamped Duffing oscillator: u'' + alpha u + beta u^3 = 0\n")

    cases = [
        ("beta=-1 (softening), alpha=1  -- APPLICABLE", -1, 1, "cos", sp.cos),
        ("beta=1 (hardening, inverted), alpha=-1  -- APPLICABLE", 1, -1, "cosh", sp.cosh),
        ("beta=1 (hardening), alpha=1  -- NOT APPLICABLE (forbidden quadrant)", 1, 1, None, None),
    ]

    for label, beta, alpha, phi_name, phi_func in cases:
        a3 = -beta
        a1 = -alpha
        print(f"--- {label} ---")
        result = check_applicability(a0=0, a1=a1, a2=0, a3=a3)
        report(result)

        if result.applicable and phi_func is not None:
            H1 = result.H1
            k = result.k
            U = H1 / phi_func(k * t)
            residual = sp.simplify(sp.diff(U, t, 2) + alpha * U + beta * U**3)
            print(f"\nCandidate solution: u(t) = {H1} / {phi_name}({k} t)")
            print(f"Symbolic check U'' + alpha U + beta U^3 = {residual}  "
                  f"({'verified' if residual == 0 else 'MISMATCH'})\n")
        elif not result.applicable:
            print("\nAs predicted: the classical hardening Duffing oscillator "
                  "admits no n=1 GEDF solution for any k.\n")

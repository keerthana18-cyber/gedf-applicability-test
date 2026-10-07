"""
gedf_applicability_checker.py
==============================

A symbolic pre-check for the applicability of the Generalized Exponential
Differential Function (GEDF) method to the cubic master ODE

    u'' = f(u) = a3*u^3 + a2*u^2 + a1*u + a0,   a3 != 0.

Implements conditions (C1)-(C4) of Theorem 4.1 (with the wavenumber k
treated as a free parameter of the traveling-wave representation), the
factorization criterion of Theorem 5.1, the integration-constant
recovery of Theorem 7.1, and the sign-quadrant obstruction of
Proposition 8.1, for a cubic master ODE obtained from a traveling-wave
reduction of a nonlinear evolution equation.

Method summary:
    The auxiliary function is taken as phi(k*Lambda) with k free, so
    condition (C4) reduces to a requirement on the SIGN of f'(R0), and
    the wavenumber is recovered afterward as
    k = sqrt(delta_phi * f'(R0)).  The amplitude carries the factor k^2,
    H1 = sqrt(eps_phi * 2 * k^2 / a3).

    A sign mismatch between a3 and f'(R0) in the "forbidden quadrant"
    a3 < 0, f'(R0) < 0 cannot be repaired by any choice of k (see
    examples/klein_gordon.py).

Usage
-----
    from gedf_applicability_checker import check_applicability, report

    result = check_applicability(a0=0.544, a1=-0.04, a2=-1.2, a3=0.5)
    report(result)

Coefficients may be Python floats/ints, SymPy Rationals, or fully
symbolic SymPy expressions.

Author's note: this implementation accompanies the paper "An Algebraic
and Geometric Applicability Theory for the Generalized Exponential
Differential Function Method" and is provided so that Theorem 4.1's
applicability test (with free wavenumber k) can be applied to a
candidate cubic master ODE without hand computation.
"""

from dataclasses import dataclass, field
from typing import Optional, List
import sympy as sp

u = sp.symbols('u')

# Auxiliary-function branch data, from Table 2 of the paper (free k):
# eps_phi (admissibility sign on a3), delta_phi (sign
# required of f'(R0), since f'(R0) = delta_phi * k^2), and a label.
BRANCHES = [
    {"name": "cos(k*Lambda)",   "eps": 1,  "delta": -1, "type": "trigonometric"},
    {"name": "-sin(k*Lambda)",  "eps": 1,  "delta": -1, "type": "trigonometric"},
    {"name": "sinh(k*Lambda)",  "eps": 1,  "delta": 1,  "type": "hyperbolic"},
    {"name": "cosh(k*Lambda)",  "eps": -1, "delta": 1,  "type": "hyperbolic"},
]


@dataclass
class GEDFResult:
    a0: object
    a1: object
    a2: object
    a3: object
    R0: object
    fR0: object                      # f(R0), should be 0 for (C1)
    fprime_R0: object                # f'(R0)
    C1: bool
    C2: bool = True                  # R0 defined via (C2); holds by construction
    matched_branch: Optional[dict] = None
    k: Optional[object] = None
    H1: Optional[object] = None
    C3: bool = False
    C4: bool = False
    applicable: bool = False
    factorization_mu2: Optional[object] = None
    diagnosis: List[str] = field(default_factory=list)
    suggested_a0: Optional[object] = None   # Theorem 7.1 recovery
    forbidden_quadrant: bool = False        # Proposition 8.1 obstruction


def _f(a0, a1, a2, a3, x):
    return a3 * x**3 + a2 * x**2 + a1 * x + a0


def _sign(expr):
    """Return +1, -1, 0, or None (undetermined, for symbolic input)."""
    expr = sp.simplify(expr)
    if expr.free_symbols:
        return None
    val = float(expr)
    if val > 0:
        return 1
    if val < 0:
        return -1
    return 0


def check_applicability(a0, a1, a2, a3, numeric_tol=1e-9) -> GEDFResult:
    """
    Check conditions (C1)-(C4) of Theorem 4.1 (free
    wavenumber k) for the cubic master ODE
    u'' = a3 u^3 + a2 u^2 + a1 u + a0.

    Returns a GEDFResult with a full diagnostic trail. Works with both
    exact SymPy values and floats; floats are compared to zero with
    `numeric_tol`.
    """
    a0, a1, a2, a3 = sp.sympify(a0), sp.sympify(a1), sp.sympify(a2), sp.sympify(a3)
    if a3 == 0:
        raise ValueError("a3 = 0: the master equation is not cubic; GEDF (n=1) does not apply.")

    def is_zero(expr):
        expr = sp.nsimplify(sp.simplify(expr))
        if expr.free_symbols:
            return sp.simplify(expr) == 0
        return abs(float(expr)) < numeric_tol

    fexpr = _f(a0, a1, a2, a3, u)
    fprime = sp.diff(fexpr, u)

    # (C2) determines R0 = -a2/(3 a3); this is the only candidate R0
    # consistent with the reduced ansatz (Theorem 3.1 / eq. 4).
    R0 = sp.simplify(-a2 / (3 * a3))

    fR0 = sp.simplify(fexpr.subs(u, R0))
    fprime_R0 = sp.simplify(fprime.subs(u, R0))

    C1_ok = is_zero(fR0)

    diagnosis = []
    if not C1_ok:
        diagnosis.append(
            f"(C1) fails: f(R0) = {fR0} != 0. The integration constant does not "
            f"currently place the equation on the applicability region."
        )

    result = GEDFResult(
        a0=a0, a1=a1, a2=a2, a3=a3, R0=R0, fR0=fR0, fprime_R0=fprime_R0,
        C1=C1_ok, diagnosis=diagnosis,
    )

    a3_sign = _sign(a3)
    fp_sign = _sign(fprime_R0)

    if fp_sign == 0:
        result.diagnosis.append(
            "(C4) fails for every auxiliary function: f'(R0) = 0, so no nonzero "
            "wavenumber k can satisfy f'(R0) = delta_phi * k^2."
        )
        result.applicable = False
        return result

    if a3_sign is None or fp_sign is None:
        result.diagnosis.append(
            "Symbolic a3 and/or f'(R0): sign cannot be determined automatically. "
            "Substitute numeric values, or check sign(a3) and sign(f'(R0)) manually "
            "against Table 2 of the paper."
        )
        return result

    # Determine which branches are admissible: need eps_phi * a3 > 0
    # (sign of a3 matches branch requirement) AND delta_phi matches the
    # sign of f'(R0) (since k^2 = delta_phi * f'(R0) must be positive).
    matched = None
    for br in BRANCHES:
        eps_ok = (br["eps"] * a3_sign) > 0
        delta_ok = (br["delta"] * fp_sign) > 0  # delta_phi and f'(R0) must share sign
        if eps_ok and delta_ok:
            matched = br
            break

    if matched is None:
        # This is the genuine, k-independent obstruction: a3 < 0 and
        # f'(R0) < 0 (Proposition 8.1 of the paper).
        result.diagnosis.append(
            f"(C3)-(C4) fail for EVERY auxiliary function and every k != 0: "
            f"sign(a3) = {a3_sign}, sign(f'(R0)) = {fp_sign}. "
            f"No wavenumber k can repair a sign mismatch between a3 and f'(R0); "
            f"this is the forbidden sign quadrant (a3 < 0, f'(R0) < 0) of "
            f"Proposition 8.1 in the paper."
        )
        result.forbidden_quadrant = (a3_sign == -1 and fp_sign == -1)
        result.applicable = False
        return result

    result.matched_branch = matched
    result.C4 = True
    result.C3 = True

    k_squared = sp.simplify(matched["delta"] * fprime_R0)
    k_val = sp.sqrt(k_squared)
    result.k = k_val

    H1_squared = sp.simplify(matched["eps"] * 2 * k_squared / a3)
    result.H1 = sp.sqrt(H1_squared)

    result.applicable = bool(C1_ok)

    if not C1_ok:
        # Theorem 7.1, Part (i): with (C2)-(C4) already satisfied via the
        # chosen k, a unique a0 = C* restores (C1), independent of a1.
        C_star = sp.simplify(-(a3 * R0**3 + a2 * R0**2 + a1 * R0))
        result.suggested_a0 = C_star
        result.diagnosis.append(
            f"(C2)-(C4) hold for the {matched['type']} branch with k = {k_val}; "
            f"Theorem 7.1 gives the required integration constant "
            f"a0 = C* = {C_star} to additionally satisfy (C1)."
        )
    else:
        # Full factorization check (Theorem 5.1), as a redundant cross-check.
        mu2 = sp.simplify(k_squared / a3) if matched["type"] == "trigonometric" \
            else sp.simplify(-k_squared / a3)
        result.factorization_mu2 = mu2
        refactored = sp.expand(a3 * (u - R0) * ((u - R0)**2 - mu2))
        if sp.simplify(refactored - fexpr) == 0:
            result.diagnosis.append(
                f"Factorization confirmed: f(u) = a3 (u-R0)[(u-R0)^2 - mu^2] "
                f"with mu^2 = {mu2}, k = {k_val}."
            )

    return result


def report(result: GEDFResult):
    print("=" * 70)
    print(f"Coefficients: a3={result.a3}, a2={result.a2}, a1={result.a1}, a0={result.a0}")
    print(f"R0 (from C2)      : {result.R0}")
    print(f"f(R0)             : {result.fR0}   -> (C1) {'holds' if result.C1 else 'FAILS'}")
    print(f"f'(R0)            : {result.fprime_R0}")
    if result.matched_branch:
        print(f"Matched branch    : {result.matched_branch['name']} "
              f"({result.matched_branch['type']}) -> (C4) holds")
        print(f"Wavenumber k      : {result.k}")
        print(f"(C3) admissibility: {'holds' if result.C3 else 'FAILS'}")
        if result.H1 is not None:
            print(f"H1                : +/- {result.H1}")
    else:
        print("(C4) FAILS for every auxiliary function, for every k != 0.")
        if result.forbidden_quadrant:
            print("--> FORBIDDEN SIGN QUADRANT (a3 < 0, f'(R0) < 0): "
                  "genuine structural obstruction, independent of k.")
    print("-" * 70)
    print(f"APPLICABLE: {result.applicable}")
    if result.diagnosis:
        print("Diagnosis:")
        for line in result.diagnosis:
            print(f"  - {line}")
    print("=" * 70)


if __name__ == "__main__":
    print("\n### Example 1: WBK-BK system (Section 8.1), R0 = 1, k = 1 special case")
    R0 = 1
    a3 = sp.Rational(1, 2) / R0**2
    a2 = sp.Rational(-3, 2) / R0
    a1 = sp.Rational(1, 2)
    a0 = sp.Rational(1, 2) * R0
    res = check_applicability(a0, a1, a2, a3)
    report(res)

    print("\n### Example 2: KD system, C = 0, general k (Section 8.2) -- APPLICABLE")
    a, m, b = 1, 0, sp.Rational(2, 5)
    a3 = sp.Rational(a**2, 2)
    a2 = sp.Rational(3*a*m - 6*b, 2)
    c = sp.Rational(64, 100)     # c = 0.64, from eq. (waveconstraint)
    a1 = c - 3*m**2
    a0 = 0
    res = check_applicability(a0, a1, a2, a3)
    report(res)

    print("\n### Example 3: nonlinear Klein-Gordon, c^2 < 1 -- FORBIDDEN QUADRANT (genuine failure)")
    a3 = sp.Rational(-1, 4)   # c^2 = 1/2  =>  a3 = a1 = 1/(c^2-1) = -2 (see examples/klein_gordon.py
    a1 = sp.Rational(-1, 4)   # for the exact worked value; this is an illustrative alternative point)
    a2 = 0
    a0 = 0
    res = check_applicability(a0, a1, a2, a3)
    report(res)

    print("\n### Example 4: fully symbolic check, R0 left free (trigonometric branch)")
    R0 = sp.symbols('R0')
    a3s = sp.symbols('a3', positive=True)
    a2s = -3 * a3s * R0
    a1s = 3 * a3s * R0**2 - sp.Rational(1, 2)   # any a1 < 3 a3 R0^2 works
    a0s = -a3s * R0**3 - a2s * R0**2 - a1s * R0
    res = check_applicability(a0s, a1s, a2s, a3s)
    report(res)

"""
Baseline: decide applicability of the n = 1 GEDF ansatz by DIRECT SUBSTITUTION.

The ansatz u = R0 + R1*phi + H1/phi is substituted into u'' = f(u); the identity
x^3 (u'' - f(u)) = 0 is a polynomial in x = phi(Lambda); equating its coefficients
gives a polynomial system in the unknowns (R0, R1, H1, K), K = k^2. The system is
solved with SymPy and a real solution with H1 != 0 and K > 0 is searched for.

Branch data (phi'' = s*K*phi,  (phi')^2 = s*K*phi^2 + c0*K):
    cos / -sin : s = -1, c0 = +1
    sinh       : s = +1, c0 = +1
    cosh       : s = +1, c0 = -1
"""
import sympy as sp

x, R0, R1, H1, K = sp.symbols("x R0 R1 H1 K")
BRANCHES = {"trig": (-1, +1), "sinh": (+1, +1), "cosh": (+1, -1)}


def coefficient_system(a0, a1, a2, a3, s, c0):
    u = R0 + R1 * x + H1 / x
    du = R1 - H1 / x**2                      # du/dphi
    d2u = du * (s * K * x) + sp.diff(du, x) * (s * K * x**2 + c0 * K)
    f = a3 * u**3 + a2 * u**2 + a1 * u + a0
    poly = sp.expand(sp.together(x**3 * (d2u - f)))
    return [c for c in sp.Poly(poly, x).all_coeffs() if c != 0]


def direct_decision(a0, a1, a2, a3):
    """Return (applicable: bool, branches: list) by solving the algebraic system."""
    found = []
    for name, (s, c0) in BRANCHES.items():
        eqs = coefficient_system(a0, a1, a2, a3, s, c0)
        sols = sp.solve(eqs, [R0, R1, H1, K], dict=True)
        for sol in sols:
            vals = [sol.get(v, v) for v in (R0, R1, H1, K)]
            if any(v.free_symbols for v in vals):
                continue
            try:
                r0, r1, h1, k2 = [sp.N(v, 30) for v in vals]
            except Exception:
                continue
            if all(abs(sp.im(v)) < 1e-20 for v in (r0, r1, h1, k2)):
                if abs(sp.re(h1)) > 1e-20 and sp.re(k2) > 1e-20:
                    found.append(name)
                    break
    return bool(found), found

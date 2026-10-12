"""Algorithm 1 of the paper in exact rational arithmetic (no SymPy required)."""
from fractions import Fraction


def decide(a0, a1, a2, a3):
    a0, a1, a2, a3 = (Fraction(v) for v in (a0, a1, a2, a3))
    if a3 == 0:
        return {"verdict": "out of scope"}
    R0 = -a2 / (3 * a3)
    t = a3 * R0
    d = a1 - 3 * R0 * t
    a0_star = -R0 * (a1 + R0 * (a2 + t))
    if d == 0 or (a3 < 0 and d < 0):
        return {"verdict": "inapplicable", "R0": R0, "d": d, "a0_star": a0_star}
    branch = "trig" if (a3 > 0 and d < 0) else ("sinh" if a3 > 0 else "cosh")
    k2 = abs(d)
    H2 = 2 * k2 / abs(a3)
    ok = (a0 == a0_star)
    return {"verdict": "applicable" if ok else "needs_a0", "branch": branch,
            "R0": R0, "d": d, "k2": k2, "H2": H2, "a0_star": a0_star}

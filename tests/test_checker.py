"""Regression tests. Run with:  python tests/test_checker.py   (no pytest required)."""
import os, random, sys
from fractions import Fraction as F
import sympy as sp

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, root)
sys.path.insert(0, os.path.join(root, "benchmarks"))
from gedf_applicability_checker import check_applicability, _sign
from decision_procedure import decide

S = lambda v: sp.Rational(v.numerator, v.denominator)
passed = 0


def check(cond, msg):
    global passed
    assert cond, "FAILED: " + msg
    passed += 1


# 1. exact decisions at extreme magnitudes (floating-point conversion would underflow)
tiny = sp.Rational(1, 10**400)
r = check_applicability(tiny, -1, 0, sp.Rational(1, 2))
check(not r.applicable and not r.C1, "a0 = 10^-400 differs from a0* = 0, so (C1) fails")
r = check_applicability(0, tiny, 0, sp.Rational(1, 2))
check(r.applicable and "sinh" in str(r.matched_branch), "f'(R0) = 10^-400 > 0 selects the sinh branch")
check(_sign(tiny) == 1 and _sign(-tiny) == -1 and _sign(sp.Integer(0)) == 0, "exact signs of tiny rationals")

# 2. the three applications of the paper
al, be, ga, R0 = sp.Rational(1), sp.Rational(1), sp.Rational(3), sp.Rational(2)    # WBK-BK
Dp = be**2 + al*ga; c = -al*R0
a3, a2, a1 = al**2/(2*Dp), 3*c*al/(2*Dp), c**2/(2*Dp)
r = check_applicability(al*R0**3*al/(2*Dp), a1, a2, a3)        # C* = alpha*R0^3, a0 = alpha*C/(2 Dp)
check(r.applicable and sp.simplify(r.H1**2 - 4*R0**2) == 0, "WBK-BK: H1 = +-2 R0 at C = C*")
r = check_applicability(0, sp.Rational(16, 25) - 0, sp.Rational(-6, 5), sp.Rational(1, 2))   # KD, a=1, m=0, b=0.4, c=0.64
check(r.applicable and sp.simplify(r.k**2 - sp.Rational(8, 25)) == 0 and sp.simplify(r.H1**2 - sp.Rational(32, 25)) == 0,
      "KD: k = 2 sqrt(2)/5, H1 = 4 sqrt(2)/5")
check(check_applicability(0, sp.Rational(1, 3), 0, sp.Rational(1, 3)).applicable, "Klein-Gordon c^2 = 4: sinh branch (a3 = a1 = 1/3)")
r = check_applicability(0, -2, 0, -2)                          # c^2 = 1/2: a3 = a1 = 1/(c^2 - 1) = -2
check((not r.applicable) and r.forbidden_quadrant, "Klein-Gordon c^2 = 1/2: forbidden sign quadrant")

# 3. agreement with Algorithm 1 on random rational input
rng = random.Random(11)
for _ in range(300):
    a1_, a2_, a3_ = [F(rng.randint(-6, 6), rng.randint(1, 4)) for _ in range(3)]
    if a3_ == 0:
        continue
    R0_ = -a2_ / (3*a3_); a0s = -R0_*(a1_ + R0_*(a2_ + a3_*R0_))
    a0_ = a0s if rng.random() < 0.6 else F(rng.randint(-5, 5), 2)
    ref = decide(a0_, a1_, a2_, a3_)
    res = check_applicability(S(a0_), S(a1_), S(a2_), S(a3_))
    check(res.applicable == (ref["verdict"] == "applicable"), f"verdict agreement for {(a0_, a1_, a2_, a3_)}")
print(f"all {passed} checks passed")

"""
Quantifier-elimination style cross-check with Z3 (nonlinear real arithmetic, nlsat).

For each branch of Table 2 the statement

  forall a0..a3 :  a3 != 0  ->  [ exists R0,R1,H1,K : K > 0, H1 != 0, coefficient system ]
                                <->  [ Delta = 0  and  sign conditions ]

is verified in two directions, where D = 3 a1 a3 - a2^2 and Delta = 27 a3^2 a0 + 2 a2^3 - 9 a1 a2 a3:
  (necessity)   no real point satisfies the coefficient system and violates the right-hand side;
  (sufficiency) for every coefficient vector satisfying the right-hand side the system has a real solution.
The coefficient system includes R1 as an unknown, so necessity also re-proves Theorem 3.1 (R1 = 0).
"""
import os, sys, time
import sympy as sp
import z3

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from direct_solver import coefficient_system, R0 as sR0, R1 as sR1, H1 as sH1, K as sK, BRANCHES

a0, a1, a2, a3 = sp.symbols("a0 a1 a2 a3")
Z = {n: z3.Real(n) for n in ("a0", "a1", "a2", "a3", "R0", "R1", "H1", "K")}


def to_z3(expr):
    expr = sp.expand(expr)
    env = {str(s): Z[str(s)] for s in expr.free_symbols}
    return eval(str(expr).replace("**", "**"), {"__builtins__": {}}, env)


D = 3*a1*a3 - a2**2
Delta = 27*a3**2*a0 + 2*a2**3 - 9*a1*a2*a3
RHS = {
    "trig": lambda: z3.And(to_z3(Delta) == 0, Z["a3"] > 0, to_z3(D) < 0),
    "sinh": lambda: z3.And(to_z3(Delta) == 0, Z["a3"] > 0, to_z3(D) > 0),
    "cosh": lambda: z3.And(to_z3(Delta) == 0, Z["a3"] < 0, to_z3(D) < 0),
}


def system(branch):
    s, c0 = BRANCHES[branch]
    eqs = coefficient_system(a0, a1, a2, a3, s, c0)
    return z3.And(*[to_z3(e) == 0 for e in eqs]), len(eqs)


def run(branch, timeout_ms=300000):
    out = {}
    sysf, n_eq = system(branch)
    base = z3.And(Z["a3"] != 0, Z["K"] > 0, Z["H1"] != 0)

    # necessity: base & system & not RHS  must be unsatisfiable
    s1 = z3.Solver(); s1.set("timeout", timeout_ms)
    s1.add(base, sysf, z3.Not(RHS[branch]()))
    t = time.perf_counter(); r = s1.check(); out["necessity"] = (str(r), time.perf_counter() - t)

    # sufficiency: for all a with RHS (and a3 != 0) there exist R0,R1,H1,K
    sa = z3.Solver(); sa.set("timeout", timeout_ms)
    ex = z3.Exists([Z["R0"], Z["R1"], Z["H1"], Z["K"]], z3.And(Z["K"] > 0, Z["H1"] != 0, sysf))
    sa.add(Z["a3"] != 0, RHS[branch](), z3.Not(ex))
    t = time.perf_counter(); r = sa.check(); out["sufficiency"] = (str(r), time.perf_counter() - t)
    return n_eq, out


if __name__ == "__main__":
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
    print("Z3", z3.get_version_string())
    for br in ("trig", "sinh", "cosh"):
        n_eq, out = run(br, limit)
        for k, (res, sec) in out.items():
            verdict = {"unsat": "VERIFIED"}.get(res, res.upper())
            print(f"{br:5s} {k:12s} ({n_eq} coefficient equations): {res:8s} {sec:8.2f} s  -> {verdict}")

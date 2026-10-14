"""
Symbolic-parameter benchmark: the three families of Section 8.

For each family the coefficients are rational functions of symbolic physical
parameters. We time
  (1) the decision procedure with symbolic coefficients (it returns the polynomials
      D and Delta of Remark 4.x and the sign conditions), and
  (2) direct substitution followed by sympy.solve of the coefficient system with the
      physical parameters left symbolic (one branch at a time, with a time limit).
"""
import multiprocessing as mp, os, sys, time
import sympy as sp

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from direct_solver import coefficient_system, R0, R1, H1, K, BRANCHES

al, be, ga, c, C, a, b, m = sp.symbols("alpha beta gamma c C a b m")
FAMILIES = {
    "WBK-BK": lambda: (lambda Dp: (al*C/Dp, c**2/Dp, 3*c*al/(2*Dp), al**2/(2*Dp)))(be**2 + al*ga),
    "KD":     lambda: (C, c - 3*m**2, (3*a*m - 6*b)/2, a**2/2),
    "KG":     lambda: (0, 1/(c**2 - 1), 0, 1/(c**2 - 1)),
}


def symbolic_decision(a0, a1, a2, a3):
    D = sp.factor(3*a1*a3 - a2**2)
    Delta = sp.factor(27*a3**2*a0 + 2*a2**3 - 9*a1*a2*a3)
    return D, Delta


def _direct(queue, coeffs, branch, with_C):
    s, c0 = BRANCHES[branch]
    eqs = coefficient_system(*coeffs, s, c0)
    unknowns = [R0, R1, H1, K] + ([C] if with_C else [])
    t0 = time.perf_counter()
    sols = sp.solve(eqs, unknowns, dict=True)
    nontrivial = [t for t in sols if t.get(H1, 1) != 0 and t.get(K, 1) != 0]
    queue.put((len(sols), len(nontrivial), time.perf_counter() - t0))


def timed_direct(coeffs, branch, limit, with_C=False):
    q = mp.Queue()
    p = mp.Process(target=_direct, args=(q, coeffs, branch, with_C))
    p.start(); p.join(limit)
    if p.is_alive():
        p.terminate(); p.join()
        return None
    return q.get()


if __name__ == "__main__":
    limit = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
    for name, make in FAMILIES.items():
        coeffs = make()
        t0 = time.perf_counter()
        D, Delta = symbolic_decision(*coeffs)
        t_dec = time.perf_counter() - t0
        print(f"\n== {name}")
        print(f"   decision procedure: {t_dec*1e3:8.1f} ms")
        print(f"   D     = {D}")
        print(f"   Delta = {Delta}")
        variants = [(False, "parameters symbolic")]
        if name != "KG":
            variants.append((True, "integration constant C also unknown"))
        for with_C, label in variants:
            for br in BRANCHES:
                res = timed_direct(coeffs, br, limit, with_C)
                if res is None:
                    print(f"   direct solve ({label}), branch {br:5s}: no result within {limit:.0f} s")
                else:
                    print(f"   direct solve ({label}), branch {br:5s}: {res[2]:6.2f} s, "
                          f"{res[0]} solution sets, {res[1]} with H1 != 0 and K != 0")

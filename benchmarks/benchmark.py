"""
Benchmark: decision procedure (Algorithm 1) versus direct substitution.

For each random rational coefficient vector the verdict of
  (1) Algorithm 1 in exact arithmetic            (decision_procedure.decide),
  (2) the SymPy implementation in the repository (check_applicability),
  (3) direct substitution and solving of the coefficient system (direct_solver)
is computed and timed. Output: benchmark_results.csv and a printed summary.
"""
import csv, platform, random, statistics, sys, time, os
from fractions import Fraction as F
import sympy as sp

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
sys.path.insert(0, os.path.dirname(here))
from decision_procedure import decide
from direct_solver import direct_decision
from gedf_applicability_checker import check_applicability


def rnd(rng, lo=-6, hi=6, den=4):
    return F(rng.randint(lo, hi), rng.randint(1, den))


def main(n_per_set=150, seed=2026):
    rng = random.Random(seed)
    rows = []
    for family in ("a0=a0*", "a0 random"):
        count = 0
        while count < n_per_set:
            a1, a2, a3 = rnd(rng), rnd(rng), rnd(rng)
            if a3 == 0:
                continue
            R0 = -a2 / (3 * a3)
            a0s = -R0 * (a1 + R0 * (a2 + a3 * R0))
            a0 = a0s if family == "a0=a0*" else F(rng.randint(-5, 5), 2)
            S = lambda v: sp.Rational(v.numerator, v.denominator)

            t0 = time.perf_counter(); r1 = decide(a0, a1, a2, a3); t1 = time.perf_counter() - t0
            t0 = time.perf_counter(); r2 = check_applicability(S(a0), S(a1), S(a2), S(a3)); t2 = time.perf_counter() - t0
            t0 = time.perf_counter(); ok3, br3 = direct_decision(S(a0), S(a1), S(a2), S(a3)); t3 = time.perf_counter() - t0

            ok1 = r1["verdict"] == "applicable"
            br1 = r1.get("branch") if ok1 else None
            agree = (ok1 == bool(r2.applicable) == ok3) and (not ok1 or br3 == [br1])
            rows.append(dict(family=family, a0=a0, a1=a1, a2=a2, a3=a3, applicable=ok1,
                             branch=br1 or "", agree=agree,
                             t_algorithm=t1, t_checker=t2, t_direct=t3))
            count += 1
    with open(os.path.join(here, "benchmark_results.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    N = len(rows)
    print(f"instances: {N}   verdict/branch agreement: {sum(r['agree'] for r in rows)}/{N}")
    print(f"applicable: {sum(r['applicable'] for r in rows)}   inapplicable: {sum(not r['applicable'] for r in rows)}")
    for b in ("trig", "sinh", "cosh"):
        print(f"  branch {b}: {sum(r['branch']==b for r in rows)}")
    for key, label in (("t_algorithm", "Algorithm 1 (exact Fractions)"),
                       ("t_checker", "repository checker (SymPy)"),
                       ("t_direct", "direct substitution (SymPy solve)")):
        v = [r[key] for r in rows]
        print(f"{label:38s} median {statistics.median(v)*1e3:10.4f} ms   mean {statistics.mean(v)*1e3:10.4f} ms   max {max(v)*1e3:10.2f} ms")
    med = statistics.median([r['t_direct'] / r['t_algorithm'] for r in rows])
    print(f"median speed-up of Algorithm 1 over direct substitution: {med:,.0f}x")
    bad = [r for r in rows if not r['agree']]
    for r in bad[:5]:
        print("MISMATCH", r)


class Tee:
    def __init__(self, *streams):
        self.streams = streams
    def write(self, text):
        for st in self.streams:
            st.write(text)
    def flush(self):
        for st in self.streams:
            st.flush()


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    with open(os.path.join(here, "benchmark_summary.txt"), "w") as fh:
        sys.stdout = Tee(sys.__stdout__, fh)
        print(f"Python {platform.python_version()}, SymPy {sp.__version__}, {platform.platform()}, "
              f"{platform.processor() or platform.machine()}, {os.cpu_count()} logical CPUs")
        main(n)
        sys.stdout = sys.__stdout__

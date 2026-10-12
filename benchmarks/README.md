# Benchmarks and verification

Environment of the reported run: Python 3.12.3, SymPy 1.14.0, Z3 5.1.0,
one logical CPU of an Intel Xeon at 2.1 GHz (Linux). Timings depend on the machine.
Rerun with `python benchmarks/benchmark.py 250` (about ten minutes).

## 1. Random rational instances (`benchmark.py`)

500 coefficient vectors with `a1, a2, a3` of the form p/q (|p| <= 6, 1 <= q <= 4,
`a3 != 0`): 250 with `a0 = a0*` (the required constant coefficient) and 250
with a random `a0`. Three decision methods are compared on every instance:

1. Algorithm 1 in exact rational arithmetic (`decision_procedure.py`),
2. the SymPy implementation (`gedf_applicability_checker.py`),
3. direct substitution of the ansatz, with `sympy.solve` applied to the
   coefficient system in the unknowns `R0, R1, H1, k^2` (`direct_solver.py`).

| quantity | result |
|---|---|
| instances | 500 (205 applicable: 80 trigonometric, 41 sinh, 84 cosh; 295 inapplicable) |
| agreement of verdicts and branches, all three methods | 500 / 500 |
| Algorithm 1 | median 0.030 ms, mean 0.033 ms, maximum 0.09 ms |
| SymPy implementation | median 4.28 ms, mean 4.74 ms, maximum 60.1 ms |
| direct substitution with `sympy.solve` | median 415.6 ms, mean 1160.1 ms, maximum 3588 ms |
| median per-instance ratio, direct / Algorithm 1 | about 15 000 |

The per-instance results are in `benchmark_results.csv`, the summary in
`benchmark_summary.txt`. The direct baseline is the generic `sympy.solve`, not a
specialized solver, so the ratio measures this baseline only.

## 2. Symbolic parameters (`benchmark_symbolic.py`, output in `benchmark_symbolic_output.txt`)

For the generalized WBK-BK system, the KD system and the nonlinear Klein-Gordon
equation, the decision procedure returns, in 4 to 20 ms, the polynomials
`D = 3 a1 a3 - a2^2` and `Delta = 27 a3^2 a0 + 2 a2^3 - 9 a1 a2 a3` that define the
applicability region. Direct solving with the physical parameters left symbolic
takes 0.3 to 1.3 s per branch and returns no solution for WBK-BK and KD (the system
is inconsistent for generic parameters, so the parametric solutions are missed).
When the integration constant `C` is added to the unknowns, it takes 0.6 to 2.6 s and
returns five solution sets, two of them nontrivial, whose reality and positivity
conditions are not stated. For Klein-Gordon it returns seven sets per branch whose
admissibility depends on the sign of `c^2 - 1`, which is left to the analyst.

## 3. Verification with Z3 (`qe_verification.py`, output in `qe_verification_output.txt`)

For each branch of Table 2, the nonlinear real arithmetic procedure of Z3 decides, for
all real `a0, ..., a3` with `a3 != 0`:

* necessity: no real solution of the coefficient system with `K > 0`, `H1 != 0`
  violates `Delta = 0` together with the sign conditions of the branch;
* sufficiency: whenever `Delta = 0` and the sign conditions hold, the coefficient
  system has a real solution with `K > 0`, `H1 != 0`.

All six statements are valid (the negated statements are unsatisfiable, each in
well under one second). Three negative controls (a wrong sign condition, `H1 != 0`
dropped, `Delta = 0` dropped) produce counterexamples, so the checks are not vacuous.
The coefficient system keeps `R1` as an unknown, so necessity also re-proves that
`R1 = 0`. This is a solver result, not a machine-checked proof.

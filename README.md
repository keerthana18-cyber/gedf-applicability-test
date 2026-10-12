# GEDF Applicability Checker

This repository accompanies the manuscript "An Algebraic and Geometric
Applicability Theory for the Generalized Exponential Differential
Function Method" (N. Keerthana, V. Umapathi, N. Annapoorani).

A symbolic pre-check for the applicability of the Generalized Exponential
Differential Function (GEDF) method to the cubic master ODE

    u'' = a3*u^3 + a2*u^2 + a1*u + a0,   a3 != 0,

implementing the necessary and sufficient applicability conditions
(C1)-(C4), with the wavenumber `k` treated as a free parameter of the
traveling-wave representation, the factorization criterion, the
integration-constant recovery mechanism, and the sign-quadrant
non-universality obstruction, all from the accompanying paper.

Current release: **1.1.0**. See `CHANGELOG.md`, `examples/README.md` and `benchmarks/README.md`.

## Requirements

- Python 3.8+
- SymPy (`pip install sympy`)
- Z3 (`pip install z3-solver`), only for `benchmarks/qe_verification.py`

Decisions are exact for rational and algebraic coefficients. Floating-point
input is rationalized with tolerance 1e-9, so a certified decision requires
exact input; signs that cannot be certified, and symbolic signs, are reported
as undetermined and are never guessed.

## Usage

```python
from gedf_applicability_checker import check_applicability, report

# KD system example, C = 0, general wavenumber k
result = check_applicability(a0=0, a1=0.64, a2=-1.2, a3=0.5)
report(result)
```

Coefficients may be Python floats, ints, `Fraction`/`Rational`, or fully
symbolic SymPy expressions (see Example 4 in `gedf_applicability_checker.py`
for a symbolic check with `R0` left free).

Running the script directly reproduces four illustrative examples,
including the WBK-BK system, the KD system under `C=0`, and a
forbidden-sign-quadrant case:

```bash
python3 gedf_applicability_checker.py
```

## What it checks

Given `(a0, a1, a2, a3)`:

1. Computes `R0 = -a2/(3*a3)` (from condition C2).
2. Checks `f(R0) = 0` (condition C1).
3. Checks the **sign** of `f'(R0)` against each auxiliary function's
   required sign (condition C4), so that the wavenumber is not
   restricted to a fixed value.
4. Checks the admissibility sign on `a3` (condition C3) and recovers
   the wavenumber `k = sqrt(delta_phi * f'(R0))` and the amplitude `H1`.
5. If (C1) fails but (C2)-(C4) hold, reports the unique integration
   constant `a0 = C*` that restores applicability.
6. If **no** auxiliary-function branch's sign requirements can be met
   for any wavenumber, reports this as the forbidden sign quadrant
   (`a3 < 0`, `f'(R0) < 0`): a genuine structural obstruction, independent of
   the wavenumber.
7. Cross-checks the result against the factorization
   `f(u) = a3(u-R0)[(u-R0)^2 - mu^2]` (trigonometric/cosh branches) or
   `f(u) = a3(u-R0)[(u-R0)^2 + nu^2]` (sinh branch) when applicable.

## Two kinds of "not applicable"

This checker distinguishes two genuinely different reasons a master
equation can be reported non-applicable:

1. **Out of scope** (`a3 = 0`): the master equation is not cubic, so
   the `n=1` GEDF ansatz does not apply by construction. This says
   nothing about conditions (C1)-(C4).
2. **Forbidden sign quadrant** (`a3 < 0` and `f'(R0) < 0`): the master
   equation is genuinely cubic, but no auxiliary function admits a
   solution for *any* wavenumber, amplitude, or integration constant.
   This is a real structural obstruction, independent of the
   wavenumber.

See `examples/README.md` for worked illustrations of both.

## Additional examples

The `examples/` directory contains standalone, runnable checks against
several nonlinear evolution equations and nonlinear oscillators from
the traveling-wave literature, beyond the three detailed applications
(WBK-BK, KD, nonlinear Klein-Gordon) presented in the paper:

- `examples/mkdv.py`: modified KdV equation, both wave-speed signs
  (applicable for `c>0`, forbidden quadrant for `c<0`)
- `examples/gardner.py`: Gardner equation (hyperbolic branch, bounded
  sech soliton)
- `examples/duffing.py`: undamped Duffing oscillator: softening and
  inverted-hardening cases (applicable), classical hardening case
  (forbidden quadrant)
- `examples/klein_gordon.py`: nonlinear Klein-Gordon equation, both
  wave-speed regimes (applicable for `|c|>1`, forbidden quadrant for
  `|c|<1`, matching Section 8.3 of the paper)
- `examples/kdv_nonapplicable.py`: plain KdV equation, demonstrating
  correct **out-of-scope rejection**: its traveling-wave reduction
  gives `a3 = 0`

Each example independently verifies its candidate solution by direct
symbolic substitution back into the governing ODE. Run any of them
with, e.g.,

```bash
python3 examples/mkdv.py
```

## Verification and benchmarks

- `tests/test_checker.py` runs the regression tests (`python tests/test_checker.py`).
- `benchmarks/decision_procedure.py` is Algorithm 1 of the paper in exact
  rational arithmetic and does not need SymPy.
- `benchmarks/benchmark.py` compares Algorithm 1, the checker, and direct
  substitution followed by `sympy.solve` on random rational coefficient vectors.
- `benchmarks/benchmark_symbolic.py` does the same for the families of Section 8
  with symbolic physical parameters.
- `benchmarks/qe_verification.py` verifies, with the nonlinear real arithmetic
  solver of Z3, that the sign conditions of Table 2 are equivalent to the
  coefficient system for every real coefficient vector.

See `benchmarks/README.md` for the results and the environment.

## License

MIT License, see `LICENSE`.

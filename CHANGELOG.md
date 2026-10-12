# Changelog

## Version 1.1.0

- Sign and zero tests are exact for rational and algebraic input. Version
  1.0.0 converted values to floating point before testing the sign, which
  misclassified rational coefficients of magnitude below about 1e-308
  (a nonzero value was read as zero). Floating-point input is rationalized
  with tolerance 1e-9, and an uncertifiable sign is reported as undetermined.
- Added `benchmarks/`: Algorithm 1 in exact arithmetic, a comparison with
  direct substitution on 500 random rational instances, a symbolic-parameter
  benchmark for the three families of Section 8, and a verification of the
  branch conditions with the nonlinear real arithmetic solver of Z3.
- Added `tests/test_checker.py` (regression tests, no extra dependencies).

## Version 1.0.0

Initial public release.

- Symbolic applicability test for the cubic master ODE
  `u'' = a3*u^3 + a2*u^2 + a1*u + a0` (`a3 != 0`), implementing
  conditions (C1)-(C4) with the wavenumber `k` treated as a free
  parameter. Condition (C4) is a sign condition on `f'(R0)`, and the
  wavenumber is recovered as `k = sqrt(delta_phi * f'(R0))`.
- Amplitude `H1 = sqrt(eps_phi * 2 * k^2 / a3)`.
- Factorization cross-check, integration-constant recovery, and
  detection of the forbidden sign quadrant (`a3 < 0`, `f'(R0) < 0`),
  reported through the `forbidden_quadrant` flag of `GEDFResult`.
  This is kept distinct from the out-of-scope case `a3 = 0`.
- Worked examples in `examples/`: modified KdV, Gardner, nonlinear
  Klein-Gordon, Duffing oscillator, and a non-applicable case (plain
  KdV). Every candidate solution is verified by direct symbolic
  substitution into the governing ODE.

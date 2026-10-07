# Changelog

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

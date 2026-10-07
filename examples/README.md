# Examples

Standalone, runnable checks of the GEDF applicability test (with the
wavenumber `k` treated as a free parameter) against several nonlinear
evolution equations and nonlinear oscillators from the traveling-wave
literature, beyond the three detailed applications (WBK-BK, KD,
nonlinear Klein-Gordon) presented in the accompanying paper.

Each script is self-contained: it performs (or states) the
traveling-wave reduction to the master ODE, calls
`check_applicability` from the parent `gedf_applicability_checker.py`,
and independently verifies the resulting candidate solution by direct
symbolic substitution back into the governing ODE.

| File | Equation | Outcome |
|---|---|---|
| `mkdv.py` | Modified KdV, wave speed `c > 0` | Applicable: hyperbolic branch |
| `mkdv.py` | Modified KdV, wave speed `c < 0` | **Not applicable**: forbidden sign quadrant |
| `gardner.py` | Gardner equation | Applicable: hyperbolic branch |
| `duffing.py` | Duffing oscillator, softening spring | Applicable: trigonometric branch |
| `duffing.py` | Duffing oscillator, hardening (inverted) | Applicable: hyperbolic branch |
| `duffing.py` | Duffing oscillator, classical hardening spring | **Not applicable**: forbidden sign quadrant |
| `klein_gordon.py` | Nonlinear Klein-Gordon, `\|c\| > 1` | Applicable: hyperbolic branch |
| `klein_gordon.py` | Nonlinear Klein-Gordon, `\|c\| < 1` | **Not applicable**: forbidden sign quadrant, no integration constant available to recover |
| `kdv_nonapplicable.py` | Plain KdV | **Not applicable**: `a3 = 0`, outside scope by construction (not a (C1)-(C4) failure) |

## Running an example

```bash
python3 mkdv.py
```

(run from inside this directory, or adjust `sys.path` as shown at the
top of each file if running from elsewhere.)

## Two kinds of "not applicable"

This set of examples distinguishes two genuinely different reasons an
equation can be reported non-applicable:

1. **Out of scope** (`kdv_nonapplicable.py`): the master ODE is not
   cubic (`a3 = 0`), so the `n=1` GEDF ansatz does not apply by
   construction: this says nothing about conditions (C1)-(C4).
2. **Forbidden sign quadrant** (`klein_gordon.py`, `duffing.py`,
   `mkdv.py`, subluminal/hardening/negative-speed cases): the master
   ODE is genuinely cubic, but `a3` and `f'(R0)` have signs that no
   auxiliary function admits, for *any* wavenumber `k` and *any*
   integration constant. This is a real structural obstruction
   (Proposition 8.1 of the paper), independent of the wavenumber.

The `mkdv.py` and `duffing.py` examples are particularly instructive:
each shows the *same* equation succeeding or failing depending only
on the sign of one physical parameter (wave speed, or spring
stiffness/orientation), directly illustrating that applicability is
governed by sign structure rather than by which named equation is
being solved.

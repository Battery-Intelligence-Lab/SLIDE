# Phase 9A audit-debt closure — 2026-07-10

This report closes AUD-1 through AUD-4 from `PLAN.md` §2.2. It distinguishes
pre-registered accuracy claims from empirical regression sentinels and records
the exact scope of two waivers. No production model equation was changed.

## AUD-1 — P2-G1 digit-parity limits

Phase 9A explicitly required both the maximum pack-voltage drift and maximum
branch-current drift to be enforced at `1e-12`. The original §5.2 text had
explicitly registered voltage and mapped-state bands at `1e-12`; it did not
explicitly mention current, so the current limit is described as the Phase-9A
audit-aligned digit-parity scale rather than retroactively mislabelled as an
original registration.

`tests/parity/P2G1_pack_test.cpp` now has no hidden 100-fold slack:

- maximum voltage drift must be `<= 1e-12 V`;
- maximum branch-current drift must be `<= 1e-12 A`;
- the existing relative state gate, with its absolute floor, is unchanged.

The deterministic 300-step 3s2p gate passed in both configurations. The local
revalidation observed:

| Build | max ΔV | voltage-band use | max ΔI | current-band use | mapped-state band use |
|---|---:|---:|---:|---:|---:|
| Debug | 5.3291e-15 V | 0.533% | 3.6948e-13 A | 36.948% | 0.739% |
| Release (`-Ofast`) | 5.3291e-15 V | 0.533% | 3.1264e-13 A | 31.264% | 0.718% |

The worst current headroom is therefore only 2.71×. Tightening the assertion is
substantive. The earlier recorded `2.95e-13 A` remains a valid historical run,
but the table above is the current cross-configuration evidence.

## AUD-2 — Phase-5 band provenance and scale

Commit history establishes that commit `09e8ec5` introduced the Phase-5
implementation, its tests, the PLAN completion claim, and all four numerical
bands together. They were not registered before the decisive implementation
run and must not be presented as independent accuracy guarantees.

The relevant solver contracts differ:

- legacy CV calls `setVoltage_iterative`, whose early-success condition is only
  `abs(Vset - V) < 1e-6 V` and which uses at most 50 false-position iterations;
- v4 CV iterates its analytical tangent until the current update is
  `<= 1e-10 A`, with at most 12 iterations;
- the parity case runs 30 s of CC and 30 s of CV with `dt = 1 s`;
- it compares only the final applied CV-interval current, not the maximum
  current difference over the trajectory.

Consequently, `0.2 µV` and `20 µA` are empirical cross-implementation
regression envelopes. Their ratio is `0.01 Ω`, a plausible differential-voltage
scale, but the test does not establish a trajectory-wide resistance bound and
the legacy contract alone cannot guarantee either value.

The `0.2 µAh` band has a dimensionally principled scale:

```text
20 µA × 30 s / 3600 s h⁻¹ = 0.1667 µAh
```

Rounding that to `0.2 µAh` leaves 1.2× headroom. It is not logically implied by
the final-current assertion, so it remains an independent empirical regression
gate on integrated trajectory drift.

The event locator bisects a step 45 times. For the tested one-second outer step,
the final time bracket is at most

```text
1 s / 2^45 = 2.8422e-14 s.
```

That makes a `2e-12 V` event residual a sensible floating-point/order-of-
magnitude sentinel for this smooth case, but there is no registered global
`|dV/dt|` bound from which it follows, and its provenance is still post-hoc.

The current Release reproduction remains inside the unchanged empirical bands:

| Quantity | Observed | Gate | Gate use |
|---|---:|---:|---:|
| final voltage difference | 0.101674 µV | 0.2 µV | 50.84% |
| final applied-current difference | 12.479451 µA | 20 µA | 62.40% |
| charge-throughput difference | 0.116939 µAh | 0.2 µAh | 58.47% |

These values are regression evidence only. The test comment now records that
provenance so a future reader cannot mistake tight agreement for an a priori
accuracy certificate.

## AUD-4 — CVODE decision

The mandatory CVODE Phase-3 arbiter is waived by D-26 for v4.0. The reason is
narrow: with frozen coefficients and piecewise-constant flux, each modal
diffusion equation is diagonal and has an exact analytic solution. A numerical
CVODE solve of that same subflow would introduce tolerance and convergence
error into a problem with a zero-reference-error closed form.

The P3-G1 oracle was strengthened before accepting the waiver. It now evaluates

```text
z(h) = exp(rate h) z(0) + expm1(rate h) forcing / rate
```

in `long double`, with the exact `rate = 0` limit, rather than copying
production's `abs(x) < 1e-7` Taylor branch and operation order. Production-
compiled A/B values are shared deliberately: they define the ODE whose time
propagator is under test. P1-G3 separately arbitrates spectral compilation.
P3-G1 covers both electrodes, nonzero forcing, all modes, and
`nch = {5, 8, 12}` at a relative/absolute `2e-12` scale in Debug and Release.
MSVC may give `long double` the same width as `double`; the portable source of
independence is therefore the distinct analytic branch and operation order,
not an assumption of extra precision.

The waiver is not broader than that evidence. The originally described full
1C/current-step, terminal-voltage CVODE comparison did not run. P3-G1 does not
validate a future generic RHS adapter, SUNDIALS interoperability, nonlinear
full-cell integration, or slow-physics splitting. Those have separate current
tests, and a shipped CVODE adapter would require its own optional-dependency and
segment-boundary gates.

## AUD-3 / Q11 — timing-protocol waiver

D-27 retroactively waives the quiet-host/Volkan-operator condition for the
already-completed v4.0 PAY-1, PAY-2, and PAY-4 checkpoints. It does not convert
the development-host runs into portable performance guarantees. Busy-host
interference can favor either numerator or denominator, so the prior statement
that a quiet machine could only improve the ratios was removed.

The phase decisions remain supportable because their correctness, equal-work,
allocation, factorisation, and memory gates passed, and PAY-1/2 alternate order
and report conservative within-run ratios. Every PAY number remains labelled
qualified development-host evidence. Any future unqualified claim requires a
named quiet host, committed raw output, repeated timings for both tools, and
load/thermal stability criteria.

The previously ignored PAY-4 JSON is now preserved at
`benchmark/results/pay4-development-host-2026-07-10.json`. It also makes an
important limitation auditable: each liionpack case has one comparator timing,
whereas SLIDE and the single-cell PyBaMM comparison record ranges. PAY-1/2 raw
stdout was not preserved, which is stated rather than reconstructed.

## Validation record

After an explicit adversarial precheck, only the registered short cases were
run. The first direct-binary attempt used the repository root and failed before
stepping because legacy data paths require CTest's configured working
directory; that harness-location failure is not counted as a numerical result.
Rerunning in the configured working directory produced:

- Debug targeted CTest: 3/3 passed (`P2G1_pack`, `core_Experiment`,
  `core_ExponentialModal`), 0.48 s total;
- Release targeted CTest: 3/3 passed, 0.69 s total;
- verbose Release Phase-5 reproduction: 14/14 assertions;
- verbose Debug and Release P2-G1: 3312/3312 assertions in each build.
- post-edit full Debug CTest: 49/49 passed, 27.57 s total;
- post-edit full Release CTest: 49/49 passed, 15.65 s total.

No benchmark was rerun: a run on the same busy host could not satisfy the
waived quiet-host/operator protocol and would add weak evidence rather than
resolve it.

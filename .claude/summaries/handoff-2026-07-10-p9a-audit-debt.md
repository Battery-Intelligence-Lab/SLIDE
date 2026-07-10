# Handoff — Phase 9A audit debt complete (2026-07-10)

## Outcome

Phase 9A is complete. AUD-1 through AUD-4 are resolved without changing any
production model equation. The next authoritative item is Phase 9B's systematic
adversarial bug-hunt.

## Changes

- P2-G1 voltage and branch-current assertions are now `1e-12`, not `1e-10`.
  Current Debug/Release maxima are `5.3291e-15 V` and, worst case,
  `3.6948e-13 A`.
- The Phase-5 comment and audit report classify the same-commit `0.2 µV`,
  `20 µA`, `0.2 µAh`, and `2e-12 V` thresholds as empirical regression
  sentinels. Release reproduction is `0.101674 µV`, `12.479451 µA`, and
  `0.116939 µAh`.
- P3-G1's expected modal state now uses an independent long-double analytic
  evaluation: `exp(rate*h)*z0 + expm1(rate*h)/rate*forcing`, with the exact
  zero-rate limit. It does not copy production's small-x Taylor branch.
- D-26 waives mandatory CVODE narrowly for that exact diagonal subflow. The
  unrun full 1C/current-step voltage comparison is not claimed.
- D-27 resolves Q11 by waiving the quiet-host/operator condition for already
  completed v4.0 PAY-1/2/4 runs. All remain qualified development-host
  evidence; busy-host noise is not assumed to have a favourable direction.
- Raw PAY-4 evidence moved from ignored `.cache` into
  `benchmark/results/pay4-development-host-2026-07-10.json`.

## Validation

- Debug targeted CTest: 3/3 passed in 0.48 s.
- Release targeted CTest: 3/3 passed in 0.69 s.
- P2-G1 verbose: 3312/3312 assertions in each configuration.
- Release Phase-5 verbose: 14/14 assertions.
- Final full Debug suite: 49/49 passed in 27.57 s.
- Final full Release suite: 49/49 passed in 15.65 s.
- The first direct-binary attempt ran from the wrong working directory and
  failed before simulation while opening legacy CSV data. It was discarded and
  rerun through CTest's configured working directory.

Detailed evidence and derivations are in
`.claude/reports/p9a-audit-debt-2026-07-10.md`.

## Next

Execute Phase 9B subsystem-by-subsystem. Every confirmed defect needs a short
failing regression before its fix; every refuted candidate belongs in the bug
ledger. Keep simulations short and perform the adversarial precheck before each.

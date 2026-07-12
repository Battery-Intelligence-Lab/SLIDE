# Handoff — M0.3 P9-G3 exact Status coverage (2026-07-12)
- Goal: prove every active `src/core` Status failure arm is measured or structurally excepted.
- Gate: Linux Clang/LLVM 18.1.3 optional-off PASS — 366 lexical, 336 active, 329 covered, 7/10 excepted, 30 inactive, 0 uncovered/unmapped.
- Evidence: every one of 52 registered tests supplied a fresh build-local profile; exports are per binary and a whole-archive anchor proves all production mappings.
- Exceptions: five defensive-only and two platform-specific, pinned by exact site plus source/context hashes.
- Hardening: bounded I/O/allocation transactions and deterministic numerical, pack, executor, parser, recorder, and corruption seams cover every reachable arm.
- Validation: coverage, native Debug, and fast-math Release passed 52/52; ASan+UBSan passed 51/51 short tests; TSan passed 3/3.
- Artifacts: `.claude/reports/p9g3-status-coverage-2026-07-12.{md,json}`; registry: `tests/coverage/status_failure_exceptions.json`.
- Caveat/debt: hosted workflow has not run; MC-1 oversized tooling/source/tests are assigned to M0.7/M0.8.
- Next: proceed immediately to M0.4 P9-G4 bug-ledger closure; no boundary pause.

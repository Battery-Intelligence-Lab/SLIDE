# Handoff — M0.4 P9-G4 adversarial bug ledger (2026-07-12)
- Goal: close every Phase-9B defect with short red evidence, a fix, and durable classification.
- Ledger: P9-B01..B50 are FIXED; 23 candidates are retained as REFUTED with proofs or tests.
- Final bugs: schedule/extent/allocation atomicity, LLVM tool selection, FTZ conductance, resistor publication, spectrum multiplicity, read amplification, I/O allocation, and source stepping.
- New hardening: integer-bit strictly-positive finite classification survives fast-math and explicit FTZ.
- Mutation reds: source stepping 870/871, FTZ conductance 874/875, compressed allocation budget 57/60.
- Validation: native Debug/Release 52/52; Linux Release PackSolver 1/1; focused ASan+UBSan 2/2.
- Exact coverage refresh: 52/52 profiles; 329 measured + 7/10 pinned exceptions, 0 failed, 30 inactive.
- Artifacts: `p9b-bug-ledger-2026-07-10.md` and `p9g4-status-coverage-refresh-2026-07-12.{md,json}`.
- Next: proceed immediately to M0.5/9C-1 single-source modal and observable physics extraction.

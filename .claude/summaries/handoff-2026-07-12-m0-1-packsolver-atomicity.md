# Handoff — M0.1 PackSolver allocation atomicity (2026-07-12)
- Goal: settle the two-file WIP and close the late-allocation reconfiguration defect.
- Changed: build every solver/executor/workspace/vector candidate locally; translate allocation/length failures; publish only with statically no-throw moves.
- Red evidence: an explicit premature-publication mutation failed the focused sentinel solve before the injected fault.
- Focused tests: Debug and fast-math Release each pass ParserAllocation 946, PackSolver 704, PackStepper 208, and P2-G1 allocation 14 assertions.
- Full tests: Debug 51/51 and fast-math Release 51/51 pass; parity and allocation gates remain green.
- Interpretation: P9-B37 is a real cold-path correctness fix; successful numerical behavior is unchanged.
- Memory: P9B ledger, CHANGELOG, LESSONS, PLAN §8, and develop/TODO are updated.
- Next: M0.2; first isolate build outputs and prove sanitizer flags reach production `slide_core`, then run full ASan+UBSan and focused WSL TSan.

# Handoff — M0.2 P9-G1 sanitizers (2026-07-12)
- Goal: make ASan+UBSan/TSan lanes real, isolated, reproducible, and green.
- Infrastructure: project-wide pre-target sanitizer flags, fatal incompatibility checks, build-local top-level executables, three explicit TSan labels, pinned Ubuntu workflow.
- Red evidence: full ASan found P9-B38, a six-byte over-read in the embedded-NUL netlist fixture; fixed with array-derived extent.
- ASan+UBSan: flags/symbols/51 build-local paths proven; final WSL Clang 18 suite 51/51 with leak/UB abort settings.
- TSan: flags/symbols proven; PackStepper/ThreadPool/AsyncRecorder 3/3. Allocation-counter test excluded because global new/delete conflict with TSan interceptors; it passes ASan.
- Regression: native Debug/Release 51/51; top-level core-only and nested-consumer smokes 1/1; workflow YAML parses.
- Hosted workflow is committed but has not run until pushed; do not claim hosted green.
- Next: M0.3 llvm-cov Status-failure-site measurement; instrument production `slide_core`, not only test wrappers.

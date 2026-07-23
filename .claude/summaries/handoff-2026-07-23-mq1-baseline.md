# Handoff — MQ.1 three-lane baseline restored (2026-07-23)
- Goal: re-establish the M0.10 Debug + fast-math Release + CUDA baseline before MQ changes.
- PASS: fresh Debug/ThinLTO, Release/IPO-off, and CUDA with host-C++ ThinLTO each pass 58/58.
- CUDA witness: 433,671 assertions / 4 cases; every repo-root SEI witness passes 57 / 3.
- FALSIFIED: `cl.exe` on `PATH` alone is insufficient; fresh CUDA needs full x64 `vcvars64.bat`.
- No source/test changed; no timing, sanitizer, coverage, hosted-CI, package, or cross-platform claim.
- New MQ.2 inputs: deprecated redundant `-Ofast` and legacy `Deep_ptr` narrowing diagnostics.
- Evidence: `.claude/reports/mq1-baseline-{preregistration,validation}-2026-07-23.md`.
- NEXT: MQ.2 — disposition all 68 survivors plus the two new build-policy findings.

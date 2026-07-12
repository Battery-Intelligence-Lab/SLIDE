# Handoff — M0.5 / 9C-1 single-source physics (2026-07-12)
- PC-10 now routes CPU, CUDA, Dual, observables, curve interpolation, and spectral diffusion through `SpmScalarKernels.hpp`.
- Independent raw `SpectralDiffusionLegacy.hpp` remains the numerical oracle and is structurally isolated.
- Capture-host exact traces match in Debug, fast-math Release, ThinLTO, and CUDA; exact hashes are explicit opt-in, not portable oracles.
- Fixed three audit discoveries: missing Dual entropic OCV tangent, ThinLTO-elided plating overflow failure, and negative spectral lane allocation.
- Mutation checks fail for raw formula copies, wrong flux/Arrhenius wiring, terminal-resistance sign, and modal coefficient duplication.
- Native Debug/Release/CUDA pass 54/54; full WSL Clang 18 ASan+UBSan passes 54/54 with no finding.
- Exact Status coverage: 339 active = 332 measured + 7 structural exceptions; zero uncovered/unmapped.
- Core-only and nested consumers pass 1/1; TSan was not rerun, and M0.2's 3/3 remains applicable because concurrency code did not change.
- Hosted CI was not run; report: `.claude/reports/p9c1-single-source-validation-2026-07-12.md`; next is M0.6/9C-2.

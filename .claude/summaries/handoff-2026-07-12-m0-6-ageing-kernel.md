# Handoff — M0.6 / 9C-2 ageing-kernel idiom (2026-07-12)
- `AgeingKernel.hpp` supplies checked scratch, ordered lane traversal, and atomic dispatch to all four bodies, plus masks/model traversal/clearing where accumulated alternatives need them.
- The pre-refactor 1,077-value all-mask A/B/A trace remains exact in Debug, fast-math Release, and ThinLTO.
- Fixed saturated crack Dual tangents, `pow(Dual,0)`, and invalid signed extents/padding across adjacent SPM scratch/cache/arena construction.
- Six AgeingKernel cases prove representative Dual FD, lane isolation, reset, mappings/order, and failure atomicity; the separate 257-lane P1G2 case proves zero warmed allocations.
- Native Debug/Release/CUDA and WSL Clang 18 ASan+UBSan pass 56/56; core-only and nested consumers pass 1/1.
- Exact Status coverage: 339 active = 332 measured + 7 pinned exceptions; 30 inactive and zero uncovered/unmapped.
- TSan was not rerun and no new TSan claim is made; hosted workflows are committed but not claimed as run.
- Report: `.claude/reports/p9c2-ageing-kernel-validation-2026-07-12.md`; next is M0.7/9C-3 cold-file splitting.

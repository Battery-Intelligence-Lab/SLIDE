# SLIDE Engineering Lessons

This file records reusable findings from the v4 implementation and validation work. Detailed design decisions remain in `PLAN.md`; session chronology remains in `.claude/discussions.md` and `.claude/summaries/`.

## Portability gates must exercise a consumer

- A successful `SLIDE_CORE_ONLY` build followed by a zero-test CTest run is not validation. Keep a dependency-light executable smoke that checks the state arena, owned thread pool, and disabled optional capabilities.
- An in-tree target cannot prove that public CMake usage works. Configure SLIDE as a nested subproject and compile headers through `target_link_libraries(slide_core)`; this exposed the incorrect use of `CMAKE_SOURCE_DIR`, which changes meaning in a superproject.
- “Dependency-free” must be scoped precisely. The v4 core is free of optional CUDA, MATLAB, zstd, Arrow, TBB, and legacy dependencies, but Eigen remains required for cold-path spectral and sparse solves. Prefer an installed package and keep a pinned source fallback.
- CI path filters are part of correctness: changes under `cmake/**` must trigger both core portability and installed-wheel jobs.

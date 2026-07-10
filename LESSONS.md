# SLIDE Engineering Lessons

This file records reusable findings from the v4 implementation and validation work. Detailed design decisions remain in `PLAN.md`; session chronology remains in `.claude/discussions.md` and `.claude/summaries/`.

## Portability gates must exercise a consumer

- A successful `SLIDE_CORE_ONLY` build followed by a zero-test CTest run is not validation. Keep a dependency-light executable smoke that checks the state arena, owned thread pool, and disabled optional capabilities.
- An in-tree target cannot prove that public CMake usage works. Configure SLIDE as a nested subproject and compile headers through `target_link_libraries(slide_core)`; this exposed the incorrect use of `CMAKE_SOURCE_DIR`, which changes meaning in a superproject.
- “Dependency-free” must be scoped precisely. The v4 core is free of optional CUDA, MATLAB, zstd, Arrow, TBB, and legacy dependencies, but Eigen remains required for cold-path spectral and sparse solves. Prefer an installed package and keep a pinned source fallback.
- CI path filters are part of correctness: changes under `cmake/**` must trigger both core portability and installed-wheel jobs.

## Documentation is an executable interface

- Keep runnable examples in one place: mark the fenced Markdown block, extract that exact text, and compile or run it away from the source tree. A copied example file can pass while the published page has already drifted.
- Validate rendered structure as well as tool exit codes. Jekyll returned success for pages without `layout: default` but emitted unthemed HTML; Doxygen returned success while its configured main page was outside `INPUT` and therefore blank. Gates now require the Jekyll front matter and inspect both rendered landing pages.
- A documentation theme is a build dependency. Pin its immutable revision and load every Liquid-tag provider explicitly; the prior remote theme referenced `github_edit_link` without enabling `jekyll-github-metadata`.
- Cross-platform commands must avoid shell-expanded wheel globs. Select the local artifact with pip's `--no-index --find-links=dist` interface instead.
- Do not turn a licensed local MATLAB result into an unverified hosted-runner claim. The same extractor supports the gate, but MEX execution remains on an explicitly provisioned licensed machine.
- Doxygen's 123 current warnings all originate in retained v3 `src/` comments. The v4 docs gate removes generator errors and new-page/configuration warnings without pretending that this legacy debt is already zero.

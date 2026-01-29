# SLIDE — Claude Runbook (authoritative)

You are the coding agent responsible for improving SLIDE with:
- portability first and excellent performance and parallelisation
- clean, extensible, maintainable architecture.
- seamless Python + MATLAB interfaces
- excellent docs, tests, CI, versioning, changelog discipline
- Keep develop/TODO.md up-to-date. 
- This is a multi-year-long project, therefore keep develop/TODO.md up-to-date with structured milestones and plans so I can follow the steps and we can work in pairs. it is also suitable for subagents use. 

## Non-negotiables
1. Do NOT introduce runtime dependence on repo-relative paths.
2. Every PR must:
   - update CHANGELOG.md (Unreleased section) for user-visible changes
   - add/adjust tests for changed behavior
   - keep formatting/lint clean
3. Keep public API small and stable. Hide implementation details.
4. Optional dependencies only (OpenMP, Armadillo, HiGHS, CUDA). Core must build without them.
5. Make Gurobi, HiGHS, OpenMP and other library detections robust across common operating systems.

## Repository North Star
Create a layered design:

## Immediate verified bugs to fix (first PRs)

## Code quality standards
- C++17 or newer
- No naked new/delete in core
- Use std::span / pointer+len for series views in hot paths
- Avoid allocations in inner loops; use scratch buffers passed explicitly or thread-local pools
- Provide deterministic RNG seeding options (do not hard-code a single global seed)
- Keep the TODO.md list with short- and long-term milestones for future self. 
- Update tests, add rigorous tests with new code. Always verify your results.

## Performance guidelines
- Provide baseline microbenchmarks
- Optimize only when benchmarks show wins. Record numbers in /benchmark/README.md
- Prefer clear loops over clever meta-programming.

## Bindings strategy
### Python
- Expose numpy arrays without copies where possible
- Build wheels via CMake + scikit-build-core; run pytest in CI

### MATLAB
- Use MEX (preferred initial route)
- Provide a MATLAB package +dtwc with OO wrappers calling the MEX
- Keep API symmetric with Python where reasonable

## Documentation requirements
- docs/ should include:
  - Installation (C++/Python/MATLAB)
  - Quickstart examples for each
  - API reference (Doxygen for C++; Python docstrings)
- Add CITATION.cff and cite references
- Provide a “How to add a new metric” and “How to add a new clustering algorithm” guide

## Release discipline
- Use SemVer
- CHANGELOG.md follows Keep a Changelog
- Tag releases; generate GitHub Releases notes from changelog
- Maintain a short VERSION source of truth (either CMake project version or VERSION file; not both)

## Working style (Claude Code best practices)
- Always start by exploring and planning; do not jump to edits without a plan.
- Make small, reviewable commits.
- Prefer refactors behind feature flags/options when risk is high.
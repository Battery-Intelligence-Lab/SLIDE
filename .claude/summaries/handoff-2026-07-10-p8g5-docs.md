# Handoff — 2026-07-10 — P8-G5 executable v4 documentation

## Goal

Close the registered P8-G5 documentation gate with truthful v4 guidance and runnable examples, while retaining the v3 material as an explicitly labelled compatibility section.

## What changed

- Added `docs/v4` installation, C++/Python/MATLAB quickstarts, dependency/PyBaMM compatibility, add-cell-model, and add-ageing-mechanism pages; replaced the site/README entry points and bannered the relevant v3 pages.
- Added `scripts/test_docs_quickstarts.py`. It checks local links and fragments, requires themed Jekyll front matter, extracts the exact marked fences, and runs selected C++, Python, or MATLAB examples.
- Rebuilt documentation CI around fresh C++ and installed-wheel Python examples, Doxygen, production Jekyll, and Pages artifacts. The pinned theme now loads its required metadata plugin; Pages/OIDC write permission exists only on master deployment.
- Added a dedicated Doxygen v4 API main page, removed the broken custom header, disabled unneeded LaTeX/Graphviz generation, and kept the API input on C++ sources.
- Documented optional dependencies, ignored PyBaMM compatibility fields, available Solution variables, supported model/discretisation/control/sensitivity/device boundaries, and the absence of a public runtime model-plugin ABI.
- Added generated Jekyll and platform MEX artifacts to `.gitignore`.

## What was tested

- Exact C++ fence as a fresh `add_subdirectory` external consumer: CTest 1/1, 7 samples.
- Exact Python fence against the isolated installed CPython 3.13 wheel: 7 samples, final voltage 3.879196 V.
- Exact MATLAB fence against the licensed R2025b MEX/package: 7 samples, final voltage 3.879196 V.
- Local link/front-matter gate: 20 local links and their Markdown heading fragments validated.
- Doxygen 1.14: exit 0, zero generator errors, non-blank v4 main page, default header with no missing custom assets. The 123 warnings all originate in retained v3 `src/` comments.
- Production GitHub-Pages/Jekyll stack: exit 0, 8 themed v4 HTML pages, `/SLIDE` base URL, and edit links targeting `master/docs`.
- Ruff check/format, yamllint, actionlint 1.7.12, and `git diff --check`.
- Independent read-only claim/CI/security review: no remaining P8-G5 blocker.

## Key results and interpretation

P8-G5 is complete. The Markdown examples are now executable interfaces, not duplicated pseudocode. C++, Python, and MATLAB agree on the bounded quickstart result. Documentation builds are content-checked rather than accepted from exit status alone. MATLAB remains honestly scoped to the licensed local gate; cross-platform hosted workflows are committed definitions and cannot be called green until pushed.

## Recommended next step

Start Phase 9A: close AUD-1, derive and record AUD-2's parity bands, and decide/record the AUD-4 CVODE waiver or short optional arbiter. Keep AUD-3 as Q11 for the user's quiet-machine decision.

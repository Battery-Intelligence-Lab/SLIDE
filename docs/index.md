---
layout: default
title: SLIDE documentation
nav_order: 0
---

# SLIDE

SLIDE is a C++20 simulator for lithium-ion cells and packs. The v4 core combines a spectral single-particle model, compile-time physics composition, contiguous structure-of-arrays state, flat compiled pack topology, Python and MATLAB interfaces, optional CUDA execution, and memory-aware recording.

> v4 is still a development preview. Python reports `4.0.0.dev0` and these docs report `4.0.0-dev`; the root CMake project deliberately remains `3.0.0` until the Phase-11 release gate aligns every public version and creates the tag together.

## Start with v4

- [Install or build the preview](v4/installation.html)
- [C++ quickstart](v4/quickstart-cpp.html)
- [Python quickstart](v4/quickstart-python.html)
- [MATLAB quickstart](v4/quickstart-matlab.html)
- [Optional dependencies and PyBaMM compatibility](v4/compatibility.html)
- [Add a cell-model family](v4/extending-cell-model.html)
- [Add an ageing mechanism](v4/extending-ageing.html)

The three quickstarts are extracted directly from these pages and executed by `scripts/test_docs_quickstarts.py`; they are not illustrative pseudocode.

## Architecture in one paragraph

User-friendly descriptions and language wrappers run on the cold path. At `build()`/`compile()`, they become homogeneous batches over a 64-byte-aligned state arena, concrete physics kernels, and a flat pack netlist. The hot path contains spans, row indices, preallocated workspace, and at most one indirect call per batch—never one virtual dispatch per cell.

## Legacy v3 documentation

The older **Getting Started**, **Usage**, and **Pack Simulation** sections describe the retained v3 façade. They remain useful for existing applications, but they are not the v4 API. v3 stays compiled and runnable through v4.x; its planned removal is v5.

## Citation

If you use SLIDE in research, cite J. M. Reniers, G. Mulder, and D. A. Howey, “Review and performance comparison of mechanical-chemical degradation models for lithium-ion batteries,” *Journal of The Electrochemical Society* 166(14), A3189 (2019), DOI [10.1149/2.0281914jes](https://doi.org/10.1149/2.0281914jes).

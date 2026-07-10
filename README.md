# SLIDE

SLIDE (Simulator for Lithium-Ion Degradation) is a C++20 library for fast lithium-ion cell and pack simulation. The v4 development core provides a spectral SPM, lumped thermal and ageing compositions, flat compiled pack topology, memory-aware recording, Python and MATLAB interfaces, forward sensitivities for PyBOP, and an optional CUDA batch backend.

> v4 is a development preview until the release gate tags `v4.0.0`. The retained v3 façade remains compiled and runnable through v4.x.

## Start here

- [v4 installation and build guide](docs/v4/installation.md)
- [C++ quickstart](docs/v4/quickstart-cpp.md)
- [Python quickstart](docs/v4/quickstart-python.md)
- [MATLAB quickstart](docs/v4/quickstart-matlab.md)
- [Dependencies and PyBaMM compatibility gaps](docs/v4/compatibility.md)
- [Architecture roadmap and validation evidence](PLAN.md)

The documentation quickstarts are executable gates. CI extracts and runs the C++ and Python fences against built artifacts; the same extractor runs the MATLAB fence on the licensed R2025b gate machine. Every case checks finite output from a short solve.

## Minimal builds

Dependency-light C++ core:

```console
cmake -S . -B build-core -DCMAKE_BUILD_TYPE=Release -DSLIDE_CORE_ONLY=ON
cmake --build build-core --config Release --parallel
```

Python wheel from source:

```console
python -m pip install --upgrade build numpy
python -m build --wheel
python -m pip install --force-reinstall --no-index --no-deps --find-links=dist "slide-sim==4.0.0.dev0"
```

Optional CUDA, MATLAB, zstd, and Arrow/Parquet features are off by default. See the [capability matrix](docs/v4/compatibility.md#dependency-and-capability-matrix) before enabling them.

## Model and sign convention

The v4 production model is a Chebyshev-spectral single-particle model with optional lumped thermal physics and SEI, surface-cracking, loss-of-active-material, and lithium-plating mechanisms. Positive current is discharge; negative current is charge. Parameters are SI internally, except capacity in Ah and cumulative energy in Wh where explicitly named.

## Research citation

If you use SLIDE in research, cite:

J. M. Reniers, G. Mulder, and D. A. Howey, “Review and performance comparison of mechanical-chemical degradation models for lithium-ion batteries,” *Journal of The Electrochemical Society* 166(14), A3189 (2019), DOI [10.1149/2.0281914jes](https://doi.org/10.1149/2.0281914jes).

SLIDE is developed at the University of Oxford's Battery Intelligence Lab and is released under the BSD 3-Clause License. See [LICENSE](LICENSE).

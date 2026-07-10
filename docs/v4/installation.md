---
layout: default
title: v4 installation
nav_order: 1
---

# Install or build SLIDE v4

## Requirements

- CMake 3.31 or newer
- a C++20 compiler (Clang, GCC, or MSVC)
- Git for a source checkout
- Eigen 3.4 headers; CMake uses an installed package when available and otherwise fetches the pinned source archive

Python requires Python 3.10 or newer. MATLAB and CUDA are opt-in and are not needed for a CPU build.

Clone the repository:

```console
git clone https://github.com/Battery-Intelligence-Lab/SLIDE.git
cd SLIDE
```

## C++ core

Build only the dependency-light v4 core:

```console
cmake -S . -B build-core -DCMAKE_BUILD_TYPE=Release -DSLIDE_CORE_ONLY=ON
cmake --build build-core --config Release --parallel
```

The preview does not yet install a standalone C++ SDK package; that export is part of the Phase-11 release gate. Until then, consume SLIDE with `add_subdirectory` and link `slide_core`, as shown in the [C++ quickstart](quickstart-cpp.html). This is also how the committed external-consumer gate validates the public headers.

For the retained v3 façade and its tests, omit `SLIDE_CORE_ONLY`:

```console
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

## Python wheel

The distribution name is `slide-sim`; the import name is `slide`. Build from the checkout and install the wheel artifact, not the source tree:

```console
python -m pip install --upgrade build numpy
python -m build --wheel
python -m pip install --force-reinstall --no-index --no-deps --find-links=dist "slide-sim==4.0.0.dev0"
```

The wheel build automatically selects `SLIDE_CORE_ONLY=ON`. NumPy is the only runtime dependency; plotting, SciPy-backed MAT-file export (`slide-sim[matlab]`), PyBaMM comparison, and PyBOP fitting are optional extras.

## MATLAB package

Point CMake at a MATLAB installation if it is not detected automatically:

```console
cmake -S . -B build-matlab -DSLIDE_CORE_ONLY=ON -DSLIDE_WITH_MATLAB=ON -DMatlab_ROOT_DIR="/path/to/MATLAB"
cmake --build build-matlab --config Release --target slide_mex --parallel
```

The build places `slide_mex` beside `matlab/+slide`. Add that `matlab` directory to the MATLAB path. The MEX target uses the interleaved-complex R2018a ABI switch; the object-oriented wrappers use newer MATLAB `arguments` syntax and are validated with R2025b, so R2018a runtime compatibility is not claimed.

## Optional backends

All optional features default off. See the [dependency matrix](compatibility.html#dependency-and-capability-matrix) before enabling CUDA, zstd, Arrow/Parquet, MATLAB, or Python bindings. An explicitly requested but unavailable toolchain fails configuration; an unrequested toolchain never breaks the default CPU configure.

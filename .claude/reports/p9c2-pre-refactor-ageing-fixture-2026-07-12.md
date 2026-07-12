# M0.6 / 9C-2 pre-refactor ageing fixture (2026-07-12)

The production source state is commit `716a1a1`, before any M0.6 ageing
scaffolding change. The fixture enables every SEI (0x0f), surface-crack (0x1f),
and LAM (0x0f) model together with lithium plating, SEI porosity loss, and crack
diffusivity coupling in the real `SpmBatch` pipeline.

Seven lanes vary temperature, signed current, modal concentration, SEI
thickness, crack surface, electrode area/thickness state, stress history, and
history interval. Three evaluations use current patterns A, B, A. The repeated
A derivative must be bit-identical, independently proving that mechanism
scratch is reset. The trace fingerprints the initial heterogeneous arena, every
full derivative arena, and every terminal-voltage vector: 1,077 doubles.

Exact hashes are opt-in through `SLIDE_ENABLE_RECORDED_SCALAR_BITS`. CMake
restricts opt-in to x64 Windows Clang 21.1.8; successful reproduction also
depends on compatible CRT/ISA because Release uses `-march=native`.

| Configuration | Values | FNV-1a | Independent mix |
|---|---:|---:|---:|
| Debug/ThinLTO (`-O0`) | 1,077 | `87119b1b6fa83b81` | `55296fa07292792f` |
| Release/fast-math | 1,077 | `dc8f92a59dd8d67f` | `59085638e06725d3` |
| Release/fast-math/ThinLTO | 1,077 | `0be580cf849e57e1` | `4c6451971b74f789` |

Portable toolchains always enforce the value count, finite outputs, repeated
scratch result, existing per-model legacy parity, and subsequent structural and
Dual finite-difference gates; they do not adopt new post-refactor hashes.

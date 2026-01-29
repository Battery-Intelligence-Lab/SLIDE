# /build-release - Optimized Release Build

Configures and builds SLIDE with full optimizations for production use or benchmarking.

## Usage

```
/build-release [--clean] [--parallel=N]
```

## Steps

### 1. Configure Release Build

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
```

### 2. Build All Targets

```bash
cmake --build build --config Release -j
```

### 3. Verify Build

```bash
# Run quick test to verify functionality
ctest --test-dir build -C Release -j2 --output-on-failure
```

## Optimization Flags by Compiler

### GCC/Clang (Linux/macOS)
```cmake
# Applied automatically in Release mode:
-O3 -DNDEBUG
# Additional flags in CMakeLists.txt:
-march=native -Ofast -ffast-math
```

### MSVC (Windows)
```cmake
# Applied automatically:
/O2 /DNDEBUG
# Additional:
/fp:fast
```

## Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | Release | Build configuration |
| `SLIDE_ENABLE_COVERAGE` | OFF | Code coverage (Debug only) |
| `ENABLE_CPPCHECK` | OFF | Static analysis |
| `ENABLE_CLANG_TIDY` | OFF | Clang static analysis |

## Performance Verification

After building, verify performance:

```bash
# Run benchmark
./build/benchmark/benchmark_Cell_SPM

# Expected times (v3.0.0):
# - 5000 1C CC cycles: < 1 minute
# - CC+CV cycles: < 2 minutes
```

## Clean Build

For a fresh build (resolves some caching issues):

```bash
# Remove build directory
rm -rf build

# Or on Windows:
rmdir /s /q build

# Reconfigure and build
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Multi-Configuration (Visual Studio)

On Windows with Visual Studio generator:

```bash
# Configure (no build type needed)
cmake -B build

# Build specific configuration
cmake --build build --config Release -j

# Run tests for Release
ctest --test-dir build -C Release
```

## Troubleshooting

### Build Fails
1. Check CMake version: `cmake --version` (requires 3.31+)
2. Check compiler: GCC 10+, Clang 12+, or MSVC 2019+
3. Check dependencies fetched correctly (Eigen, Boost, etc.)

### Slow Performance (MSVC)
Known issue: MSVC ~3x slower than Clang due to vectorization differences.
Ensure `/O2` and `/fp:fast` are applied.

### LTO Errors
Link-Time Optimization may fail on some systems:
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF
```

## Output Locations

| Target | Location |
|--------|----------|
| Main executable | `build/slide` or `build/Release/slide.exe` |
| Library | `build/src/libsrc.a` or `build/src/Release/src.lib` |
| Benchmarks | `build/benchmark/` |
| Tests | `build/tests/` |

## Dependencies

Managed via CPM (auto-downloaded):
- Eigen 3.4.0
- Boost 1.89.0 (header-only parts)
- range-v3 0.12.0
- fmt 11.0.2
- NLopt 2.10.0
- Catch2 3.6.0 (tests only)

Third-party licenses generated to: `build/bin/third_party.txt`

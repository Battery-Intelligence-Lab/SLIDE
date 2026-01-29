# SLIDE C++ Style Guide

This document describes the C++ coding conventions used in the SLIDE project, derived from analyzing the existing codebase.

## C++ Standard

- **Required:** C++20
- **Features used:** concepts, ranges, std::span, constexpr improvements
- **CMake minimum:** 3.31

## Naming Conventions

### Classes and Structs
- **PascalCase** (UpperCamelCase)
- Examples: `StorageUnit`, `Cell_SPM`, `Module_p`, `Battery`, `Cycler`, `State_SPM`
- Cell types use underscores: `Cell_SPM`, `Cell_ECM`, `Cell_Bucket`

### Functions
- **camelCase** or **snake_case** (codebase has mixed usage)
- Examples: `setCurrent`, `getStates`, `timeStep_CC`, `validStates`
- Prefer consistency within a file/class

### Member Functions
- **camelCase** preferred
- Examples: `setCurrent()`, `getOCV()`, `timeStep_CC()`, `storeData()`
- Getters: `V()`, `I()`, `T()`, `Cap()` (short names for frequently-called methods)

### Variables
- **snake_case** or **camelCase** (mixed in codebase)
- Examples: `Vcell_valid`, `nch`, `dt`, `tlim`
- State variables often use abbreviated names: `I`, `V`, `T`, `SOC`

### Member Variables
- Public members: descriptive names (`Rcontact`, `Vmin`, `Vmax`)
- Private members: trailing underscore optional
- Constants in classes: UPPER_CASE (`VMIN`, `VMAX`)

### Constants
- **UPPER_SNAKE_CASE** for compile-time constants
- Examples: `MODULE_NSUs_MAX`, `DATASTORE_CELL`, `T_ENV`
- Use `constexpr` over `#define` where possible

### Template Parameters
- Descriptive names or single letters
- Examples: `settings::nch`, `N_RC`, `cell_t`

### Namespaces
- Primary namespace: `slide`
- Nested: `slide::settings`, `slide::util`
- Use `namespace fs = std::filesystem;` for aliases

### Enums
- Enum class names: **PascalCase**
- Enum values: **PascalCase** or descriptive
- Examples: `Status::Success`, `Status::Vmin_violation`

## File Organization

### Header Files
1. Copyright/license block
2. `#pragma once`
3. Project-local includes
4. Standard library includes
5. Third-party includes (Eigen, Boost, etc.)

Example:
```cpp
/*
 * SLIDE - Simulator for Lithium-Ion Degradation
 * Copyright (c) 2024, University of Oxford
 */

#pragma once

#include "State_SPM.hpp"
#include "../StorageUnit.hpp"

#include <span>
#include <vector>

#include <Eigen/Dense>
```

### Source File Organization
- Cell implementations split across multiple files:
  - `Cell_SPM.hpp` - class definition
  - `Cell_SPM.cpp` - constructor, initialization
  - `Cell_SPM_dstate.cpp` - ODE derivatives
  - `Cell_SPM_diffusion.cpp` - diffusion solver
  - `Cell_SPM_degradation.cpp` - aging models
  - `Cell_SPM_thermal.cpp` - thermal calculations

## Formatting

### Indentation
- **2 spaces** (no tabs)
- Configured in `.clang-format`

### Braces
- Opening brace on same line for functions
- Allman style acceptable for class definitions

```cpp
class Cell_SPM : public Cell
{
public:
  double V() override {
    // implementation
  }
};
```

### Line Length
- No hard limit (ColumnLimit: 0 in clang-format)
- Break long parameter lists for readability

### Pointer/Reference Alignment
- Attached to type: `const std::vector<double> &x`

## Modern C++ Features

### Smart Pointers
- `Deep_ptr<T>` for StorageUnit ownership (enables polymorphic deep copy)
- `std::unique_ptr<T>` for exclusive ownership
- Avoid raw `new`/`delete` in core code

```cpp
// Deep_ptr enables polymorphic copies
Deep_ptr<StorageUnit> su = makeBattery(...);
auto copy = su;  // Deep copy via su->copy()
```

### std::span for Array Views
```cpp
// Zero-copy state access
void getStates(std::span<double> s) const;
Status setStates(std::span<double> s);
```

### State Classes
- Inherit from `std::array` for performance
- Provide named accessors for clarity

```cpp
struct State_SPM : public std::array<double, N>
{
  double &I() { return (*this)[0]; }
  double &V() { return (*this)[1]; }
  // ...
};
```

### constexpr
- Use for compile-time configuration:
```cpp
namespace settings {
  constexpr int nch = 5;  // Chebyshev nodes
  constexpr double T_ENV = 298.15;  // Kelvin
}
```

### User-Defined Literals
```cpp
using namespace slide::literals;
double temp = 25.0_degC;  // Converts to Kelvin
```

## Memory Management

### No Raw new/delete
- Use standard containers (`std::vector`, `std::array`)
- Use smart pointers when dynamic allocation needed
- Use `Deep_ptr` for polymorphic battery hierarchy

### Thread-Local Storage
- Use for per-thread scratch buffers in hot paths:
```cpp
thread_local std::vector<double> scratch_buffer;
```

### Avoid Allocations in Hot Paths
- Pre-allocate vectors with known sizes
- Pass output parameters instead of returning vectors
- Use `std::span` for views into existing data

## Error Handling

### Status Codes
- Use `Status` enum for recoverable errors
- Return `Status::Success` on success
- Specific violations: `Status::Vmin_violation`, `Status::Vmax_violation`, etc.

```cpp
Status setCurrent(double I, bool checkV) {
  if (V() < VMIN()) return Status::Vmin_violation;
  // ...
  return Status::Success;
}
```

### Exceptions
- Use for unrecoverable errors or programming errors
- Catch at appropriate boundaries

### Assertions
- Use Catch2 `REQUIRE()` in tests, not raw `assert()`
- `assert()` acceptable for internal invariants in debug builds

## Documentation

### Doxygen Style
```cpp
/**
 * @brief Set the current flowing through the storage unit.
 * @param I Current in Amperes (positive = discharge)
 * @param checkV If true, validate voltage limits
 * @return Status indicating success or violation type
 */
Status setCurrent(double I, bool checkV);
```

### Inline Comments
- Use `//` for single-line comments
- Explain non-obvious logic, not obvious code

## Performance Considerations

### Hot Path Guidelines
- Minimize virtual function calls where possible
- Use `std::span` instead of copying vectors
- Pre-compute values that don't change per time step
- Profile before optimizing

### Benchmarking
- Benchmarks in `benchmark/` folder
- Record results with timestamps in `benchmark/README.md`
- Compare Debug vs Release builds

## Parallelization

### OpenMP
- Must be optional (guarded by `#ifdef _OPENMP`)
- Provide serial fallback
- Use `settings::isParallel` to control at runtime

## Formatting Tool

Run clang-format before committing:
```bash
clang-format -i src/**/*.cpp src/**/*.hpp
```

Or check without modifying:
```bash
clang-format --dry-run --Werror src/**/*.cpp
```

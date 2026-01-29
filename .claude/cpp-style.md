# DTWC++ C++ Style Guide

This document describes the C++ coding conventions used in the DTWC++ project, derived from analyzing the existing codebase.

## C++ Standard

- **Minimum:** C++17
- **Target:** C++17 (avoiding C++20 features for broader compiler compatibility)

## Naming Conventions

### Classes and Structs
- **PascalCase** (UpperCamelCase)
- Examples: `Problem`, `Data`, `Range`, `Index`, `DataLoader`

### Functions
- **snake_case**
- Examples: `readFile`, `load_folder`, `dtwFull`, `dtwBanded`
- Exception: DTW algorithm variants may use mixed case for readability (`dtwFull_L`)

### Member Functions
- **snake_case**
- Examples: `set_numberOfClusters`, `get_name`, `p_vec`, `centroid_of`

### Variables
- **snake_case**
- Examples: `p_vec`, `p_names`, `clusters_ind`, `centroids_ind`

### Member Variables
- **snake_case** with trailing underscore for private members
- Examples: `name_`, `solver_`
- Public members may omit the underscore: `Nc`, `distMat`, `band`

### Constants
- **UPPER_SNAKE_CASE** for compile-time constants
- Examples: `DEFAULT_BAND_LENGTH`, `DEFAULT_MIP_SOLVER`

### Template Parameters
- Descriptive names with `_t` suffix or `T` prefix
- Examples: `data_t`, `Tfun`, `Tpath`

### Namespaces
- **snake_case** or short lowercase names
- Primary namespace: `dtwc`
- Nested: `dtwc::settings`, `dtwc::init`, `dtwc::scores`

### Enums
- Enum class names: **PascalCase**
- Enum values: **PascalCase**
- Examples: `Method::Kmedoids`, `Method::MIP`, `Solver::Gurobi`

## File Organization

### Header Files
1. Doxygen documentation block at top
2. `#pragma once` (not traditional include guards)
3. Project-local includes
4. Standard library includes
5. Third-party includes (e.g., Armadillo)

Example:
```cpp
/**
 * @file Problem.hpp
 * @brief Brief description
 * @date 19 Oct 2022
 * @author Author Name
 */

#pragma once

#include "Data.hpp"
#include "settings.hpp"

#include <string>
#include <vector>

#include <armadillo>
```

### Namespace Aliases
Use `namespace fs = std::filesystem;` at namespace scope.

## Formatting

### Indentation
- **2 spaces** (no tabs)
- Configured in `.clang-format`

### Braces
- **Allman style** for classes, structs, functions
- Opening brace on new line for class/function definitions

```cpp
class Problem
{
public:
  void doSomething();
};

void Problem::doSomething()
{
  if (condition) {
    // short blocks may use attached braces
  }
}
```

### Line Length
- No hard limit enforced (ColumnLimit: 0 in clang-format)
- Prefer readable line breaks for long parameter lists

### Pointer/Reference Alignment
- Attached to type (right alignment): `const std::vector<data_t> &x`

### Template Formatting
- Space after `template` keyword: `template <typename T>`
- No spaces inside angle brackets: `std::vector<int>`

## Documentation

### Doxygen Style
```cpp
/**
 * @brief Short description.
 * @details Longer explanation if needed.
 * @tparam data_t Description of template parameter.
 * @param name Description of parameter.
 * @return Description of return value.
 */
```

### Member Variable Comments
- Use `//!<` or `/*!< */` for inline documentation
```cpp
int Nc{ 1 }; /*!< Number of clusters. */
bool is_distMat_filled{ false }; //!< Whether distance matrix is computed.
```

### Inline Comments
- Use `//` for single-line comments
- Place on same line for short explanations
```cpp
if (&x == &y) return 0; // Same data, distance is 0
```

## Modern C++ Features

### Initialization
- Prefer brace initialization for members:
```cpp
int Nc{ 1 };
bool flag{ false };
```

### Type Inference
- Use `auto` for complex return types and iterators
- Be explicit when type clarity matters

### Move Semantics
- Use `std::move` for efficiency in constructors and functions
```cpp
p_vec = std::move(p_vec_new);
p_vec.push_back(std::move(p));
```

### Structured Bindings (C++17)
```cpp
auto &[short_vec, long_vec] = getVectors();
```

### constexpr
- Use for compile-time constants:
```cpp
constexpr int DEFAULT_BAND_LENGTH = -1;
constexpr data_t maxValue = std::numeric_limits<data_t>::max();
```

### Type Aliases
- Prefer `using` over `typedef`:
```cpp
using distMat_t = arma::Mat<double>;
using path_t = std::decay_t<decltype(settings::resultsPath)>;
```

## Memory Management

### No Raw new/delete
- Use standard containers (`std::vector`, `std::string`)
- Use smart pointers when dynamic allocation is needed

### RAII
- All resources managed through constructors/destructors
- File handles, locks, etc. cleaned up automatically

### Thread-Local Storage
- Use `thread_local` for per-thread scratch buffers:
```cpp
thread_local arma::Mat<data_t> C;
thread_local static std::vector<data_t> buffer(10000);
```

## Error Handling

### Exceptions
- Use `std::runtime_error` for runtime errors
- Provide descriptive error messages
```cpp
throw std::runtime_error("Error opening file: " + path.string());
```

### Validation
- Validate inputs at public API boundaries
- Use assertions for internal invariants in debug builds

## Performance Considerations

### Avoid Allocations in Hot Loops
- Reuse buffers via thread-local storage
- Reserve capacity for vectors when size is known

### Use References for Large Objects
```cpp
void process(const std::vector<data_t> &data);
```

### Prefer operator[] Over .at() in Hot Paths
- `.at()` has bounds checking overhead
- Use `operator[]` when indices are guaranteed valid

## Parallelization

### OpenMP
- Must be optional (guarded by `#ifdef _OPENMP`)
- Provide serial fallback for all parallel code

## Formatting Tool

Run clang-format before committing:
```bash
clang-format -i path/to/file.cpp
```

Or check without modifying:
```bash
clang-format --dry-run --Werror path/to/file.cpp
```

# /add-cell-type - Create New Cell Type

Creates a new cell type inheriting from the Cell base class.

## Usage

```
/add-cell-type NAME [--physics|--empirical]
```

## Required Files

1. **Header**: `src/cells/Cell_NAME/Cell_NAME.hpp`
2. **Implementation**: `src/cells/Cell_NAME/Cell_NAME.cpp` (optional if header-only)
3. **State** (if needed): `src/cells/Cell_NAME/State_NAME.hpp`
4. **Test**: `tests/unit/Cell_NAME_test.cpp`

## Steps

### 1. Create Cell Header

```cpp
// src/cells/Cell_NAME/Cell_NAME.hpp
#pragma once

#include "../Cell.hpp"
#include "State_NAME.hpp"  // if custom state needed

namespace slide {

class Cell_NAME : public Cell
{
public:
  // Constructor
  Cell_NAME();

  // Required virtual methods
  double V() override;                    // Terminal voltage [V]
  double getOCV() override;               // Open circuit voltage [V]
  double getRtot() override;              // Total resistance [Ohm]

  Status setCurrent(double I, bool checkV) override;
  void timeStep_CC(double dt, int nstep) override;

  // State management
  void getStates(std::span<double> s) const override;
  Status setStates(std::span<double> s) override;
  bool validStates() override;
  size_t getNstates() const override;

  // Thermal
  double T() override;
  double getThotSpot() override;

  // Copy
  Cell_NAME* copy() override;

private:
  State_NAME st;  // State container
};

} // namespace slide
```

### 2. Implement Required Methods

Key methods to implement:

- `V()`: Calculate terminal voltage from state
- `getOCV()`: Open circuit voltage (often from lookup table)
- `setCurrent()`: Apply current, check voltage limits
- `timeStep_CC()`: Advance simulation by dt seconds
- `getStates()/setStates()`: Bulk state access via span
- `validStates()`: Check state variables are within physical bounds

### 3. Create State Class (if needed)

```cpp
// src/cells/Cell_NAME/State_NAME.hpp
#pragma once

#include <array>

namespace slide {

struct State_NAME : public std::array<double, N_STATES>
{
  // Named accessors
  double& I() { return (*this)[0]; }
  double& V() { return (*this)[1]; }
  double& T() { return (*this)[2]; }
  double& SOC() { return (*this)[3]; }
  // ... additional states

  const double& I() const { return (*this)[0]; }
  // ... const versions
};

} // namespace slide
```

### 4. Register in cells.hpp

```cpp
// src/cells/cells.hpp
#include "Cell_NAME/Cell_NAME.hpp"
```

### 5. Write Tests

```cpp
// tests/unit/Cell_NAME_test.cpp
#include <catch2/catch_test_macros.hpp>
#include "cells/Cell_NAME/Cell_NAME.hpp"

using namespace slide;

TEST_CASE("Cell_NAME basic operations", "[Cell_NAME]")
{
  Cell_NAME cell;

  SECTION("Initial state is valid") {
    REQUIRE(cell.validStates());
    REQUIRE(cell.V() >= cell.VMIN());
    REQUIRE(cell.V() <= cell.VMAX());
  }

  SECTION("Discharge reduces voltage") {
    double v_initial = cell.V();
    cell.setCurrent(1.0, true);  // 1A discharge
    cell.timeStep_CC(100, 1);    // 100 seconds
    REQUIRE(cell.V() < v_initial);
  }
}
```

### 6. Add to CMakeLists.txt

```cmake
# In src/cells/CMakeLists.txt
target_sources(cells PRIVATE
  Cell_NAME/Cell_NAME.cpp
)
```

## Checklist

- [ ] Header file with class definition
- [ ] All virtual methods implemented
- [ ] State class (if custom states needed)
- [ ] Registered in cells.hpp
- [ ] Unit tests with Catch2
- [ ] Added to CMakeLists.txt
- [ ] Doxygen comments on public methods

## Examples

Existing cell types to reference:
- `Cell_SPM`: Physics-based Single Particle Model
- `Cell_ECM`: Empirical Equivalent Circuit Model (templated)
- `Cell_Bucket`: Simple capacity-tracking cell

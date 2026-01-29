# /add-degradation-model - Add Degradation Mechanism

Adds a new degradation model to Cell_SPM or extends existing mechanisms.

## Usage

```
/add-degradation-model NAME --mechanism=SEI|LAM|CS|LiPlating
```

## Background

SLIDE supports multiple degradation mechanisms:
- **SEI**: Solid Electrolyte Interface growth
- **LAM**: Loss of Active Material
- **CS**: Surface Cracks (crack surface area growth)
- **LiPlating**: Lithium plating on anode

Each mechanism can have multiple model variants selected via `DEG_ID`.

## Steps

### 1. Add Parameters

Create or extend parameter struct in `src/cells/Cell_SPM/param/`:

```cpp
// src/cells/Cell_SPM/param/NEWparam.hpp
#pragma once

namespace slide {

struct NEWparam
{
  double param1{0.0};    //!< Description [units]
  double param2{1.0};    //!< Description [units]
  int model_type{0};     //!< Model selection (0=disabled, 1=model A, etc.)

  // Default constructor
  NEWparam() = default;

  // Parameterized constructor
  NEWparam(double p1, double p2, int type)
    : param1(p1), param2(p2), model_type(type) {}
};

} // namespace slide
```

### 2. Update DEG_ID

Add model flag to `DEG_ID` struct in `src/cells/Cell_SPM/param/DEG_ID.hpp`:

```cpp
struct DEG_ID
{
  // Existing fields...
  int SEI_id{0};
  int SEI_porosity{0};

  // Add new mechanism
  int NEW_id{0};         //!< NEW model: 0=off, 1=model A, 2=model B

  // Update DegArray if needed for iteration
};
```

### 3. Implement Model Equations

Add to appropriate degradation file or create new one:

```cpp
// In Cell_SPM_degradation.cpp or new file

void Cell_SPM::NEW_degradation(double I, double dt)
{
  // Skip if disabled
  if (deg_id.NEW_id == 0) return;

  const auto& p = newparam;  // Parameter reference

  switch (deg_id.NEW_id) {
    case 1: {
      // Model A implementation
      double rate = p.param1 * std::abs(I);
      st.NEW_state() += rate * dt;
      break;
    }
    case 2: {
      // Model B implementation
      // ...
      break;
    }
    default:
      break;
  }
}
```

### 4. Add State Variables (if needed)

Extend `State_SPM` if the model tracks new state:

```cpp
// In State_SPM.hpp
// Add index constant
static constexpr size_t i_NEW_state = 25;  // Example index

// Add accessor
double& NEW_state() { return (*this)[i_NEW_state]; }
const double& NEW_state() const { return (*this)[i_NEW_state]; }
```

### 5. Call from dstate

Integrate into `Cell_SPM_dstate.cpp`:

```cpp
void Cell_SPM::dstate(double dt)
{
  // ... existing code ...

  // Call new degradation model
  NEW_degradation(I(), dt);
}
```

### 6. Add Tests

```cpp
// tests/unit/Cell_SPM_degradation_test.cpp

TEST_CASE("NEW degradation model", "[Cell_SPM][degradation]")
{
  Cell_SPM cell;

  // Enable the model
  cell.deg_id.NEW_id = 1;
  cell.newparam = NEWparam(0.001, 1.0, 1);

  // Run cycles
  double initial_state = cell.st.NEW_state();
  cell.setCurrent(1.0, false);
  cell.timeStep_CC(3600, 360);  // 1 hour

  // Verify degradation occurred
  REQUIRE(cell.st.NEW_state() > initial_state);
}
```

### 7. Document in DEG_ID

Update documentation for model selection:

```cpp
/**
 * @brief NEW_id - NEW degradation model selection
 *
 * Values:
 * - 0: Disabled (default)
 * - 1: Model A - description [Reference]
 * - 2: Model B - description [Reference]
 */
int NEW_id{0};
```

## Checklist

- [ ] Parameter struct created/extended
- [ ] DEG_ID updated with model flag
- [ ] Model equations implemented
- [ ] State variables added (if needed)
- [ ] Called from dstate() or timeStep_CC()
- [ ] Unit tests written
- [ ] Doxygen documentation added
- [ ] Reference paper cited in comments

## Existing Models Reference

| Mechanism | File | Models |
|-----------|------|--------|
| SEI | Cell_SPM_degradation.cpp | SEI_id: 1-4 |
| LAM | Cell_SPM_degradation.cpp | LAM_p, LAM_n models |
| CS | Cell_SPM_degradation.cpp | CS models |
| LiPlating | Cell_SPM_degradation.cpp | Plating models |

See existing implementations for patterns and conventions.

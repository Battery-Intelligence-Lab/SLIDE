# SLIDE Python Style Guide

This document describes the Python coding conventions for the SLIDE project bindings and scripts.

> **Note:** Python bindings are planned but not yet implemented. This guide establishes conventions for when they are built.

## Python Version

- **Minimum:** Python 3.9
- **Recommended:** Python 3.11+

## Style Standard

Follow [PEP 8](https://peps.python.org/pep-0008/) with the following project-specific guidelines.

## Naming Conventions

### Modules
- **snake_case**, lowercase
- Examples: `slide_core.py`, `battery.py`, `cell_spm.py`

### Classes
- **PascalCase**
- Mirror C++ class names where appropriate
- Examples: `StorageUnit`, `CellSPM`, `CellECM`, `ModuleS`, `ModuleP`, `Battery`, `Cycler`

### Functions and Methods
- **snake_case**
- Examples: `set_current`, `get_states`, `time_step_cc`, `run_cccv`

### Variables
- **snake_case**
- Examples: `cell_voltage`, `time_step`, `n_cycles`

### Constants
- **UPPER_SNAKE_CASE**
- Examples: `DEFAULT_TEMPERATURE`, `VMIN`, `VMAX`

### Private/Internal
- Single leading underscore: `_internal_method`

## Type Hints

Use type hints for function signatures (Python 3.9+ style):

```python
import numpy as np
from numpy.typing import NDArray

def set_current(self, current: float, check_voltage: bool = True) -> bool:
    ...

def get_states(self) -> NDArray[np.float64]:
    ...
```

For complex types:
```python
from typing import Optional, Union

def run_cycle(
    self,
    current: float,
    v_limit: float,
    t_limit: float,
    dt: float = 1.0,
    n_data_points: Optional[int] = None
) -> CycleResult:
    ...
```

## Imports

### Order
1. Standard library imports
2. Related third-party imports
3. Local application imports

Separate each group with a blank line:

```python
import os
from pathlib import Path

import numpy as np
import pandas as pd

from slide import Battery, CellSPM, Cycler
```

## Docstrings

Use NumPy-style docstrings for consistency with scientific Python packages:

```python
def run_cc(
    self,
    current: float,
    v_limit: float,
    t_limit: float,
    dt: float = 1.0
) -> CycleResult:
    """
    Run constant current (CC) phase on the battery.

    Parameters
    ----------
    current : float
        Applied current in Amperes. Positive = discharge.
    v_limit : float
        Voltage limit in Volts (min for discharge, max for charge).
    t_limit : float
        Time limit in seconds.
    dt : float, optional
        Time step in seconds, by default 1.0.

    Returns
    -------
    CycleResult
        Object containing time, voltage, current, and temperature arrays.

    Examples
    --------
    >>> cell = CellSPM()
    >>> cycler = Cycler(cell)
    >>> result = cycler.run_cc(current=1.0, v_limit=2.7, t_limit=3600)
    >>> print(f"Final voltage: {result.voltage[-1]:.3f} V")

    See Also
    --------
    run_cv : Run constant voltage phase.
    run_cccv : Run combined CC-CV cycle.

    Notes
    -----
    The sign convention follows the C++ library:
    - Positive current = discharge (current flows out)
    - Negative current = charge (current flows in)
    """
    ...
```

## Class Design

### PyBaMM-Compatible API

Target compatibility with PyBaMM Experiment interface:

```python
class Battery:
    """Lithium-ion battery simulation using Single Particle Model.

    Parameters
    ----------
    capacity : float
        Nominal capacity in Ah.
    v_min : float
        Minimum voltage limit in V.
    v_max : float
        Maximum voltage limit in V.
    temperature : float, optional
        Initial temperature in Kelvin, by default 298.15 K.

    Attributes
    ----------
    voltage : float
        Current terminal voltage in V.
    current : float
        Current flowing through battery in A.
    temperature : float
        Cell temperature in K.
    soc : float
        State of charge (0.0 to 1.0).
    """

    def __init__(
        self,
        capacity: float,
        v_min: float = 2.7,
        v_max: float = 4.2,
        temperature: float = 298.15
    ):
        self._handle = _slide_core.create_battery(...)

    @property
    def voltage(self) -> float:
        """Current terminal voltage in Volts."""
        return self._handle.V()

    @property
    def current(self) -> float:
        """Current in Amperes (positive = discharge)."""
        return self._handle.I()

    def set_current(self, current: float) -> bool:
        """Apply a current to the battery.

        Parameters
        ----------
        current : float
            Current in Amperes (positive = discharge).

        Returns
        -------
        bool
            True if successful, False if voltage limits violated.
        """
        return self._handle.setCurrent(current)

    def step(self, dt: float) -> None:
        """Advance simulation by one time step.

        Parameters
        ----------
        dt : float
            Time step in seconds.
        """
        self._handle.timeStep_CC(dt)
```

## NumPy Integration

### Array Handling
- Accept both lists and numpy arrays
- Convert to numpy internally if needed
- Return numpy arrays for consistency
- Use zero-copy where possible via pybind11/nanobind

```python
def get_states(self) -> NDArray[np.float64]:
    """Get all internal state variables.

    Returns
    -------
    NDArray[np.float64]
        Array of state variables (see State_SPM for order).
    """
    # Zero-copy view into C++ state array
    return np.asarray(self._handle.getStatesView())
```

### Memory Efficiency
- Use `np.ascontiguousarray()` before passing to C++ if needed
- Document when zero-copy is possible
- Avoid unnecessary copies in hot paths

## Error Handling

### Exceptions
- Use built-in exceptions when appropriate
- Create custom exceptions for domain-specific errors

```python
class SlideError(Exception):
    """Base exception for SLIDE-related errors."""
    pass

class VoltageLimitError(SlideError):
    """Raised when voltage limits are violated."""
    pass

class TemperatureLimitError(SlideError):
    """Raised when temperature limits are violated."""
    pass
```

### Validation
- Validate inputs at public API boundaries
- Provide clear error messages

```python
def set_current(self, current: float) -> bool:
    if not isinstance(current, (int, float)):
        raise TypeError(f"current must be numeric, got {type(current)}")
    if abs(current) > self.max_current:
        raise ValueError(
            f"current magnitude {abs(current):.2f} A exceeds "
            f"maximum {self.max_current:.2f} A"
        )
```

## Testing

### pytest Style
```python
import pytest
import numpy as np
from slide import CellSPM, Cycler


class TestCellSPM:
    @pytest.fixture
    def cell(self):
        """Create a default SPM cell for testing."""
        return CellSPM()

    def test_initial_voltage(self, cell):
        """Cell should have valid initial voltage."""
        assert 2.7 <= cell.voltage <= 4.2

    def test_discharge_reduces_voltage(self, cell):
        """Discharging should reduce voltage."""
        initial_v = cell.voltage
        cycler = Cycler(cell)
        cycler.run_cc(current=1.0, v_limit=2.7, t_limit=100)
        assert cell.voltage < initial_v

    @pytest.mark.parametrize("current", [0.5, 1.0, 2.0])
    def test_discharge_currents(self, cell, current):
        """Test various discharge currents."""
        cycler = Cycler(cell)
        result = cycler.run_cc(current=current, v_limit=2.7, t_limit=60)
        assert len(result.time) > 0
```

### Fixtures
```python
@pytest.fixture
def sample_battery():
    """Create a sample battery configuration."""
    return Battery(capacity=3.0, v_min=2.7, v_max=4.2)

@pytest.fixture
def degradation_params():
    """Default degradation parameters for testing."""
    return {
        "sei_enabled": True,
        "lam_enabled": False,
        "crack_enabled": False
    }
```

## Formatting Tools

### Black
Use Black for automatic formatting:
```bash
black python/
```

### isort
Use isort for import sorting (compatible with Black):
```bash
isort python/
```

### Configuration
In `pyproject.toml`:
```toml
[tool.black]
line-length = 88
target-version = ['py39', 'py310', 'py311', 'py312']

[tool.isort]
profile = "black"
```

## Project Structure

```
python/
├── slide/
│   ├── __init__.py          # Package init, version, public API
│   ├── _core.pyi            # Type stubs for C++ bindings
│   ├── battery.py           # Battery class wrapper
│   ├── cell.py              # Cell classes (SPM, ECM)
│   ├── module.py            # Module classes (series, parallel)
│   ├── cycler.py            # Cycler and procedures
│   └── plotting.py          # Visualization utilities
├── tests/
│   ├── __init__.py
│   ├── conftest.py          # pytest fixtures
│   ├── test_cell.py
│   ├── test_cycler.py
│   └── test_battery.py
├── examples/
│   ├── basic_cycling.py
│   ├── degradation_study.py
│   └── pybamm_comparison.py
└── src/
    └── py_main.cpp          # pybind11/nanobind bindings
```

## Package Metadata

In `__init__.py`:
```python
"""SLIDE - Simulator for Lithium-Ion Degradation.

A fast lithium-ion battery simulation library with degradation modeling.
"""

from ._version import __version__
from .cell import CellSPM, CellECM
from .module import ModuleS, ModuleP
from .battery import Battery
from .cycler import Cycler

__all__ = [
    "CellSPM",
    "CellECM",
    "ModuleS",
    "ModuleP",
    "Battery",
    "Cycler",
    "__version__",
]
```

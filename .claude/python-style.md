# DTWC++ Python Style Guide

This document describes the Python coding conventions for the DTWC++ project bindings and scripts.

## Python Version

- **Minimum:** Python 3.7 (for broad compatibility)
- **Recommended:** Python 3.9+

## Style Standard

Follow [PEP 8](https://peps.python.org/pep-0008/) with the following project-specific guidelines.

## Naming Conventions

### Modules
- **snake_case**, lowercase
- Examples: `setup.py`, `py_main.cpp`, `ucr_benchmark.py`

### Classes
- **PascalCase**
- Examples: `CMakeExtension`, `CMakeBuild`, `DTW`, `KMedoids`

### Functions and Methods
- **snake_case**
- Examples: `build_extension`, `get_ext_fullpath`, `fit`, `predict`

### Variables
- **snake_case**
- Examples: `ext_fullpath`, `cmake_args`, `n_clusters`

### Constants
- **UPPER_SNAKE_CASE**
- Examples: `PLAT_TO_CMAKE`, `DEFAULT_BAND_WIDTH`

### Private/Internal
- Single leading underscore: `_internal_method`
- Double underscore for name mangling (rare): `__private`

## Type Hints

Use type hints for function signatures (Python 3.7+ style):

```python
def build_extension(self, ext: CMakeExtension) -> None:
    ...

def fit(self, X: np.ndarray) -> "KMedoids":
    ...
```

For complex types, use `typing` module:
```python
from typing import List, Optional, Union

def cluster(
    self,
    data: List[np.ndarray],
    n_clusters: int,
    initial_medoids: Optional[List[int]] = None
) -> np.ndarray:
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
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from dtwcpp import DTW, KMedoids
```

### Import Style
- Prefer explicit imports over `from module import *`
- Group related imports on same line if short

## Docstrings

Use NumPy-style docstrings for consistency with scientific Python packages:

```python
def dtw_distance(x: np.ndarray, y: np.ndarray, band: int = -1) -> float:
    """
    Compute Dynamic Time Warping distance between two time series.

    Parameters
    ----------
    x : np.ndarray
        First time series (1D array).
    y : np.ndarray
        Second time series (1D array).
    band : int, optional
        Sakoe-Chiba band width. -1 means no constraint (default).

    Returns
    -------
    float
        The DTW distance between x and y.

    Examples
    --------
    >>> x = np.array([1, 2, 3, 4])
    >>> y = np.array([1, 2, 2, 3, 4])
    >>> dtw_distance(x, y)
    1.0

    See Also
    --------
    distance_matrix : Compute pairwise DTW distances.
    """
    ...
```

## Class Design

### Scikit-learn Compatible API

For clustering classes, follow scikit-learn conventions:

```python
class KMedoids:
    """K-Medoids clustering using DTW distance.

    Parameters
    ----------
    n_clusters : int
        Number of clusters.
    metric : str, default="dtw"
        Distance metric to use.
    random_state : int, optional
        Random seed for reproducibility.

    Attributes
    ----------
    labels_ : np.ndarray
        Cluster labels for each sample.
    cluster_centers_ : np.ndarray
        Indices of medoid samples.
    inertia_ : float
        Sum of distances to closest medoid.
    """

    def __init__(
        self,
        n_clusters: int,
        metric: str = "dtw",
        random_state: Optional[int] = None
    ):
        self.n_clusters = n_clusters
        self.metric = metric
        self.random_state = random_state

    def fit(self, X: np.ndarray) -> "KMedoids":
        """Fit the K-Medoids model.

        Parameters
        ----------
        X : np.ndarray of shape (n_samples, n_features)
            Training data.

        Returns
        -------
        self
            Fitted estimator.
        """
        ...
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Predict cluster labels for samples."""
        ...

    def fit_predict(self, X: np.ndarray) -> np.ndarray:
        """Fit and return cluster labels."""
        return self.fit(X).labels_
```

## NumPy Integration

### Array Handling
- Accept both lists and numpy arrays
- Convert to numpy internally if needed
- Return numpy arrays for consistency

```python
def distance_matrix(X):
    X = np.asarray(X)  # Convert if needed
    # ... compute ...
    return result  # Return numpy array
```

### Memory Efficiency
- Use `np.ascontiguousarray()` before passing to C++ if needed
- Document when zero-copy is possible

## Error Handling

### Exceptions
- Use built-in exceptions when appropriate
- Create custom exceptions for domain-specific errors

```python
class DTWError(Exception):
    """Base exception for DTW-related errors."""
    pass

class InvalidBandWidthError(DTWError):
    """Raised when band width is invalid."""
    pass
```

### Validation
- Validate inputs at public API boundaries
- Provide clear error messages

```python
def fit(self, X):
    if X.ndim != 2:
        raise ValueError(f"X must be 2D, got {X.ndim}D")
    if X.shape[0] < self.n_clusters:
        raise ValueError(
            f"n_samples={X.shape[0]} must be >= n_clusters={self.n_clusters}"
        )
```

## Testing

### pytest Style
```python
import pytest
import numpy as np
from dtwcpp import DTW, KMedoids


class TestDTW:
    def test_same_series_zero_distance(self):
        """Distance to self should be zero."""
        x = np.array([1, 2, 3, 4, 5])
        dtw = DTW()
        assert dtw.distance(x, x) == 0.0

    def test_symmetric(self):
        """DTW distance should be symmetric."""
        x = np.array([1, 2, 3])
        y = np.array([1, 2, 2, 3])
        dtw = DTW()
        assert dtw.distance(x, y) == dtw.distance(y, x)

    @pytest.mark.parametrize("band", [1, 5, 10])
    def test_banded_dtw(self, band):
        """Banded DTW should not exceed full DTW."""
        x = np.random.randn(100)
        y = np.random.randn(100)
        dtw_full = DTW().distance(x, y)
        dtw_banded = DTW(band=band).distance(x, y)
        assert dtw_banded >= dtw_full
```

### Fixtures
```python
@pytest.fixture
def sample_data():
    """Generate sample time series data."""
    np.random.seed(42)
    return np.random.randn(50, 100)  # 50 series, length 100
```

## Formatting Tools

### Black
Use Black for automatic formatting:
```bash
black path/to/file.py
```

### isort
Use isort for import sorting (compatible with Black):
```bash
isort path/to/file.py
```

### Configuration
In `pyproject.toml`:
```toml
[tool.black]
line-length = 88
target-version = ['py37', 'py38', 'py39', 'py310', 'py311']

[tool.isort]
profile = "black"
```

## Project Structure

```
python/
├── dtwcpp/
│   ├── __init__.py      # Package init, version, public API
│   ├── _core.py         # Internal wrapper around C++ bindings
│   ├── dtw.py           # DTW class
│   ├── clustering.py    # KMedoids, CLARA classes
│   └── metrics.py       # Evaluation metrics
├── tests/
│   ├── __init__.py
│   ├── test_dtw.py
│   ├── test_clustering.py
│   └── conftest.py      # pytest fixtures
└── py_main.cpp          # pybind11 bindings
```

## Package Metadata

In `__init__.py`:
```python
"""DTWC++ - Dynamic Time Warping Clustering Library."""

from ._version import __version__
from .dtw import DTW
from .clustering import KMedoids, CLARA

__all__ = ["DTW", "KMedoids", "CLARA", "__version__"]
```

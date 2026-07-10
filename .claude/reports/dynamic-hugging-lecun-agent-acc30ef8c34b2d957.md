# SLIDE Research Report: Literature, Usage, Strengths, Limitations, and Comparisons

> **Note**: WebSearch and WebFetch tools were both unavailable during this research session. This report combines information extracted from the local repository (README, JOSS paper draft, CHANGELOG, docs, TODO, etc.) with knowledge from training data (cutoff: May 2025). URLs are provided where known, but could not be verified live.

---

## 1. Original SLIDE Papers

### Primary Paper (Cite This)
- **Title**: "Review and performance comparison of mechanical-chemical degradation models for lithium-ion batteries"
- **Authors**: Jorn M. Reniers, Grietus Mulder, David A. Howey
- **Journal**: Journal of The Electrochemical Society, 166(14), A3189
- **Year**: 2019
- **DOI**: [10.1149/2.0281914jes](https://doi.org/10.1149/2.0281914jes)
- **Key contribution**: Comprehensive review and benchmark of multiple degradation models (SEI growth, LAM, surface cracking, Li-plating) within a unified SPM framework. SLIDE was released as the companion open-source code.

### Digital Twin Paper
- **Title**: "Digital twin of a MWh-scale grid battery system for efficiency and degradation analysis"
- **Authors**: Jorn M. Reniers, David A. Howey
- **Year**: 2022
- **Status**: Referenced in JOSS paper bibliography; demonstrates SLIDE applied to real grid-scale battery systems.

### JOSS Paper (Draft/In Progress)
- **Title**: "SLIDE: A C++ software for high-fidelity simulation of lithium-ion battery energy storage system degradation"
- **Authors**: Jorn M. Reniers, Volkan Kumtepeli, David A. Howey
- **Affiliation**: Department of Engineering Science, University of Oxford; Brill Power (Reniers)
- **Date**: 18 May 2022 (draft)
- **Status**: The paper.md contains "Disclaimer: paper writing is still ongoing; please do not use this version as a reference." -- likely not yet published.
- **Location in repo**: `joss/paper.md`

### Foundational SPM Implementation
- **Spectral Single Particle Model** by Adrien Bizeray, Jorn Reniers, and David Howey
- **GitHub**: https://github.com/davidhowey/Spectral_li-ion_SPM
- SLIDE extends this with degradation mechanisms and an eigenspace transformation for computational efficiency.

### Zenodo DOI
- Latest release DOI: https://zenodo.org/badge/latestdoi/185216614

---

## 2. How SLIDE Has Been Used in Research

Based on training knowledge (papers citing SLIDE up to May 2025):

### Known Uses
1. **Grid-scale battery digital twins** -- Reniers & Howey used SLIDE to build a digital twin of a MWh-scale grid battery system, modeling degradation under real market dispatch profiles.

2. **Degradation model comparison** -- The original paper itself is a benchmark study, comparing mechanical-chemical degradation models. Researchers have used SLIDE to reproduce and extend these comparisons.

3. **Optimal battery operation** -- Studies on optimal charging protocols and cycle-life extension have used SLIDE to simulate thousands of cycles rapidly, testing different C-rates, voltage windows, and rest periods.

4. **Pack-level simulation** -- After the SLIDE-pack merger (v3.0), the tool supports multi-cell battery pack degradation simulation with thermal coupling, used for studying cell-to-cell variability effects.

5. **Energy storage system design** -- The JOSS paper references use in renewable energy grid balancing, frequency control, and wholesale energy trading simulation contexts.

### Related/Competing Software Referenced in SLIDE's JOSS Paper
- Kumtepeli et al. (2020) - Energy arbitrage optimization with battery storage
- SimSES (Naumann et al., 2017; Moller et al., 2022) - Stationary energy storage simulation
- liionpack (Tranter et al., 2022) - Python pack simulation built on PyBaMM

---

## 3. SLIDE's Documented Strengths

### Speed (Primary Selling Point)
- **5000 1C CC cycles in under 1 minute** (CC only)
- **Under 2 minutes with CV phase** added
- 100 1C CC cycles in ~0.9 seconds
- Written in C++ for maximum computational performance
- Eigenspace transformation of the spectral SPM gives diagonal state-space matrix, improving both speed and numerical accuracy

### Comprehensive Degradation Modeling
- Multiple degradation mechanisms in a single coupled model:
  - SEI growth (multiple variants)
  - Loss of Active Material (LAM)
  - Surface cracking
  - Lithium plating
- Users can toggle individual models on/off via `DEG_ID`
- Fitting parameters are user-configurable

### Battery Tester Paradigm
- Programmed like a physical battery tester: CC, CV, CCCV, current profiles
- Pre-built procedures for calendar aging, cycle aging, drive cycle aging
- Reference performance tests: capacity measurement, OCV curves, pulse discharge
- Data stored at configurable time intervals, mimicking real tester output

### Pack-Level Simulation (v3.0+)
- Hierarchical architecture: Cell -> Module (series/parallel) -> Battery
- Polymorphic design: external code works identically with a single cell or 500 kWh battery
- Thermal coupling between cells, modules with cooling systems (HVAC, open-loop)
- Cell-to-cell variability modeling

### Physics-Based Model
- Single Particle Model with coupled bulk thermal model
- Spectral Chebyshev discretization (high accuracy with few nodes)
- Entropic heat generation
- Temperature-dependent parameters

### Modern C++ (v3.0)
- C++20 standard with concepts, ranges, std::span
- Deep_ptr for polymorphic deep copy
- Status codes for error handling
- Multi-threaded parallelization for independent degradation simulations

---

## 4. SLIDE's Documented Limitations

### From the TODO/Known Issues (develop/TODO.md)
1. **`redistributeCurrent_new` requires up to 2500 iterations** -- major performance bottleneck for parallel module current distribution
2. **MSVC ~3x slower than Clang** -- vectorization differences
3. **Thermal model bugs**: `therm.Qcontact` 6x overestimated; T_MODEL==2 causes thermal runaway
4. **CV implementation issues** -- "CV is not doing intended thing for both SLIDE and slide-pack"
5. **Series module voltage limit handling unclear** -- what to do if one cell reaches max?
6. **No Python bindings** -- planned but not yet implemented
7. **No MATLAB MEX interface** -- MATLAB scripts are for post-processing only (reading CSV)
8. **determineOCV is inefficient** -- search space can be reduced
9. **Cell_SPM holds all aging model parameters** -- should use composition pattern
10. **Global settings** (`settings::isParallel`, `T_MODEL`, `T_ENV` as constants) limit flexibility
11. **Data I/O bottleneck** -- writing CSV files during long degradation simulations can take hours, even though calculations take minutes
12. **Large memory for degradation states** -- `double degState[CELL_NSTATE_MAX]` causes 700 kB unnecessary memory

### Inherent Model Limitations
- **Single Particle Model only** -- no full Pseudo-2D (P2D/DFN) model, limiting accuracy for high C-rates or thick electrodes
- **Lumped thermal model** -- no spatially-resolved thermal simulation within cells
- **No electrolyte dynamics** -- SPM assumes uniform electrolyte concentration
- **No multi-physics coupling** beyond electrochemistry + thermal + degradation (no mechanical stress at continuum level)

### Usability Limitations
- Requires C++ knowledge to extend or modify
- Results output as CSV files, requiring separate post-processing
- Documentation partially outdated (references TDM-GCC, Eclipse IDE)
- No GUI or interactive interface
- Chebyshev discretization was historically MATLAB-dependent (now moved to C++ in v3)

---

## 5. Feature Requests (from TODO.md and Community)

### From develop/TODO.md
- **Python bindings** via pybind11/nanobind (high demand)
- **MATLAB MEX interface**
- **PyBaMM-compatible Experiment interface** (drop-in replacement)
- **Julia wrapper** (suggestion from "Brady")
- **SUNDIALS solver option** (suggestion from Martin Robinson, Oxford RSE)
- **GPU acceleration** (optional)
- **JSON configuration files** for modules (nested structure)
- **GITT function** (Galvanostatic Intermittent Titration Technique)
- **SOC/Temperature dependent RC pairs** for Cell_ECM
- **Snapshot testing framework**
- **Automated tests against PyBaMM**
- **Matio and parquet data format support**
- **Higher-order integration** (RK4) and adaptive time stepping
- **Header-only library** portions for easy compilation
- **Visitor pattern** for hierarchical structures
- **Factory methods** for battery creation
- **Config file model** (suggestion from Martin Robinson)

### From Code Review (Martin Robinson, Oxford RSE)
- Config file + model file + options separation
- Standardized class naming (CamelCase)
- Template project structure following OxRSE conventions

---

## 6. Comparison with PyBaMM

### Architecture Differences
| Aspect | SLIDE | PyBaMM |
|--------|-------|--------|
| Language | C++ (C++20) | Python (with CasADi/SUNDIALS backends) |
| Primary Model | SPM only | SPM, SPMe, DFN/P2D, and more |
| Speed | Very fast (5000 cycles < 1 min) | Slower for equivalent SPM, but flexible |
| Degradation | Built-in coupled degradation models | Degradation via submodels (growing ecosystem) |
| Pack simulation | Built-in (v3.0+) | Via liionpack (separate package) |
| Extensibility | Requires C++ knowledge | Python-based, easier to extend |
| Solver | Custom forward Euler / RK4 | SUNDIALS (IDA/CVODES), CasADi |
| Symbolic math | None | Full symbolic model definition |
| Parameter sets | Hardcoded C++ structs | BPX format, parameter sets library |
| Community | Small (Oxford Battery Intelligence Lab) | Large (PyBaMM team + global community) |
| License | BSD 3-clause | BSD 3-clause |

### SLIDE's Advantages Over PyBaMM
1. **Raw computational speed** -- C++ with eigenspace-transformed SPM is significantly faster for large-scale degradation sweeps
2. **Integrated pack-level degradation** -- single codebase handles cell-to-cell variability + thermal coupling + degradation
3. **Battery tester paradigm** -- more intuitive for experimentalists who think in terms of CC/CV/CCCV procedures
4. **Lower overhead** -- no Python interpreter, no symbolic compilation step

### PyBaMM's Advantages Over SLIDE
1. **Model diversity** -- SPMe, DFN/P2D, lead-acid, and custom models
2. **Symbolic model definition** -- users can modify equations without recompiling
3. **Python ecosystem** -- Jupyter notebooks, matplotlib, scipy integration
4. **Active community** -- larger development team, more frequent releases, extensive documentation
5. **Parameter management** -- BPX standard, Chen 2020, Marquis 2019 parameter sets
6. **Advanced solvers** -- adaptive time stepping via SUNDIALS, DAE support
7. **Experiment API** -- human-readable experiment definitions

### From SLIDE's TODO.md
The SLIDE team explicitly plans to create a "PyBaMM-compatible Experiment interface (drop-in replacement)" and "Automated tests against PyBaMM" -- acknowledging PyBaMM as the de facto standard interface.

---

## 7. SLIDE's Unique Selling Points

1. **Fastest open-source degradation simulator** -- C++ implementation with eigenspace-optimized SPM. No other open-source tool can simulate 5000 degradation cycles in under 1 minute.

2. **Unified pack + degradation simulation** -- SLIDE (post v3.0 merge with slide-pack) is unique in coupling physics-based cell degradation with pack-level electrical and thermal simulation in a single C++ codebase. PyBaMM requires a separate liionpack package.

3. **Multiple coupled degradation mechanisms** -- The original 2019 paper's contribution was the first comprehensive open-source implementation comparing and coupling SEI, LAM, surface cracking, and Li-plating models in one simulator.

4. **Hierarchical polymorphic architecture** -- The StorageUnit -> Cell -> Module -> Battery hierarchy with Deep_ptr enables arbitrary nesting of series/parallel configurations with automatic electrical and thermal consistency.

5. **Digital twin capability** -- Designed and demonstrated for real-time or faster-than-real-time simulation of MWh-scale systems, enabling digital twin applications for grid batteries.

6. **Thermal management system simulation** -- Built-in CoolSystem hierarchy (HVAC, open-loop, passive) with control strategies, unique among open-source battery degradation tools.

---

## 8. Oxford Battery Modelling Group (Battery Intelligence Lab)

### Overview
- **Lab**: Battery Intelligence Lab (formerly Oxford Battery Modelling Group)
- **PI**: Prof. David A. Howey, Department of Engineering Science, University of Oxford
- **Website**: https://howey.eng.ox.ac.uk

### Key People (Related to SLIDE)
- **Jorn M. Reniers** -- Original creator of SLIDE and SLIDE-PACK. Now at Brill Power, Oxford. Wrote several papers using SLIDE including the digital twin paper.
- **Volkan Kumtepeli** -- Current maintainer. Merged SLIDE + SLIDE-PACK, modernized to C++20, wrote JOSS paper.
- **David A. Howey** -- Project administrator, supervision, co-author on all papers.
- **Adrien Bizeray** -- Developed the foundational Spectral SPM that SLIDE builds upon.

### Related Tools from the Group
1. **Spectral_li-ion_SPM** (https://github.com/davidhowey/Spectral_li-ion_SPM) -- MATLAB implementation of spectral SPM, foundation for SLIDE
2. **Brill Power** (https://www.brillpower.com) -- Spin-out company from the lab, Jorn Reniers now works there. Active battery management and intelligence products.

### Research Focus Areas
- Lithium-ion battery degradation modeling and prediction
- Battery state estimation (SOC, SOH)
- Battery management systems
- Grid-scale energy storage optimization
- Digital twins for battery systems
- Fast charging optimization
- Cell-to-cell variability

### External Collaborators Referenced in SLIDE
- **Martin Robinson** (Oxford RSE) -- Code review suggestions for SLIDE
- **Brady** -- Suggested Julia wrapper

---

## Key URLs

| Resource | URL |
|----------|-----|
| SLIDE GitHub | https://github.com/Battery-Intelligence-Lab/SLIDE |
| Battery Intelligence Lab | https://howey.eng.ox.ac.uk |
| Primary Paper DOI | https://doi.org/10.1149/2.0281914jes |
| Zenodo DOI | https://zenodo.org/badge/latestdoi/185216614 |
| Spectral SPM | https://github.com/davidhowey/Spectral_li-ion_SPM |
| SLIDE Documentation | https://Battery-Intelligence-Lab.github.io/SLIDE/ |
| Volkan Kumtepeli GitHub | https://github.com/ElektrikAkar |
| David Howey GitHub | https://github.com/davidhowey |

---

## Limitations of This Report

**WebSearch and WebFetch tools were both unavailable** during this session, so I could not:
- Verify current GitHub star/fork counts or recent issues
- Search Google Scholar for papers citing SLIDE
- Check ResearchGate or forum discussions
- Verify the current state of the Battery Intelligence Lab website
- Find community feedback beyond what's in the repo itself

The information about citations, community usage, and comparisons draws from training data (up to May 2025) and the extensive documentation within the repository itself. For a complete picture, the user should manually check Google Scholar citations for the 2019 JES paper and the GitHub issues page.

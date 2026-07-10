---
layout: default
title: C++ quickstart
nav_order: 2
---

# C++ quickstart

This program compiles Chen2020 parameters, builds one `nch=12` SPM batch, parses a short experiment, and runs six 10-second intervals. Positive current means discharge.

<!-- doc-test:cpp -->
```cpp
#include <core/Experiment.hpp>
#include <core/ParameterSet.hpp>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main()
{
  using slide::Status;
  namespace core = slide::core;

  core::ParameterSet parameters;
  if (core::ParameterSet::chen2020(parameters) != Status::Success
      || parameters.set("Initial state-of-charge", 0.8) != Status::Success)
    return 1;

  core::SpmFactoryInput input;
  if (parameters.toSpmInput(input) != Status::Success)
    return 2;

  core::SpmBatch batch;
  const core::SpmModelOptions options{ .nch = 12 };
  if (core::buildSpmBatch(input, options, 1, batch) != Status::Success)
    return 3;

  const std::vector<std::string> steps{
    "Discharge at 1 C for 60 seconds",
  };
  core::Experiment experiment;
  core::ParseDiagnostic diagnostic;
  if (core::Experiment::parse(steps, experiment, diagnostic) != Status::Success)
    return 4;

  core::CyclerV2 cycler;
  core::ExperimentSolution solution;
  if (cycler.configure(batch) != Status::Success
      || cycler.run(experiment, 10.0, solution) != Status::Success)
    return 5;
  if (solution.time.size() != 7 || !std::isfinite(solution.voltage.back()))
    return 6;

  std::cout << "samples=" << solution.time.size()
            << " final_voltage=" << solution.voltage.back() << '\n';
}
```

Use SLIDE as a subdirectory while the v4 SDK export is pending:

```cmake
cmake_minimum_required(VERSION 3.31)
project(my_slide_app LANGUAGES CXX)

set(SLIDE_CORE_ONLY ON CACHE BOOL "" FORCE)
add_subdirectory(path/to/SLIDE slide-build)

add_executable(my_slide_app main.cpp)
target_link_libraries(my_slide_app PRIVATE slide_core)
```

The documentation gate extracts the C++ block above into a separate superproject, compiles it against the public target, and runs it.

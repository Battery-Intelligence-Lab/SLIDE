---
layout: default
title: Coding conventions
nav_order: 3
---

## Coding Conventions:

- Do not use "using" in global namespace.
- Include project-specific headers before standard headers. See: https://stackoverflow.com/questions/2762568/c-c-include-header-file-order
- Use descriptive names for variables.
- Break down complex conditionals into multiple statements.
- Minimise nested statements.

### Code editor

We recommend using [Visual Studio Code (VS Code)](https://learn.microsoft.com/dotnet/core/tutorials/with-visual-studio-code). Enable `format on save` to automatically enforce `.clang-format`.

## Copyright
- Verify that the code and libraries do not violate any copyright conditions.
- Our software is licensed under a non-copyleft license, and we aim to avoid incorporating copyleft libraries.
- Report any copyright violations by opening an issue immediately. We will remove any infringing content as soon as possible.
- Using AI tools to help write code is acceptable as long as it is thoroughly reviewed by a human developer and does not infringe on copyrights.


## Versioning

We use semantic versioning with MAJOR.MINOR.PATCH where:

- PATCH: Backwards compatible bug fix without introducing a new feature. 
- MINOR: Backwards compatible new feature. 
- MAJOR: Incompatible changes, especially on API. 

For a detailed tutorial on versioning, see [this video](https://www.youtube.com/watch?v=xvPiZyx0cDc).

## Issue management

As a small library developed by a small team, we review issues informally. To avoid duplicate efforts, please comment on the issue you are working on.

### Commenting: 

We are transitioning to Doxygen-compatible comments. Feel free to update any incompatible comments you encounter. Refer to the [Doxygen guidelines](https://doxygen.nl/manual/docblocks.html) for commenting best practices.
# Changelog {#changelog}

[TOC]

This changelog contains a non-exhaustive list of new features and notable bug-fixes (not all bug-fixes will be listed). 

<br/><br/>
# SLIDE v3.0.0 (aka slide-pack)

## New features


### Battery pack simulation support is added. 
* SLIDE is now now capable of simulation large battery packs. 

### Other updates
* Different battery cell types: Cell_Bucket, Cell_ECM and Cell_PbA are added. 
* [CHANGELOG.md](https://github.com/davidhowey/SLIDE/blob/master/CHANGELOG.md) is added. 
* `*.CMake` files are added for a better build. 
* CPack support for unit tests are added. 

## Notable Bug-fixes

## API changes
* Most functions now use smart pointers. 
* Template / auto parameter deduction is used. 

## Dependencies
  * We require at least CMake 3.17.
  * Now a compiler with C++20 support is required. 

<br/><br/>
# SLIDE v2

Modern version of slide. ([\#8](https://github.com/davidhowey/SLIDE/pull/8)) Deprecated but you can find a copy of Slide V2 in [this branch](https://github.com/davidhowey/SLIDE/tree/SLIDE_v2).

## New features / updates
- Most of the functions, especially `estimateOCVparameters()` should be faster now. 
- SLIDE should be able to run on different platfroms/compilers through CMake.
- Name capitalisation bugs are fixed. 
- Documentation is updated. 
- Auxillary functions and data types are added (such as `fixed_data`, `linspace`, `logspace`, `linstep`, etc.).
- Parallelisation is wrapped into `slide::run` function. 
- Arrangement of files is changed, `data`, `results`, and `docs` folders are created. 
- Modern C++ features are now incorporated into the code. 
  - Instead of increasing stack memory to hold large arrays, now we use `std::vector`. 
  - Most of the raw arrays are changed with `std::array`.
  - `#define` directives are changed with `constexpr` values.
  - `<filesystem>` library is used for handing file/directory operations. 

## Dependencies
  * Required C++ standard is upgraded from C++11 to C++17. 

<br/><br/>
# SLIDE v1

This is the initial release of SLIDE. Deprecated but you can find a copy in [this branch](https://github.com/davidhowey/SLIDE/tree/SLIDE_v1).

## Dependencies
  * A compiler with C++11 support. 
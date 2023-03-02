# Changelog {#changelog}

[TOC]

This changelog contains a non-exhaustive list of new features and notable bug-fixes (not all bug-fixes will be listed). 

<br/><br/>
# SLIDE v3.0.0 (aka slide-pack)

## New features


### Battery pack simulation support is added
* SLIDE is merged with slide-pack. Now, SLIDE is now now capable of simulation large battery packs.

### Other updates
* Different battery cell types: Cell_Bucket, Cell_ECM and Cell_PbA are added. 
* [CHANGELOG.md](https://github.com/davidhowey/SLIDE/blob/master/CHANGELOG.md) is added. 
* `*.CMake` files are added for a better build. 
* CPack support for unit tests are added. 
* `slide::Clock` class is created. `<chrono>` for timing is adapted for correct timing on Mac. 
* `Cell_KokamNMC.cpp` is deleted and it became an header. 
* DEG_ID improved, zero is not counted anymore.
* Free functions to simplify other classes are added: `free::check_Cell_states`
* Catch outside `LiPlating` is removed since it is not needed to be handled since it only throws when id is wrong. 
* `slide.hpp` for including the library is added. 
* `NULL` -> `nullptr`
* `shared_ptr` classes are turned to `unique_ptr`. Copying `shared_ptr` is eliminated. 
* print arguments for some functions are removed (e.g., `V(bool print), getOCV(bool print)`)


### C++20 upgrades
* `std::span` for state assignments. `XYdata_ss` for span. 


## Notable Bug-fixes
* `if(abs(Ii < 0.1))` in Cycler. 
* Typo `bockDegAndTherm` is fixed. 

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
- Most of the functions, especially `estimateOCVparameters()` should be faster now. Unnecessary file reading inside classes is removed. `fitAMnAndStartingPoints` was calling readOCVinput all the time. It is fixed. 
- SLIDE should be able to run on different platfroms/compilers through CMake.
- Name capitalisation bugs are fixed. 
- Documentation is updated. Creating wiki page from `*.pdf`s. 
- False position method for current finding is implemented for SLIDE. 
- For formatting code `clang-format` is adopted. 
- `std::endl` are replaced by `"\n"`. 
- Auxillary functions and data types are added (such as `fixed_data`, `XYdata`, `linspace`, `logspace`, `linstep`, etc.).
- Parallelisation is wrapped into `slide::run` function and `<thread>` library is started being used. 
- File structure is changed, `data`, `results`, and `docs` folders are created. 
- `loadCSV_1col` can take `std::array` or `std::vector` as argument now. It also reuses `loadCSV_mat`
- Binary search (with `std::lower_bound`) for `linInt` is implemented. Detecting if the data is fixed step. Remove linear search from interpolation! O(n) -> O(log(n)) or O(1) for fixed-step data. Nan initialization for interpolation is removed. 
- Used linspace, logspace when searching for parameters! 
- Remove empty destructors etc. Rule of zero. 
- changed pow to scientific notation  `7*pow(10,-10)` ->  `7e-10`. 
- Try-catch block based program flow is mostly removed due to performance penalty. 
- Model Cp, Cn, Vn, Vp, Q are 2D std::arrays, and shorter code. 
- OCVcurves is being moved to a struct. 
- State default constructor is removed since it sets all states to zero, thanks to uniform initialiser. `State() = default`
- State `setStates(State& s)`, `State::setIniStates` are removed. 
- name inputs are deleted from `fitAMnAndStartingPoints`.
- ValidOCV, calculateError, readOCVinput, discharge, fitAMnAndStartingPoints size input removed. 
- `linInt_noexcept`, `discharge_noexcept`, `writeCharacterisationParam`,  `linstep_fix and`, `logstep_fix`, `calcSurfaceConcentration`, `calcDiffConstant`, `calcMolarFlux`, `calcOverPotential` are added. 
- `findCVcurrent_bin` (false position method) is added. 
- `log(x + sqrt(1+x^2))` is changed with `asinh(x)`
- `slide::State` now inherits from `std::array`
- Interpolation is in OCVcurves now. Simplified.  
- New cost function for `determine_OCV`, it is now much faster. 
- RK4 is implemented and FWEuler is seperated. 
- Unit test to control results is written in MATLAB. But only folders, so make also for files. 
- overwriteCharacterisationStates, overwriteGeometricStates are moved to cell, use cell versions. 
- validState is moved to Cell. 
- void `checkModelparam();` is created to check MATLAB param. 
- void `initialise(slide::State &s_ini)` is removed.
- `linInt_matrix` is removed. 
- `oneLevelCharacterisationFit` is removed. Its contents moved to `hierarchicalCharacterisationFit`. 
- A struct created for stress parameters.
- Modern C++ features are now incorporated into the code. 
- Ahpi_vec and Ahni_vec (100k double heap vectors) are removed. `slide::fixed_data` type is used. 
- Reading into double and converting whole vector in integer -> reading into integer.
- Header guards use `#pragma once` now. 
- Change pointers to references whereever possible. 
- printlevels are converted to enum
- Removed many magical numbers:  273 -> `PhyConst::Kelvin`
- Same parameters are defined multiple places are cleaned. For example Kokam parameters.
- paths are moved to settings.
- Data recording system is converted to push_back.
  - Instead of increasing stack memory to hold large arrays, now we use `std::vector`. 
  - Most of the raw arrays are changed with `std::array`; hence, the need for passing number of array elements is eliminated. 
  - `#define` directives are changed with `constexpr` values.
  - `<filesystem>` library is used for handing file/directory operations (`mkdir` and other C-style functions are removed). 
  - TDM-GCC/compiler specific features are eliminated (e.g., dynamic array on stack). 
  - Release mode is added to CMake. 
  - C-based header files are replaced by C++ versions (e.g., `<stdio.h>` to `<cstdio>`)
  - Range-based for loops adopted whereever possible. 

## Notable Bug-fixes
- Using `#define ns` for number of states was causing problems with `<chrono>` library literal `ns`; therefore, it was not possible to compile SLIDE other than TDM-GCC with C++1 standard. Bug was fixed by removing preprocessor directives.  
- Now heap allocation is used instead of stack area expansion. \
- `using namespace std;` from global scope is removed. `std::` is added all necessary functions. 
- warning: variables `j` and `timeCheck` used in loop condition not modified in loop body. `cycler.cpp mode==1, for (int j = 0; j < timeCheck; i++)` They are changed accordingly. 
- Important bug on linux, only use `std::abs` otherwise it may call int `abs(int)`.
- Object slicing due to copying cells is removed. 
- `abs(Icell - I) < pow(10, -10)` -> here `pow(10,-10)` may be ZERO. So `pow(10.0, -10.0)` is better. And `1e-10` is even better. 
- Remove used-defined classes from std namespace (such as state). A namespace called "slide" is created. `State` and `Cell_user` are moved from `std` to `slide`. 
- Others:
  - Typo `CalendarAgeig` -> `CalendarAgeing` fixed. 

## Tried but not implemented: 
- Fixed time step `*.csv` files were tried but no measurable effect in -O0.


## Dependencies
  * Required C++ standard is upgraded from C++11 to C++17. 

<br/><br/>
# SLIDE v1

This is the initial release of SLIDE. Deprecated but you can find a copy in [this branch](https://github.com/davidhowey/SLIDE/tree/SLIDE_v1).

## Dependencies
  * A compiler with C++11 support. 
# SLIDE Completed Items Archive

This file contains completed development items archived from TODO.md for historical reference.

---

## Code-Review by Martin Robinson
- [x] It is both a library and application in between -> Resolved: Now a library
- [x] Catch2 library for testing - Implemented
- [x] codecov - Integrated
- [x] Doxygen - Configured

## Open Issues (Resolved)
- [x] sp - sn question - Resolved
- [x] cps - cns question: dOCV_neg and dOCV_tot are taking zp_surf as input due to our data - Documented

## Current Priority (Completed)
- [x] Cycler `getSafetyVmin` and `getSafetyVmax` removed for redundancy and expensiveness
- [x] Cycler `initialise` removed, now constructor is used
- [x] nOnce scaling removed; will be added in future if necessary
- [x] check_safety removed since VMAX and VMIN are checked
- [x] regulate diagnostic variable removed; now act on individual cell limits
- [x] Instead of Battery class holding converter pointer, it is now holding a converter object
- [x] Cannot set current for Module_p<Module_p> - Fixed (see `test_Hierarchichal_p`)
- [x] Some tests were hard-coded for Cell_type which has been changed to template
- [x] Module_p::getRtot() check for contact resistances removed for simplification
- [x] V(bool print), getOCV(bool print) print argument removed
- [x] dV > settings::MODULE_P_V_ABSTOL && dV > settings::MODULE_P_V_RELTOL * (*V_max_it) -> Changed && to ||
- [x] getVi recalculation inefficiency fixed - Removed and getVall added
- [x] Module_p mean voltage calculation fixed (was causing 3.6 to 3.58 voltage difference)
- [x] redistributeCurrent inside rebalance deleted since not required

## From SLIDE v2 (Completed)
- [x] Both Cell setT and State setT checks if temperature is between limits - Cell setT removed
- [x] Error in BasicCycler converted to warning with time step adjustment (Throw 1003 eliminated)

## slide_pack Changes (Completed)
- [x] Cycler.CC was charging a bit if already satisfied voltage limit - Fixed
- [x] Module shared_ptr<StorageUnit> SUs[MODULE_NSUs_MAX] -> Converted to vector
- [x] checkUp_getCells and checkUp_getModules deleted and replaced with visitors

## New Features (Completed)
- [x] SmallVector added. DegArray now derived from SmallVector
- [x] if (succ != 1) after CV phase changed with more meaningful limit reaching condition
- [x] validSUs functions combined
- [x] getVariations() deleted since not necessary to hold these variables inside cells
- [x] storeData() pattern removed
- [x] removing getVi for series module and getSUVoltages for all
- [x] getSUCurrents removed - Was only used in module_p for a very short function

## Architecture Changes (Completed)
- [x] All shared pointers converted to unique pointers (shared_ptr -> unique_ptr)
- [x] For setSUs and validSUs span is used
- [x] OCVcurves -> OCVcurves* (12 seconds to 0.22 seconds improvement)

## Header-Only Library Progress (Completed)
- [x] Cell, Cell_Bucket, Cell_ECM: Made header-only (80 kB -> 77 kB)
- [x] Model_SPM, State_SPM: Made header-only
- [x] interpolation.cpp, read_CSVfiles.cpp, slide_aux.cpp, util.cpp, util_error.cpp (77 kB -> 76 kB)

---

*Last updated: 2026-01-29*
*Archived from TODO.md during codebase reorganization*

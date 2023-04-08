# Information about this branch of this library: 
- Only a test version, not complete. 

# Disclaimer: 
- This TODO list is pretty much informal text of what is in my mind. 

## TODO: 


### Code-review by Martin Robinson: 
- [ ] It is both a library and application in between. 
  - [ ] One way: have a config file model file and options. 
  - [ ] Brady: Julia wrapper? 
  - [ ] Martin: Lot more Python audience. 
  - [ ] Class naming inconsistent: standard format. Capital letter Camel case. 
  - [ ] There are classes with completely capital? Which one? 
  - [ ] SUNDIALS is pretty quick. 
  - [ ] Automated tests against PyBAMM 
  - [ ] Catch2 library for testing. 
  - [ ] Snapshot testing. 
  - [ ] Template projects in OxRSE template-project-cpp
  - [ ] gui_starter_library Jason Turner. 
  - [ ] clang-tidy 
  - [ ] codecov. 
  - [ ] Doxygen. 
  - [ ] CPack


### Comparing against slide-pack:
- [ ] CC, CV, CCCV functions seem to be working. However CV is not doing the intended thing for both SLIDE and slide-pack.
- [ ] paperCode::Vequalisation seems to be working. 
- [ ] paperCode::thermalModel does not work for slide-pack.

### Open issues: 
- [x] sp - sn question 
- [x] cps - cns question: dOCV_neg and dOCV_tot are taking zp_surf as input due to our data. Very important. 
- [ ] Should capacity check also contain a CV phase? 


### Current priority:
- [ ] Unnecessary printing statements with `prdet` is removed to reduce cluttering. 
- [ ] check_safety seems useless. Why do we have that when we have VMIN and VMAX? 
- [ ] regulate diagnostic variable: if diagnostic on the individual cell limits are respected. Otherwise system level. 
- [ ] T_MODEL and T_ENV etc. should not be constants!!!!
- [ ] Add a GITT function. 
- [ ] Instead of Battery class holding converter pointer, it is now holding a converter object. 
- [ ] Test cases testing shared_ptr logic (i.e., testing if the pointer points to another thing when Cell is changed.) should be removed because when you change cells the previous cell is removed due to unique_ptr logic. 
- [ ] std::span<double> & for set states do not work when a vector given. Need a better idea to set and get. Maybe another class. 
- [ ] Cannot set current for Module_p<Module_p>. See `test_Hierarchichal_p`. Very important.
- [ ] copy functions are commented out. 
- [ ] validSUs actually not necessary! Removing. Write free functions to check voltage/current equality!
- [ ] T_MODEL == 2 causes thermal runaway in test_CyclerVariations_high
- [ ] redistributeCurrent_new and other iterative algorithms require many iterations (up to 2500)!
- [ ] test_Cycler_CoolSystem passes only when T_MODEL==2
- [ ] Some tests were hard-coded for Cell_type which has been changed to template.
- [ ] operator[] is added to Module class to reach SUs. 
- [ ] getNSUs() -> size() so module size should be number of SUs.
- [x] Module_p::getRtot() check for contact resistances is removed for simplification. 
- [ ] Why do we need getRtot? Only for algorithms? Maybe getThevenin would be better. 
- [x] V(bool print),  getOCV(bool print) print argument is removed. 
- [ ] therm.Qcontact is not calculated properly it is 6x overestimated! 
- [x] dV > settings::MODULE_P_V_ABSTOL && dV > settings::MODULE_P_V_RELTOL * (*V_max_it) ->  && to ||
- [x] getVi is calculating the SU voltage seen at the terminal by subtracting voltage drops on resistances. However, since it is called multiple times it recalculates everyting multiple times. Therefore, not very efficient. It is removed and getVall added. 
- [x] Previously Module_p was calculating the mean voltage; since voltage of all cells were not reliable; however it is not needed. Also it was wrong causing 3.6 to 3.58 voltage difference! Now, it is SUs[0]->V() - I() * Rcontact[0];
- [x] redistributeCurrent inside rebalance is deleted since it should not be required since all have same voltage. It will be eliminated completely soon. 
- [ ] setVoltage function is being added for a better CV period. 
- [x] "#if TIMING" is removed, profiler should be used if needed. 
- [x] "setI_iterative" is removed. 
- [ ] redistributeCurrent() -> PI Control does not work well causing high error in current. 
- [ ] Add snapshot tests for 
- [ ] SOC/Temperature dependent RC pairs for ECM. 
- [x] Making Cell_ECM template to remove Cell_Bucket.
- [x] Bugfix: ECM had - in the equation, corrected. 
- [ ] Fixed data function argument should be reconsidered! 
- [ ] Status member functions like  status.good()
- [ ] It should be decided if we throw an error in interpolation or not for testing invalid states. 
- [ ] Consider using std::variant for some data types. 
- [ ] Procedure: 
  - [ ] Markers to incidate end of this section deleted. As well as seperators. 
  - [x] checkUp_getCells and checkUp_getModules are deleted and replaced with visitors. 
  - [ ] Reduce dynamic pointer cast in Procedure and let Polymorphism work. 
  - [ ] Capacity checking protocol in Procedure::CheckUp is removed.  
- [x] TIME_INF is created for a large time value. 
- [x] getDaiStress is simplified by removing unnecessary R multiplication and division. 
- [ ] State classes should have const members. 
- [ ] Make more methods const : Vmin(), Vmax(), VMIN(), VMAX(), Cap()
- [ ] Stressparam vs. classes may encapsulate calculating some things. 
- [x] cube and sqr functions added for utility. 
- [ ] Module requires number of cells to construct which is unnecessary. 
- [ ] Cycler CV is not working???? CCCV works if it does not have current beforehand it does not work... 
- [ ] Create a small class config for limits reached. It should have CheckLimits.
- [ ] Some variables like ncheck and nbal were defined double but they are int. 
- [ ] Instead of taking unique_ptr or raw ptr, functions should take object references. 
- [ ] Literal operators for units are being added. 
- [ ] Add CTest support for tests
- [ ] Add static analysers: include-what-you-use, valgrind, etc. 
- [ ] Doxygen integration - adapting commenting style. 
- [ ] CPack and installation improvements. 
- [ ] Figure out why OBJECT libraries cause linker error. We converted everything to STATIC
- [ ] Check #CHECK and #TODO tags in the code.  
- [ ] Add variable data storage (Vk is working on)
  - [ ] Create enum for storable data. 
  - [ ] Why are we using different vector<struct> anyway? We can use just a big vector + deserializer? 
  - [ ] Add a generic state term to cover time, Ah, Wh. 
  - [ ] In future it is possible to use enums as a counter for vector based states. 
- [ ] Probably there is a bug in slide-pack where time, Ah, Wh values are not resetted after a throw. Try two cell in series config with one has smaller capacity so it reaches its capacity earlier. Then charge with 1C + large time step. Probably first cell will store Ah and second won't. and there will be difference in their Ah. Or even it does not, it will be different than real Ah. Just charge and discharge. 
- [ ] Cycler should be able to take unique_ptr and convert. 
- [ ] setSUs and assigning unique pointers then testing individually is very difficult. Clearly a design problem.
- [ ] getNSUs may slow down time to time. 
- [ ] std::vector<double> Iolds in Module_p.cpp


### From SLIDE v2: 

- [ ] Convert strings in function parameters to string references. 
- [ ] Also add a particular data type for data, so warmstart and interpolators can be hold. 
- [ ] Inline small functions. 

- [ ] Why does getstates call the setstates twice?
- [ ] Convert if(x<x_min); if(x>x_max) to   if (x<x_min) else if (x>x_max); So less loop. 
- [ ] Use more initializer lists: 
  - [ ] State constructor now uses initializer list. 
- [ ] Some instances of slide states are converted to pass-by-reference. 
- [ ] Model.Dn and Dp are created as nch elements but require nch+1 elements! 
- [x] Both Cell setT and State setT checks if temperature is between limits. So cell setT is removed. 
- [ ] There are some places where double is compared to zero. Such as: 
- [ ] Cycler followI does not require "DegradationData_batteryState.csv" and "DegradationData_OCV.csv" but it is created anyway. 
- [ ] Convert printing function to template function.  -> Did not work. Be careful about string initialisation.
- [ ] BasicCycler.cpp is very large. 
- [ ] Convert verbose to settings::verbose? -> May reconvert. 
- [ ] Create a special class for data so you can use a warmstart. 
- [ ] But need to shift all paths to bin/debug etc. Maybe a binary path.
- [ ] Remove physical constants etc. non-battery dependent values outside of the class. 
- [ ] Do not use exceptions for normal system operation: example BasicCycler::CC_halfCell_full
- [ ] When using vector and array, remove unnecessary size information from function arguments. 
- [ ] Instead of using references, use structural binding and std::pair or tuple to return multiple elements. 
- [ ] getters should be const reference and try to add const to all non-changing functions. 
- [ ] create mutable getters as private functions for internal usage and const ref getters for external usage. 
- [ ] Convert states to array representation. 
- [ ] States has initstates which is checked all the time. Move it into cell. It is not used except validstates 
- [ ] Make validState() bool not void.
- [ ] Write for loop to simplify validState
- [ ] To make OCVcurves std::array, they have to appear on derived classes. 
- [ ] Why estimateOCVparameters() cannot find any parameters? 
- [ ] for (int ap = 0; ap < nAMstep / 3 + 1; ap++) in determine_OCV.cpp with ap3 = 3 * ap + 2 means  nAMstep+2 steps. is it correct? 
- [ ] Look at .m files for problems like _SOC. 
- [ ] Write unit tests.
- [ ] For starting and terminating things. 
- [ ] Automatise SOC -> OCV conversion (see: for other cells the user has to derive the conversion from the OCV curve)
- [ ] Write parser for cycling things. 
- [ ] In parallel, output texts are mixed. Write a string stream to pass. 
- [ ] Use actual error classes and predefine errors. 
- [ ] Python and MATLAB interfaces. 
- [ ] Can we use __func__ to substitute function names in some warnings? 
- [ ] Cell::setVlimits does not control VMIN < VMAX. 
- [ ] determine_characterisation and cycling both has "CCCV" functions. 
- [ ] BasicCycler::returnCyclingData try to see if we can get away with constant reference. -> yes we can. 
- [ ] Important: Cycler copies cell, so changing cell state probably does not work. Check characterisation. -> It is cell parameters. 
- [ ] Documentation of the new stuff. 
- [ ] In determine_characterisation unnecessary try-catch when error calculation. Flag is added. Also it was not breaking the loop when there is large error in one data set. 
- [ ] setCyclingDataTimeResolution does not throw anymore, it just makes time resolution 0 if it is negative.
- [ ] Actually nrStep+3 steps are taken! Therefore, nr 6 -> 9 (also same in AMp -> +3 iterations!)
- [ ] Improve const correctness.
- [ ] Discharge functor. 
- [ ] In determineOCV, discharge function calculates sp, sn and Ah one more time before saving. 
So they are one step ahead of ocvpi, ocvni and V. Is it a bug? Or is it saving exceeding values with correct voltages? (Original code) 
- [ ] Iterator for fixed data is implemented but not completely working. 
- [ ] static thread_local keywords are being used. 
- [ ] More range-based for loops.
- [ ] prev(T x_current) and next(T x_current) are added to fixed_data to satisfy generic interface. 
- [ ] Add package manager installation option. 
- [ ] Why is integer T is used in Cycler::checkUp_pulse? 
- [ ] Higher order integration Runge Kutta 4 and adaptive time stepping. -> Fixed for slide
- [ ] The fitting algorithm does not currently support fitting the activation energy directly.
- [ ] "check what would happen if the current is set to 0 for one time step" why? no needed. 
- [ ] better findCurrent algorithm for slide-pack
- [ ] There is not much difference between setting I and having small time step and having dt + dt_I time step because we do not store I. 
- [ ] Check if sparam.s_dai_p and others are correctly used in multi-step integration. 
- [ ] updateDaiStress and updateLaresgoitiStress may not need to be calculated if blockDegradation is true.
- [x] Error in BasicCycler::.... The total time is not a multiple of the time step dt, -> Converted to warning. And time step adjustment is added. Throw 1003 is eliminated.
- [ ] BasicCycler::storeResults now does not check if it is storing the same point. Therefore, you should check. 
- [ ] When creating xy data put something not to check if it is fixed. 
- [ ] determineOCV is still inefficient, we may reduce search space lot more by looking at the sp and sn requirements. 
- [ ] Explation function for all classes? 

slide_pack integration:
- [x] SPMModel.h, SPMModel.cpp

### slide_pack changes: 
- [x] Cycler.CC was charging a bit if already satisfied voltage limit is given.
- [ ] StorageUnit parent shared_pointer -> raw poiner. 
- [x] Module shared_ptr<StorageUnit> SUs[MODULE_NSUs_MAX] -> vector.
- [ ] Rcontact is removed from Module and it will be moved into Module_s and Module_p
- [ ] (verbose_gl > v_noncrit_gl)  ->  (settings::verbose >= printLevel::printNonCrit) 
- [ ] (verbose_gl >= v_crit_gl)   (settings::verbose >= printLevel::printCrit)
- [ ] #if DATA_STORE is being removed and module is now a template. 
- [x] getIndex seems to be useless. 
- [ ] getting Vmax everytime is a problem. 
- [ ] getNstates is not needed. It is only used for checking sizes and it won't be needed. 
- [ ] Use default instead of empty destructor. 
- [ ] If base class does what needs to do, do not redefine everything. 
- [ ] Histogram should be in another class. 
- [ ] DEG_ID and OCVcurves may need to be shared pointers
- [ ] More hierarchy for ECM and SPM 
- [ ] Factory methods. 
- [ ] Add MSVC things 
- [ ] Can we make ID static since it is same in all classes? Do we change? Look at it 
- [ ] StorateUnit -> StorageUnit
- [ ] getSUTemperatures in Module creates unnecessary storage then -> getSUTemperature(i).
- [ ] typeid(*SUs[i]) == typeid(Module) || typeid(*SUs[i]) == typeid(slide::Module_p) || typeid(*SUs[i]) == typeid(Module_s) -> is not good dynamic_pointer_cast already gives nullptr. 
- [ ] All shared pointers are converted to unique pointers. 
- [ ] validStates should not backup/should not change states! SetStates should not back up but can change. 
- [ ] getNstates is now useless. Removed from all functions. 
- [ ] Cell_Bucket setSOC was repeating some code. Therefore, we combined it with setStates. However, the logic should be divided better.
- [ ] setStates to support r-values. 
- [ ] do not copy shared pointers in for loops.
- [ ] Why there are so many ifs instead of else-ifs
- [x] validSUs functions are combined. 
- [ ] Cycler was checking if storeData is full but storeData itself should check and call writeData. Or better, Cycler should collect everything and write. 
- [ ] su->getIndex() should not be used  in cycler because it is removed. 
- [ ] timeStep_CC_i -> seems to be used only for binding; therefore, not needed anymore. 
- [ ] par variable inside module should be constexpr somehow. 
- [ ] holding a vector of instances in a class? 
- [ ] Creating a thermal model. 
- [ ] Bugfix: Module_s::getI() tries to return value if (getNSUs() >= 0) but it should be >
- [x] getSUCurrents is removed. It was only used in module_p for a very short function. 
- [ ] Although not very necessary, eliminating internal getNSUs() and size().
- [ ] Add small free / helper functions. 
-   [ ] mean. 
-   [x] is_zero.
- [ ] getSUVoltages can replace some of the functions. 
- [ ] getSUTemperature does not seem to be useful, should we allow indiced calling? 
- [ ] Why thermal model takes 1-dim arrays? 
- [ ] Writing functions must be outside of Procedure. 
- [ ] Intermediate design decision: setStates -> span,  getStates -> vector. Also hold the states as `std::array<double,N>` in actual object. In future we can allocate alltogether. 
- [ ] EmptyStorageUnit to remove if(nullptr) parent? 
- [ ] Both soft and hard limits for battery is not good. 
- [ ] Some of the vectors can be eliminated 
-   [ ] via small vector optimisation. 
-   [ ] via unique_ptr(T[]). Use span for getting 
- [ ] Larger member variables should defined first in classes due to padding.
- [ ] use #include `<source_location>` to simplify error/warning/explanation messages.
- [ ] Delegated constructors? 
- [ ] Make this 'verb' things compile time things!!!!
- [ ] do not put get to everything   getV() -> V(),  getT() -> T()
- [ ] getVcheck -> checkV. 
- [x] isCharging() -> I() < 0  and isDischarging() -> I() > 0 are added to enforce consistent current representation. 
- [ ] setI should not return voltage, setI should not throw for no reason.
- [ ] why setI returns voltage? it should not. It is also not used anywhere. 
- [ ] setI -> setCurrent
- [ ] output of the setstates is also not used!!!! 
- [ ] error codes should be universal. Functions throw or return random codes. 
- [ ] cycler has set current and what we did setI is mixed with it :O 
- [ ] v < VMIN() && isDischarging()  this if else statements can be converted into summation and multiplication. 
- [ ] viewStates, viewVariations returns span or some other lightweight object. 
- [ ] CellData, CellDataStorage and CellDataWriter are created. We need policy design. 
- [ ] Cell::checkVoltage does not throw anymore. We encapsulate it in try-catch instead of in upper function. 
- [ ] getCSurf to return flag + pair? Optional?
- [ ] st.setZ can have a move object. 
- [ ] why setCS does not set Vcell_valid=false? 
- [ ] Don't forget to add timing functions into some particular functions. 
- [ ] Why case2 of CS starting with csparam.CS2alpha is different in old and new models? Also Case 1?
- [ ] Why case1 of LAM is different?
- [ ] LAM and CS are changed to be same as new versions.
- [ ] double degState[CELL_NSTATE_MAX] is useless. Causes use of 700 kB unnecessary memory. 
- [ ] Create another class for thermal model.
- [ ] What happens considering only one electrode? 
- [ ] Why Qrev definition is different? In new one it is -I() * T() * dOCV  but old one -i_app / L * st.T() * dOCV, because of change in units? 
- [ ] SetSOC does not change any concentration in SPM cell. However, it is effective in others. So is it an estimation thing or not? 
- [ ] Why does getOCV does not include entropic coefficient? 
- [ ] Use template functions instead of repeating st dependent but same functions. 
- [ ] setCurrent output is set to int. 
- [ ] Multiple CMakeLists is created which also allows me to find more path errors.
- [x] Geometry_SPM is created. 
- [ ] DEG_ID variables are changed from int to uint_fast8_t -> 144 byte to 36 byte size reduction. We also need to improve print function.
- [ ] DegArray is added to DEG_ID to remove deg_id.SEI_n > deg_id.len type of controls at each step.
- [ ] Degradation ID's should be enums. 
- [ ] Cannot use const with DEG_ID
- [ ] Range based for loops are added. 
- [ ] SmallArray is required. 
- [ ] i < getNSUs() - 1 is a bug! Fix immediately. 
- [ ] setStates validStates all calling each other copying states for nothing! We need to remove redundancy. We need a state pointer thing to do that.
- [ ] Remember that in copy(), usage stats are not copied across to the new module.
		 * So this function will write all zeros if called with a copy of the SU   -> is this a bug or feature? 

- [x] for setSUs and validSUs span is used. 
- [ ] We are reading file lots of times when individually creating cells! In release mode it takes 30 seconds. 
  - [x] slide::Model_SPM -> slide::Model_SPM*  and makeModel();  from 34 seconds to 12 seconds.  18900x326 double = 47 MB RAM is also saved.
  - [ ] OCVcurves -> OCVcurves*   12 seconds to 0.22 seconds. Do not forget things otherwise than NMC. 
- [ ] Benchmarking?
  - [ ] Now EPFL battery 1 hour timestep CC is like 2 seconds for and 3.5 with ageing. 
- [ ] su->getIndex() in Cycler is not necessary. Cycler should keep its index. 
- [ ] Why cycler CC is much slower then time step CC ? 
- [ ] Not necessary constants to remove:
  - [ ] CELL_NCH, CELL_NOCV_MAX
- [x] remainder for integers are eliminated. 
- [x] storeData() pattern is removed.  
- [ ] Ncells for module is created.
- [ ] etacell_valid -> is removed. Should be checked in future if it really matters for performance. 
- [x] Subfolders should also include linker options. 
- [x] removing getVi for series module. and getSUVoltages for all. 
- [ ] Consider making test functions friend and getVi protected. 
- [ ] "HVAC coolsystem for active cooling with the environment" obligation should be  removed. 
- [ ] Minor updates: 
  - [x] i in Cycler::rest is removed. 
- [ ] VMAX and VMIN in modules calculated VMAX and VMIN by summing then checks if each cells has a violation. That would give the same result. 
- [ ] We need to have constant values for VMIN, VMAX for all SUs which are initialised at beginning. 
- [ ] #CHECK how module calculates isCharging? 
- [ ] Add typename to all classes for printing purposes. 
- [ ] checkCurrent only checks voltage!
- [ ] SetCurrent has unusual error codes! 
- [ ] Gradually setting current should be there. 
- [ ] Cycler checks Vlow many times. Maybe it is better to cache? 
- [ ] Unnecessary voltage check in SetCurrent of module_s is removed. 
- [ ] SEI, CS, LAM, LiPlating, getDaiStress are used once so can make them inline. 
- [ ] Cell_SPM::setSOC does not actually set SOC.
- [ ] Check why in Battery.cpp Tbatt check is repeated. Should Tbatt be something else? 
- [ ] Cell::setT(Tnew) does not check temperature anymore, be careful! 
- [ ] redistributeI returns number of iteration but sometimes it is asked for voltage. Minor bug. Its number of iteration is only used for test.
- [ ] Check issues Jorn listed in redistributeI()
- [ ] Try to find a way to solve all parallel circuits. 
- [ ] #CHECK getVi for parallel. 
- 
- [x] SmallVector is added. DegArray is now derived from SmallVector.
- [x] double Rdc is removed from Cell.hpp
- [ ] Memoize Cap. 
- [x] Bugfix: dt / 3600.0;  should be dti / 3600.0;
- [x] if (succ != 1) after CV phase is changed with more meaningful limit reaching condition. 
- [ ] Change const string& with string_view.
- [ ] Why do we cheeck Vini in Cycler::setCurrent? 
- [ ] Reached Voltage limit for CC and reached current limit for CV were both same  = 1 so we distinguised. 
- [ ] Should we include entropic effect in OCV or not? 
- [ ] SOC -> why do we use columb counting? 
- [ ] Definitely create a file type to compactly save files and retrieve. 
- [ ] Remove Error IDs.xlsx 
- [ ] "${CMAKE_CURRENT_LIST_DIR}/" is mostly eliminated since it is not needed in newer CMake versions. 
- [ ] "develop" folder is added for developer-related matters. 
- [ ] License files of individual libraries are moved into the corresponding folders. 
- [ ] tests folder is created and unit tests are moved into that folder. 
- [ ] separator variables in checkUp functions seem to be unnecessary; therefore, removing. 
- [ ] Changing constructor delegation. Module_p and Module_s constructors are combined. 
- [ ] ID should be unique. 
- [ ] Module functions are being combined. 
- [ ] std::algorithms and free functions for modules.
  - [ ] transform_sum is added. 
- [ ] storeData(getNcells()) pattern is not good. 
- [x] getVariations() is deleted since it is not necessary to hold these variables inside cells. 
- [x] For some reason slide-pack had different LAM parameters (new fitting?). I will test the difference. Because: there are different set of fitting data.
- [ ] Cell_SPM should not hold all ageing model parameters. 
- [ ] Determine sn and AMp from boundary conditions. They are very sensitive. 

### C++20 changes (yay!):
- [x] std::span for state assignments. 

### Some new ideas to implement: 
- [ ] For XY data read to vector but then create a specific-sized data structure with all heap allocated as if make_shared.
- [ ] begin and end functions for StorageUnit to traverse the children. 
- [ ] Classes to hold static vector of their elements for make_X;
- [ ] Constructor chaining. 
- [ ] Free functions to call member functions. 
- [ ] Status class to hold error codes. 
- [ ] std variant with regular pointer and unique pointer OR a boolean to indicate deleter. 
- [ ] Cycler kind of things should be able to take things other than SU pointer. A template pointer could make things faster. But let's see. 
- [ ] Making SLIDE a header-only library for easy compilation. (Maybe use a proper CMake config?)
  - [x] Cell, Cell_Bucket, Cell_ECM: Instead of increase in *.exe, there is a decrease. 80 kB -> 77 kB
  - [x] Model_SPM, State_SPM. 
  - [x] interpolation.cpp, read_CSVfiles.cpp, slide_aux.cpp, util.cpp, util_error.cpp -> 77 kB -> 76 kB
- [ ] addData parameter in timeStep_CC functions is now useless since it only stores throughputs which are now included in states. So it is removed. 
- [ ] Threadlocal vectors for dynamic-sized stack for some functions.

### Formatting: 
- [ ] Configure clang-format, cmake-format etc. 

### Developer changes: 
- [ ] CMake folder and some files are added. 

### JOSS: 
- [ ] Added JOSS folder and Github workflow. 




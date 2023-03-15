---
layout: default
title: Degradation
nav_order: 3
---


# Degradation
## New degradation model
If you find a new degradation model for an existing mechanism (e.g. a new SEI-growth model) and want to implement it, follow these steps:

- Go to the struct with the parameters for that mechanism in `cell.hpp` (e.g. SEIparam).
- Add the definitions of the fitting parameters you want for that mechanism (e.g. sei4k with a rate constant). Set the value in the constructors of all the child classes of `cell.hpp`.
- If you need new physical parameters (e.g. the volume fraction of the SEI layer), it is recommended you code them in the Cell class (e.g. `e_sei`). You can also do this in the struct if you want but it is recommended to put the physical constant in the Cell class directly. Add the parameter definition in the class definition in `cell.hpp` and give its value in the constructors of all the child classes of `cell.hpp`.
- Go to the relevant function in `cell.cpp` (e.g. SEI)
- Add a new case in the switch deciding which model to use (based on `deg_id.SEI_id[i]`). There you can code your new model (e.g. as case 4), using the newly defined fitting parameters (e.g. `seiparam.sei4k`) and the newly defined physical cell constants (e.g. `e_sei`). The battery states are in the field s of the `Cell` (so e.g. the cell’s temperature can be accessed by s.getT()). If this model affects a battery state which so far was not affected by the degradation mechanism, this extra state must become an output parameter, do the following:
    - Adapt the function definition in the `*.hpp` and `*.cpp` file, e.g. add `double* thickn` as output if you want to decrease the electrode thickness due to SEI growth.
    - You then have to set the value of this new output variable in your new model, and set it to 0 for all other models.
    - In `dstate`, when the degradation mechanism is called, add the extra output variable, e.g. `double dthickn_sei`; and SEI(`OCVnt`, ..., `&dthickn_sei`)
    - Then at the last code block of `dstate` (where the time derivatives are calculated), you have to add the effects from the new degradation models (e.g. `dstates[2*nch + 4] = dthickn_sei + dthickn;`).
- It is a good idea to update the comments for the degradation models (the `DEG_ID` struct is defined in `cell.hpp` and its model identifiers are repeated in main in `main.cpp`).
- Then you can simply select this new model in the main function in `main.cpp` (e.g. `deg.SEI_id[0] = 4`) and the code will use this new degradation model.


## New degradation mechanism

There are currently 4 mechanisms implemented (SEI growth, surface cracking, LAM and lithium plating). This is a purely *logical* split (i.e. you can put a model for SEI growth in the function which normally does LAM growth), so in theory you can add your new model for the new mechanism under one of the existing mechanisms (see the section [New degradation model](#new-degradation-model)). This is however a bit messy, so it is better to add a fully new degradation mechanism (e.g. binder decomposition or BD). To do that, follow these steps:

- In `cell.hpp`, go to the definition of the struct `DEG_ID`. Add the identifiers for the models for this mechanism (e.g. `int BD[10]` which will allow you to use 10 models for this new mechanism at the same time and `int BD_n`).
- In main in `main.cpp` add values for the new fields in `DEG_ID` you have made (e.g. `deg.BD[0] = 1;` and `deg.BD[1] = 4;` and `deg.BD_n = 2;` to use models 1 and 4 for binder decomposition)
- In the Cell class you will have to define a new function which will implement the various models for this new mechanism (e.g. `BD`). You have to add this function definition in the header file and implement it in the `*.cpp` file. You can give the function the inputs you need. The outputs must be the time derivatives of the battery states on which this degradation model has effect (e.g. `double* den`, `double * dep` if you want to change the volume fraction of active material in both the anode and cathode).
- Implement the models for this new mechanism. You can follow the same structure as the existing function. E.g. using a loop for each model to use (from `0` to `deg_id.BD_n`) and a switch (on `deg_id.BD[i]`) to implement the various models.
- In the function `dstate` in `cell.cpp` you have to call the new degradation function (e.g. `double dep_bd, den_bd;` and `BD(&den_bd, &dep_bd`). Then at the last code block of `dstate` (where the time derivatives are calculated), you have to add the effects from the new degradation models (e.g. `dstates[2*nch + 5] = dep + dep_bd;` and `dstates[2*nch + 6] = den + den_sei + den_bd;`).

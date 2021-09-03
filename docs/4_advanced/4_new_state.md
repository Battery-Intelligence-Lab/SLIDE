---
layout: default
title: State
nav_order: 4
---


# Adding a new state-variable

This is a very complicated procedure, please only do this if you understand the code.

A battery has a certain number of states, which are grouped in the class State. State-variables are all independent parameters of Cell which vary over time (things like the voltage and stress are not independent, they can be calculated if the battery states are known). If you want to implement a new model or a new mechanism, you might have to add new state-variables. If it are totally new variables, you can add them directly. In most cases, the variables will be constants in the `Cell` class and they have to be removed there first.

In the example used below, you want to implement a new degradation model which will decrease the rate constant of the li-insertion reaction at the graphite (`kn`) over time (e.g. due to SEI growth).

In the example used below, you want to implement a new degradation model which will decrease the rate constant of the li-insertion reaction at the graphite (`kn`) over time (e.g. due to SEI growth).
- On the top of `state.hpp`, increase the value of `ns` to reflect the new number of state variables e.g. to (`2*nch+15`)
- In the class definition of State in `state.hpp` add the following
    - A field with the variable in the private part of State, e.g. `double kn`
    - A getter and setter in the public part of State, e.g. `double getKn();` and `void setKn();`
    - Add an extra input parameter in the function `iniStates`, e.g. `double kni`
- In the class implementation in `state.cpp`, add the following
    - In the constructor, `State()`, initialise the new state to `0`
    - In the function `iniStates(..)` add the new input parameter in the function definition (`double kni`) and set the value of the state-variable (`kn = kni`)
    - In the function `getStates(int nin, double states[])` add the new variable in the output array, e.g. states[2*nch+14] = kn). Also update the assert-statement which checks the number of states (this is done such that you don’t forget to add the state variable here), e.g. assert(ns == 2*nch + 15)
    - Add the getter and setter functions which returns / set the value of kn
    - In the function `setStates(double states[])`, add the new variable, e.g. `kn = states[2*nch+14]`.
    - In the function `setStates(State si)` add the new variable, e.g. `kn = si.getKn()`.
    - In the function `validState()` add a check on the allowed values, e.g. `bool knmin = kn < 0`, print an error message if the value is not allowed and update the condition when the error is thrown, e.g. `if (Tmin || .. || knmin)`
- Remove the definition of kn in the class definition of `Cell`, in `cell.hpp` (at this point a lot of errors will appear in `cell.cpp`). If your newly added state was not a class-variable you can skip this step
- Go to the implementation of the Cell class in cell.cpp and do the following
    - Everywhere where the original class-variable kn was called, there is now an error. Replace the calls to kn with the getter for the new state `s.getKn()`. If your newly added state was not a class-variable you can skip this step
    - Add the degradation model/mechanism which affects the new state as explained in the previous chapters, e.g. add SEI growth model 4 with the model for the decrease in kn.
        - If the decrease is conditional on the amount of SEI growth predicted by the other models, you can add it as a new case of `deg_id.SEI_porosity`, e.g. else if `deg_id.SEI_porosity == 1 {*dkn = seiparam.sei_newParameter * *isei}` (and you can make a `case == 2` with both cases 0 and 1 if needed, or you an add a new identifyer for the decrease in `kn` in `DEG_ID`).
        - You will have to add new output values for `kn` in the degradation mechanism which will decrease kn, see the section [New degradation model](../4_advanced/3_new_degradation.html).
    - In `dstate` you have to add the time derivatives of this new state at two locations
        - In the code block *if we ignore degradation in this time step, we have calculated everything we need* you have to set the time derivative of the new state to 0 (assuming the new state is only affected by degradation, if it is a cycling state then you have to give the time derivative here) , e.g. `states[2*nch+14] = 0`.
        - In the final code block titled ‘time derivatives’ you have to give the value of the time derivative for the new state (if it was a cycling-state then you have to give the same value as before; if it is a degradation state then here you give the nonzero time derivative) , e.g. `states[2*nch+14] = dkn`. Also update the assert-statement which checks the number of states
- In the function `checkUp_batteryStates(...)` in `cycler.cpp` the value of the new state will be automatically recorded. I.e. in all the output `*.csv` files, there will be an extra column with the evolution of the new state parameter, e.g. the one but last column in the `*.csv` file shows how `kn` decreases over time. You don’t have to change anything here, it’s mentioned so you know how the new state will be recorded.
- In the MATLAB script `readAgeing_BatteryState.m` you have to account for this extra column (which is not going to be the last column, as you can see in `checkUp_batteryStates(...)`, the last column is the battery resistance while the new state variable is the one but last column). E.g. in MATLAB, add `state{i}.kn = A(:,30);` and `state{i}.R = A(:,30);`. If you want to display the evolution of the new state variable, you have to add a subplot in the figure (or make a new figure)
- In main in `main.cpp`, call the degradation model/mechanism which will affect the new state and simulate it.

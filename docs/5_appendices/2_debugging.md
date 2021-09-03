---
layout: default
title: Debugging & Common errors
nav_order: 2
---

# Cycling

This section gives some advice when errors appear in the code.

The code is written with a defensive programming style. This means that exceptions are thrown in case of unexpected events (the `throw` statements in the code). When calling functions which might throw errors, the code is usually wrapped in a try-catch block ( try {...} catch(int e){...}) to handle the errors. Some errors can be corrected by the calling function, while others can’t and are thrown on to the higher functions. When an error reaches the top level (i.e. main), the simulation crashes.

Normally, there are print statements every time is thrown or caught. Users can choose from the `verbose` levels which determine how many of these statements are actually printed to the terminal. The value of verbose is set in the main-function (in `main.cpp`)

- 0: almost no error messages are printed. Only illegal inputs are reported (e.g. negative voltages). If you choose this level of verbosity, the code might crash without printing a message.
- 1: only *fatal* errors are printed. I.e. if we think this error is going to make the code crash, an error message is printed. If we think the error can be corrected by higher functions, no message is printed. It isn’t always known at the lowest level if an error is going to be fatal or not, so this is not exact. I.e. the code might crash without printing a message and you might see an error message, but the code can recover from this error. It is recommended you use this setting as standard.
- 2: all error messages are printed. Every time an error is thrown, a message is printed irrespectively of whether the error can be corrected or not. This setting might produce a lot of output in certain cases. If the code crashes and you don’t know why, it is recommended you change verbosity to 2 so you can see what the last error was (and this error will cause the code to crash).
- 3: Apart from printing all the error messages, also a message is printed every time a function in the Cycler or BasicCycler is started and terminated. This will allow you to find in which function the error happens. Note that this generates a lot of print statements, especially from the function `findCVcurrent` and `findCVcurrent_recursive` in `BasicCycler.cpp` (this function solves the nonlinear equation to find the current required to keep a given voltage).
- 4: On top of the previous output, also some information about the high-level flow of the program is printed by the functions in the Cycler, of the type “in function xxx we are going to charge the cell”. This might help to logically understand why an error occurs and where exactly it is (e.g. only after 1500 cycles into the degradation profile). As said before, there will be a lot of output from `findCVcurrent` and `findCVcurrent_recursive` so it might be a good idea commenting out the two print statements in that function (put a double back slash in front of the print statements at the start and end of `findCVcurrent` and `findCVcurrent_recursive` in `BasicCycler.cpp`).
- 5: On top of the previous output, detailed information is printed by the functions in the BasicCycler, of the type ‘in time step xx the cell voltage is yy V and the cell current is zz A’.
- 6: On top of the previous output, findCVcurrent_recursive prints details of the iterations when it is trying to find the current required to reach a given voltage.
- 7: On top of the previous output, a message is printed every time a function in the Cell is started. This will help locate the error if it happens in a function from Cell.

If you are calculating with multiple threads, error messages of different threads might interrupt each other. I.e. you get first some words from one message from thread x, then some words for an error message from thread y, etc. If you are debugging, it is recommended you isolate the thread which causes the problem by setting `isParallel{true}` in `constants.hpp` to `isParallel{false}`. Then you can clearly see the error messages from the one thread. This is especially the case if you set verbosity to 2 or higher.
All the errors also carry an identification code (i.e. a number). The excel document `Error IDs.xlsx` lists the numbers, the location where this error originates from and a very short description of the error. This information should already be in the error message which is printed but in case of confusion you can look the codes up.
To locate the error, it is recommended you slowly increase the value of verbose in the main function. Start with 2, if you can’t see the problem set it to 3, then to 4, etc. A lot of output will be generated from 3 onward (due to the print statements in `findCVcurrent_recursive`), too much to follow what is going on. Higher values will produce even more output, so they should only be used to exactly locate an error and understand why this is happening.
Due to the high output (from 3 or higher), the simulation will go significantly slower. Especially from 5 onward, the simulation slows down by orders of magnitude. So only set these higher levels when debugging.

## Common errors

Below is a list of some errors which might appear in the code if you change certain settings.

### Illegal state (error 15)

This error is thrown if an invalid battery state is detected (by the function `validState` in `state.cpp`. It is treated like a fatal error and the simulation stops. This is done because illegal states might lead to unexpected errors later in the code. E.g. if the temperature becomes too high, the numerical time integration might become unstable and you get infinite or NaN (see the last section of document [Characterisation parametrisation](../3_using/5_character_parametrisation.html) or the section on numerical problems below). This error will also be thrown if a cell has degraded too much (e.g. the electrodes are too thin). This check is done to avoid infinite loops, numerical problems etc. which will happen if the electrodes get very thin (or some other state has degraded too much).
However, an illegal state does not immediately lead to an error, and the boundaries of what is ‘legal’ are arbitrarily set, apart from the non-negativity of most variables (negative values will lead to errors).
If you get this error but want to continue simulating, you can increase the allowed range of the variables (e.g. increase the maximum temperature to 75 degrees if you want to simulate a cell of 60 degrees). You can also disable the throw-statement (i.e. add a double back slash in front of the statement `throw 15`) but this is not recommended.

### Illegal voltages or currents with the wrong signs, (errors 1004, 1005, 1006 and 1007)

These errors are the most difficult to solve. Users can define two ranges.

- The minimum and maximum allowed voltage of the cell (`Vmax` and `Vmin`), defined in the constructors of the Cell subclasses
- The upper and lower voltage between which you want to load the cell (`Vupp` and `Vlow`). E.g. doing a CC discharge to a Vlow > Vmin or following a drive cycle current profile while staying
between certain voltage limits. This variable is an input to the functions in the BasicCycler (and is therefore determined by the various functions in the Cycler).

Usually the problem is that in one way or another, the cell has gotten to a voltage which is just below the lower voltage limit, or just above. When the error message is printed, it gives the actual voltage and the limit which is exceeded. Often, these numbers might look identical but that is due to rounding. In reality, the actual voltage might be a few micro-Volts above the maximum or below the minimum (so if you print up to the mV range, the numbers look the same).
There is no easy solution for this problem. Simply increasing the allowed range typically doesn’t solve it because the limits are shared for all functions (e.g. you first charge to `Vupp`, then do something which causes the voltage to go just above `Vupp` (e.g. changing the cell temperature or current) and then you try to charge again: this last charge leads to an error in `CC_t_V` (error 1005) and increasing `Vupp` won’t solve it (because then you also charge a bit more). Some functions can recover from this error (e.g. `CC_t_CV_t` which does both a CC and CV will try to skip the CC phase and go straight to the CV phase which sometimes is possible but not always). Other functions will just throw on the error and the code crashes.

In this case, the only solution is to manually find how the cell got outside of the allowed range and then correcting that line. Isolate the thread which causes the error (i.e. comment out all others by adding a double backslash in front of the line defining and joining them) and increase the verbosity.


### File could not be opened (error 1001)

The code tries to create a folder with the user-defined name, consisting of the prefix (in main) and the identifier of the `Cycler` (in the functions in `cycling.cpp` and `degradation.cpp`). This means that the values of those strings must be valid names for folders. If they are not, the subfolder can’t be created and you get an error when you try to open a file in this folder.

Try to simplify the name and omit symbols and spaces (e.g. backslashes or forward slashes are not usually allowed, underscored should be fine).

### Arrays have the wrong length (error 107)

You have to increase the value of `len` in the struct `DEG_ID`, defined in `cell.hpp`.

### Numerical problems

Changing certain parameters (the radius, the diffusion constant, and all the thermal parameters (cooling, density, heat capacity, etc.) or using extreme states (too high temperature, too much degraded cells, etc.) might lead to errors in the code. There are two potential reasons

#### Numerical stability

See the section [Characterisation parametrisation](../3_using/5_character_parametrisation.md#note-on-numerical-stability-of-changing-the-radius-diffusion-constant-or-thermal-parameters).


#### Discretisation problems (e.g. error 101)

In this case, the problem is that because you take discrete time steps, the states vary in discrete steps too. This might mean that one state is valid, while the state in the next time step is illegal.

Consider a case where the cell is being discharged to the minimum voltage (e.g. 2.7V). Over the time steps, the lithium concentration in the anode increases (and it decreases in the cathode), which will cause the voltage to decrease. Ideally this goes in small steps (the lithium fractions change by a few percentages and the voltage by a few millivolts). In this case you can stop when you reached a step just below the limit (e.g. 2.69V).

But if the current or time step is very large, or the amount of active material or the diffusion constant is very small, the changes in surface concentration (and cell voltage) over one time step become very large. Suppose that the voltage at time t is 2.71V with a corresponding lithium fraction in the anode of 0.96, so the code takes another time step to further discharge the cell. But due to the large change in surface concentration, the lithium fraction in the anode now becomes 1.01, which is illegal (the lithium fraction has to be smaller than 1). This will lead to error 101.

To solve this, try the following.
- First check the parameters you are using (e.g. you got the wrong units and are now trying to run a 1000C current instead of a 1C). Idem for the diffusion constant, electrode thickness, etc.
- Secondly you can decrease the time step, which will decrease the change per time step.
- Another thing you can try is to decrease the rate of change in the cell’s current (`dIcell` in the constructors of the Cell-subclasses). This will solve the case where the problem happened when you were suddenly changing the applied current which lead to a crazy surface concentration (but not the case where you were just running a large current and then the error happened)
- Finally, you can be more conservative with your limits. In the case above, increasing the minimum voltage to 2.75V will solve the problem because you stop one time step earlier when the concentration is still in the allowed range.

## Other errors

For other problems, the error message hopefully explains what is wrong, which allows you to correct the error (e.g. you specified a negative current limit for the CV phase). To debug, try the following:
- Isolate the thread by setting `isParallel{true}` in `constants.hpp` to `isParallel{false}`.
- Increase the verbosity (in the functions in `cycling.cpp` and `degradation.cpp` which cause the problem) with the variable verbosity.
- Find where in the code the error happened, and why it happened.
- Try to solve it
- If you have no clue, try disabling the throw-statement. Certain errors are not really a problem, just things which should not happen or things which might cause problems later in the code. This is the case for the ‘illegal state’-error (15), the ‘initial state’-error (14) and a whole bunch of other errors. If the code runs fine after you disabled the error, then ignore the error. If you get all sorts or other errors, you will have to correct the original error.
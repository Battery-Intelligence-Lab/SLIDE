---
layout: default
title: Overall Approach
nav_order: 2
---

# Explanation of Battery-code

<p style='text-align: justify;'>
This section gives a brief overview of the code implemented in the pack version of `SLIDE`. It builds further upon the Module-project, which has a more in-depth explanation of certain parts of the code.
</p>

## Overall Approach

***Object Oriented Programming***


Currently, there is code implemented for a cell to do a few CC CV cycles, and for a cell to follow a predefined current profile. In both cases, the current, voltage and temperature of the cell is stored every few seconds in `.csv` files, which can be displayed using MATLAB scripts.

The code is written in the object-oriented programming framework. This means there are classes defining concepts (e.g., a cell must have a current, voltage and state of charge) and instances of these classes give their realisations for a simulation (e.g. a cell with a 1 A current, 3.7 V voltage and 50 % SoC). This is ideal for a large simulation since you can easily make more instances without having to change the code.

It also implies that every class is responsible for dealing with everything it involves. For instance, Cells calculate both their electrical or electrochemical model and their thermal model and they have a function to step forward in time which will update the voltage and temperature. Similarly, modules are responsible for keeping their cells consistent: ensuring they have the same current / voltage (depending whether it is series or parallel), ensuring thermal coupling between adjacent cells, etc.

So there is no separate `electrical`, `thermal` or `chemical` model which are calculated individually and then somehow coupled. All models are implemented in the same equations, and all are resolved together. This situation may change in future to make SLIDE more flexible. 

***Polymorphism***


The code relies very heavily on polymorphism. StorageUnit is an `Interface`, defining a class and its functions but without implementing anything. Other classes (Cell, Module, Battery) then implement all these functions. This means that all instances of these classes can be treated as StorageUnits. When a function from StorageUnit is called (e.g. getting the voltage), the correct implementation of that function is found in the most specific class of which this instance is a member (i.e. dynamic binding). So if other code works with StorageUnits, then it will work with cells, modules and batteries. In other words, for external code it doesn’t matter whether you work with a single cell or a battery of 500 kWh. However, specialised versions will play a role in future. 

The most complicated example of these are Modules. We typically think of modules as a stack of cells connected in series or parallel. But in large-scale batteries, these modules would then be connected in series or parallel to form a rack. And in turn, racks can be connected in series or parallel to give a container, etc.
Instead in the code, a Module is defined as a stack of StorageUnits connected in series or parallel. This means that you can make a Module consisting of a number of cells, but equally you can make a Module consisting of other Modules (which in turn might consist of Cells). The dynamic binding will ensure consistency (e.g. getting the voltage from the top level module will add up the voltages from the bottom modules [if series], and the voltages from the bottom modules are calculated as the sum of the voltages of all cells [if series], such that the total voltage is correctly calculated as the sum of the voltages of all cells.
Therefore, modules, racks, etc are all represented by the Module class.

Note that the polymorphism is achieved using smart pointers to `StorageUnit`, along with the keyword virtual in the function definition. 

***Class Diagram***

For more information about classes hieararchies and diagrams plase refer to our  [API documentation in Doxygen](../Doxygen/index.html)


***Time integration***


Every class inheriting from `StorageUnit (SU)` has a function to integrate in time, since every SU is responsible for its own time integration. Therefore, time integration is _cascaded down_, i.e. the `Cycler` will set a current, and then call the time integration function of the `StorageUnit` it involves, which will call the corresponding function of `Module`, which will call the time integration on each of its connected `StorageUnits` (either `Modules` which will again pass it on or `Cells` which can compute their own model). Therefore, time keeping is very strict to ensure consistency (all elements must take a time step of exactly the same duration).

There is one feature for adaptive time stepping: `Cycler` can instruct to take N time steps at the same time. This has two effects:

-	`Modules` instruct their connected `SUs` to take N time steps before returning to the Module. So e.g. in a parallel module, you will only correct the cell currents to equalise the voltages after N time steps. The time integration of the connection `SUs can be done on separate threads to accelerate the computation.
-	In SPM cells, the thermal model and degradation model is only resolved every N time steps. This is because temperature and degradation work on a smaller time scale, but note that N has to be small ( < 10) to ensure accuracy and stability, especially for the thermal model.


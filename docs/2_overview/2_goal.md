---
layout: default
title: Goal
nav_order: 2
---

# Goal of Slide

<p style='text-align: justify;'>
The main goal of Slide (Simulator for Lithium-Ion Degradation) is to allow fast simulation of degradation of li-ion batteries. Simulating 5000 full 1C cycles with a time resolution of 2 seconds (1C charge, 1C discharge) for one cell takes about 40 seconds. Including a CV phase on charge to a current limit of 0.05C increases the calculation time to 85 seconds.

The project uses the single particle model to simulate the main reaction in the battery. Various physics-based degradation models from literature have been implemented, and the user can select which degradation models to use. The parameters of the battery model, and of the different degradation models can be changed by the user.

Some degradation procedures have already been implemented, such as calendar ageing (resting a cell for a long time), cycle ageing (continuously charging and discharging) and drive cycle ageing (continuously repeating a given current profile). The settings of these basic procedures can be easily changed 
(e.g. the user can choose whether to do a CC or CCCV charge and the voltage window to be used).

The user can program his/her own degradation procedures similarly to the way a battery tester is programmed. 
The code has functions implemented to do a CC or CV (dis)charge (or a combination of both) where the parameters such as the current, or the set voltage can be determined by the user. By calling these functions with different parameters, the user can specify his/her own degradation experiments. 
Similarly, functions have been implemented to measure the (remaining) cell capacity, the OCV curves, do a pulse discharge test, etc., which the user can call do to regular check-ups during the degradation simulation.

The C++ code writes its results in *.csv files. The results of the check-ups during the degradation experiments (e.g. capacity measurements) are written in separate files. When a cell is being cycled, the user can choose to record the cell current, voltage and temperature at a constant time interval (e.g. 5 seconds) in which case this data is also stored in csv files. MATLAB functions to read these results have been implemented as well. E.g. to read the results of the pre-defined calendar ageing function, 
the user can run readCalendarAgeing.m, which will plot the outcomes. 
</p>

## Speed of calculation & data

The main advantage this code is that it fast and relatively flexible to use. The calculation speed depends on what you are simulating. Below are a few tips.

-  <p style='text-align: justify;'> Calculating a CV takes much longer than a CC. As reference, calculating 100 1C cycles with only CC takes about 0.9 second while calculating 100 1C cycles with also a 
    CV phase on charge (with a current limit of 0.05C) takes 1.9 second. This depends of course on how long the CV phase takes (a lower diffusion constant will give a longer 
    CV and therefore increase calculation time). </p>
- <p style='text-align: justify;'> The code can record periodic values of current, voltage and temperature. This slows down the code and creates a lot of data, especially when you simulate long term degradation. 
    Writing the data to files takes long, with the exact time depending on the speed of your hard disk. As a reference, the same 100 1C cycles (CC only) as before but with a 5s 
    time recording interval takes 18 seconds, with the vast majority of the extra time spent on writing 17 MB of data to the hard disk. Simulating the 100 cycles with an additional CV 
    charge and 5s data collection takes 22 seconds and 21 MB is written. When simulating degradation experiments, gigabytes of data will be generated and it can take up to hours to write all this data 
    (even though the actual calculations take only a couple of minutes). Especially for profile ageing (drive cycles), the amount of data generated (and time to write this data) is huge. 
    But it can be used to compare the simulation with lab data. </p>

- <p style='text-align: justify;'> Some degradation models significantly increase the calculation time. Most models only add about 10% to the calculation time but some models have more impact. E.g. The degradation models with the stress model from Dai (LAM model 1 or CS model 2) double the calculation time. </p>

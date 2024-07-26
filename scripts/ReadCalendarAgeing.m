% This script plots the results in the csv files written by the function
% CalendarAgeing in main.cpp
%
%
% Copyright (c) 2019, The Chancellor, Masters and Scholars of the University 
% of Oxford, VITO nv, and the 'Slide' Developers.
% See the licence file LICENSE for more information.

clc
close all
clear

%% User input: identification of the simulation
% In the code below, the user has to copy the input given in the
% main()-function in Main.cpp

% The value of the prefix (pref in c++)
pref = '0';                   

% The vlaue of the struct which defined which degradation models were used
% (deg in c++). Note that in Matlab we can directly input the arrays (e.g.
% if we used 2 models for one mechanism as is the case for LAM below)
sei_id = [2];                           % which SEI model(s) was (were) used
sei_por = 0;                            % whether the porosity was reduced (SEI_porosity from c++)
CS_id = [0];                            % which crach growth model(s) was (were) used
CS_diff = 0;                            % whether the diffusion constant was reduced (CS_diffusion from c++)
LAM_id = [2 3];                         % which LAM model(s) was (were) used
pl_id = [1];                            % which li-plating model(s) was (were) used


%% Identify which files should be read
% Each simulated experiment has its own identifier describing that
% experiment. This identifier is the last part of the name of the folder in
% which the results are written.
% These identifiers are defined in the variable name, which is redefined
% for every simulated experiment in CalendarAgeing in Degradation.cpp in
% the c++ code (note: the line says e.g. name = pref + "Cal-T5_SoC100"; In
% that case, the identifier is Cal-T5_SoC100 (while pref has the prefix and
% degradation model identification)
% In Matlab we can make a cell array with all the identifiers at once
% (rather than redefining it every time again)
IDs = {'Cal-T45_SoC50','Cal-T45_SoC90','Cal-T45_SoC100',...
    'Cal-T25_SoC50','Cal-T25_SoC90','Cal-T25_SoC100',...
    'Cal-T5_SoC50','Cal-T5_SoC90','Cal-T5_SoC100'};

% Degradation identifier
ageingID = printDEGID(sei_id, sei_por, CS_id, CS_diff, LAM_id, pl_id);
pathVar;

%% read and plot the battery states
% Makes 2 graphs:
%   one with the remainig capacity
%   one with the degradation details such as LLI
%       one subplot per degradation detail

FECx = false;                       % use 'time' as x-axis of the plots
                                    % if true, 'full equivalent cycles' are
                                    % used for the x-axis but since the
                                    % cell is resting, this doesn't make
                                    % sense
readAgeing_BatteryState

%% Read and plot the OCV curves
% Makes 1 graph with the OCV curves 
%   has a subplot per calendar regime and shows how the half-cell OCV
%   curves change due to degradation
readAgeing_OCV

%% Read and plot the voltage during the CCCV check-up cycles
% makes a large number of graphs, one per calendar regime
%   each figure has one subplot per cycle in the CCCV checkup
%   in the release of this code, it has a 6 subplots (0.5C, 1C and 2C
%   discharge and charge)
% Each plot shows how the voltage curve changes, due to changes in the OCV
% curve, resistance increase, and other degradation effects.
readAgeing_CCCV

%% Read and plot the voltage during the pulse charge
% makes one figure with the simulated voltage when a pulse profile is
% applied to the cell
%   the figure has one subplot per calendar regime
% Each plot shows how the voltage curve changes, due to changes in the OCV
% curve, resistance increase, and other degradation effects.
readAgeing_pulse

%% Read and plot the cycling data
% makes a large number of graphs, one per calendar regime
%   each figure has two subplots, one with the simulated current and
%   voltage, and a second one with the simulated temperature.
% It show the simulated cycling data during the degradation experiment. It
%   always shows the cycling data of the 'actual experiment' (i.e. the
%   calendar period). If the value of 'proc.includeCycleData' in
%   CalendarAgeing in Degradation.cpp was true, the cycling data will also
%   include the check-ups (i.e. the current/voltage/temperature while a
%   check-up is done). If the value of the variable was false, the cycling
%   data will only have the data from the experiments without the check-ups
readAgeing_cyclingData

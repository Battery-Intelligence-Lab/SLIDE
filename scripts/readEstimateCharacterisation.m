% Script to read the results of the characterisation fitting tool from the C++ code
% (estimateCharacterisationAtReferenceT, which calls all the other functions).
% It is ran by estimateCharacterisation() in Main.cpp
%
%
% Copyright (c) 2019, The Chancellor, Masters and Scholars of the University 
% of Oxford, VITO nv, and the 'Slide' Developers.
% See the licence file LICENCE.txt for more information.

clc
close all
clear

%% User input
% The user needs to tell us the names of the csv files with the OCV curves
% The variables below must have the same value as in the C++ code (in code
% block '1 USER INPUT' in estimateCharacterisation in determineCharacterisation.cpp

names = {'Characterisation_0.2C_CC_discharge.csv', ...      % names of the files with the measured voltage curves
    'Characterisation_0.5C_CC_discharge.csv',...
    'Characterisation_1C_CC_discharge.csv',... 
    'Characterisation_2C_CC_discharge.csv',...
    'Characterisation_3C_CC_discharge.csv'};
nameCCCVfit = 'characterisationFit_';                       % prefix appended before the name of the output data files with the simulations for the best fit

pathVar; % get necessary paths. 


%% Read and plot the best fit

% Find how to organise the subplots
nregimes = length(names);
if nregimes <= 15
    nrows = 3;                      % number of rows
else 
    nrows = 4;                      % number of columns
end
ncols = ceil(nregimes/nrows);


figure()

% Loop for all the cycles
for i=1:length(names)
    na1 = names{i};                                         % file with the measured data
    na2 = strcat(nameCCCVfit,names{i});                     % file with the simulated curves
    Vcell = csvread( fullfile(pathvar.results_folder, na1) );
    Vsim = csvread( fullfile(pathvar.results_folder, na2) ); 
    
    subplot(nrows,ncols,i)
    plot(Vcell(:,1),Vcell(:,2),'.-',Vsim(:,1), Vsim(:,2))
    legend('measured','simulated')
    xlabel('discharged Ah')
    ylabel('voltage [V]')
    title(['Fitted voltage curve for ' names{i}])
end
% Script to read the results of the OCV fitting tool from the C++ code
% (determineOCV.cpp, in which the function estimateOCVparameters calls
% everything). It is ran by estimateOCVparameters() in Main.cpp
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
% block '1 USER INPUT' in estimateOCVparameters in determineOCV.cpp

namecell = 'OCVfit_cell.csv';       % the user-supplied measured OCV curve of the cell (supplied by the user)
nameneg = 'OCVfit_anode.csv';       % the user-supplied anode OCV curve
namepos = 'OCVfit_cathode.csv';     % the user-supplied cathode OCV curve
nameOCV = 'OCVfit_sim.csv';         % the simulated OCV curve (simulated by the code)

pathVar;

%% plot the best fit

OCVcell = csvread( fullfile(pathvar.data_folder, namecell) );   % the user-supplied measured OCV curve of the cell (supplied by the user)
OCVn = csvread( fullfile(pathvar.data_folder, nameneg) );       % the user-supplied anode OCV curve
OCVp = csvread( fullfile(pathvar.data_folder, namepos) );       % the user-supplied cathode OCV curve
OCVsim = csvread( fullfile(pathvar.results_folder, nameOCV) );     % the simulated OCV curve (simulated by the code)
    % discharged charge [Ah]
    % cell OCV [V]
    % anode li-fraction [-]
    % anode voltage [V]
    % cathode li-fraction [-]
    % cathode voltage [V]

% Compare the simulated and measured OCV curve of the cell. And the second
% subplot shows the half-cell OCV curves
figure()
subplot(2,1,1)
    plot(OCVcell(:,1),OCVcell(:,2),'.')
    hold on
    plot(OCVsim(:,1), OCVsim(:,2),OCVsim(:,1),OCVsim(:,4),OCVsim(:,1), OCVsim(:,6))
    legend('measured OCV','simulated OCV', 'simulated OCV_{anode}','simulated OCV_{cathode}')
    xlabel('discharged Ah')
    ylabel('OCV [V]')
    title('Fitted OCV curve')
subplot(2,1,2)
    plot(OCVn(:,1),OCVn(:,2),'.-b')
    hold on
    plot(OCVsim(:,3), OCVsim(:,4),'-b','linewidth',4)
    plot(OCVp(:,1),OCVp(:,2),'.-r')
    plot(OCVsim(:,5), OCVsim(:,6),'-r','linewidth',4)
    legend('anode full', 'anode used','cathode full', 'cathode used')
    xlabel('lithium fraction [-]')
    ylabel('OCV [V]')
    title('half-cell curves indicating the region used in the full cell')
    % The dots show the full half-cell curves (li-fractions 0 to 1)
    % The lines show which segment of the full curves is used (li-fractions
    %   which actually occur in the cell)



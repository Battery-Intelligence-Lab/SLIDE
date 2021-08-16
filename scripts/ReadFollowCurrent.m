% script to read the results of FollowCurrent, where a cell has to follow a
% certain current load profile
%
%
% Copyright (c) 2019, The Chancellor, Masters and Scholars of the University 
% of Oxford, VITO nv, and the 'Slide' Developers.
% See the licence file LICENCE.txt for more information.

clc
clear
close all

%% User input: identification of the simulation
% In the code below, the user has to copy the input given in the
% main()-function in Main.cpp

% The value of the prefix (pref in c++)
pref = '0';        
pathVar; % get necessary paths. 


% The vlaue of the struct which defined which degradation models were used
% (deg in c++). Note that in Matlab we can directly input the arrays (e.g.
% if we used 2 models for one mechanism as is the case for LAM below)
sei_id = [2];                           % which SEI model(s) was (were) used
sei_por = 0;                            % whether the porosity was reduced (SEI_porosity from c++)
CS_id = [0];                            % which crach growth model(s) was (were) used
CS_diff = 0;                            % whether the diffusion constant was reduced (CS_diffusion from c++)
LAM_id = [2 3];                         % which LAM model(s) was (were) used
pl_id = [1];                            % which li-plating model(s) was (were) used

% Specify the name of the current profile which was followed (defined in
% the function FollowCurrent in Cycling.cpp
orig = csvread(fullfile(pathvar.data_folder, 'Current Profile drive cycle UDDS.csv'));

%% Read the simulation results

% Identifiers for the folders and file
ID = 'FollowCurrent';                   % the identification string in the folder name (as specified in the function CCCV)
fileName = 'CyclingData';               % We want to read the cycling data of the results
ageingID = printDEGID(sei_id, sei_por, CS_id, CS_diff, LAM_id, pl_id); % Degradation identifier
fol = strcat(pref,'_',ageingID,'_',ID); % Folder in which to find the results
fol = fullfile(pathvar.results_folder, fol);

% A loop to reach each data batch from the simulation
fi = dir (fol);                         % list the files in the subfolder
A = [];
for jj = 1:length(fi)                   % loop through all files
    ni = fi(jj).name;    
    if contains(ni,fileName)            % check if it is a file with cyclingData, i.e. if the name contains 'cyclingData'
        % older versions of Matlab don't have the 'contains' command. In
        % that case, use strfind. Replace 'if contains(ni,fileName)' with
        % 'if strfind(ni,fileName)'

        % Read the file
        na = fullfile(fol,ni);          % Full name of the file (including the path)
        try
            Ai = csvread(na);
        catch
            Ai = nan(1,15);             % If we can't read the file, print a warning to the user, and fill the data with NaNs
            warning(['warning no file ' na ' could be found'])
        end

        % Add up the cumulative variables (e.g. total time, total charge
        % throughput, etc.). In the C++ code, they are reset to 0 after a
        % data batch is written, so here we need to add them up to insert
        % the new data batch 'after' the previous one.
        % These cumulative variables are in columns 1-3,9-15
        cumulative = [1:3 9:15];
        if ~isempty(A)
            Ai(:,cumulative) = Ai(:,cumulative) + A(end,cumulative);
        end

        % Append the new data behind the old
        A = [A ; Ai];
    end
end

% Make a struct with the data in named fields
sim.timetot     = A(:,1);               % total time since the start [s]
sim.Ahtot       = A(:,2);               % total charge throughput
sim.Whtot       = A(:,3);               % total energy throughput
sim.I           = A(:,4);               % current, < 0 for charge [A]
sim.V           = A(:,5);               % cell voltage [V]
sim.OCVp        = A(:,6);               % cathode potential [V]
sim.OCVn        = A(:,7);               % anode potential [V]
sim.T           = A(:,8);               % cell temperature [K]
sim.timeCha     = A(:,9);               % time spent on charging [s]
sim.AhCha       = A(:,10);              % charged charge [Ah]
sim.WhCha       = A(:,11);              % charged energy [Wh]
sim.timeDis     = A(:,12);              % time spent on discharging [s]
sim.AhDis       = A(:,13);              % discharged charge [Ah]
sim.WhDis       = A(:,14);              % discharged energy [Wh]
sim.timeRest    = A(:,15);              % time spent on rest [s]

%% read the original current profile

% Make a structure with the current profile the battery was supposed to
% follow such that we can compare the simulation with the profile.
prof.I      = orig(:,1);           % current [A]
prof.time   = cumsum(orig(:,2));   % cumulative time since the start [s]

%% Make a figure showing the cycling data

% First compare the current profile with the current the battery actually
% delivered. If a voltage limit was reached, there will be a discrepancy
% between both.
% The second and third subplot show the simulated voltage and temperature
figure()
subplot(3,1,1)
    plot(sim.timetot,sim.I,prof.time,prof.I,'.')
    legend('real current','current from the profile')
    ylabel('[A]')
    xlabel('[sec]')
    title('current')
subplot(3,1,2)
    plot(sim.timetot,sim.V)
    ylabel('[V]')
    xlabel('[sec]')
    title('cell voltage')
subplot(3,1,3)
    plot(sim.timetot,sim.T-273)
    ylabel('[degrees]')
    xlabel('[sec]')
    title('temperature')
    
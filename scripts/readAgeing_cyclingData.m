% Script to read and plot the cycling data from the degradation simulation.
% The cycling data contains the current, voltage and temperature of the
% cell at a fixed time interval.
% If 'proc.includeCycleData' in the degradation function in Degradation.cpp
% is true, the cycling data contains both the periods from the actual
% degradation experiment, and from the check-ups.
% If 'proc.includeCycleData' in the degradation function in Degradation.cpp
% is false, the cycling data contains only the periods from the actual
% degradation experiment.
%
% This script should not be executed on its own, but is called by one of
% three higher-level scripts: 
%   readCalendarAgeing
%   readCycleAgeing
%   readProfileAgeing
%
%
% Copyright (c) 2019, The Chancellor, Masters and Scholars of the University 
% of Oxford, VITO nv, and the 'Slide' Developers.
% See the licence file LICENCE.txt for more information.

fileName = 'CyclingData';                           % Name of the file which contains the cycling data

%% Read the csv files & plot the cycling data
% Note we read and plot in one loop because this avoids you have to store
% all the cyling data from all simulations first (and then plot all of it). 
% The cycling data can be very large (many GigaBytes) so we want to 
% minimise the memory requirements.

% A loop to read & plot the files of each ageing regime
for i=1:length(IDs)
    
    fol = strcat(pref,'_',ageingID,'_',IDs{i});     % Folder in which to find the file
    fol = fullfile(pathvar.results_folder, fol);
    fi = dir (fol);                                 % list the files in the subfolder
    A = [];
    
    % Check that the folder exist and has data files
    if length(fi)<1
        warning(['warning no files exist in folder ' fol])
        A = nan(1,15);
    end
    
    % A loop to read each cycling data batch from this ageing regime
    for jj = 1:length(fi)                           % loop through all files from this folder
        ni = fi(jj).name;    
        if contains(ni,fileName)                    % check if it is a file with cyclingData, i.e. if the name contains 'cyclingData'
            % older versions of Matlab don't have the 'contains' command. In
            % that case, use strfind. Replace 'if contains(ni,fileName)' with
            % 'if strfind(ni,fileName)'
           
            % Read the file
            na = fullfile(fol,ni);                  % Full name of the file
            try
                Ai = csvread(na);
            catch
                Ai = nan(1,15);                     % If we can't read the file, print a warning to the user, and fill the data with NaNs
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
            
            % Append the new data behind the old data
            A = [A ; Ai];
        end        
    end
    
    % Store the data in a struct
    capini = state{i}.cap(1);          % initial capacity of this cell [Ah] (from readAgeing_BatteryState)
    cycdat.timetot     = A(:,1);       % total time since the start [s]
    cycdat.Ahtot       = A(:,2);       % total charge throughput
    cycdat.Whtot       = A(:,3);       % total energy throughput
    cycdat.I           = A(:,4);       % current, < 0 for charge [A]
    cycdat.V           = A(:,5);       % cell voltage [V]
    cycdat.OCVp        = A(:,6);       % cathode potential [V]
    cycdat.OCVn        = A(:,7);       % anode potential [V]
    cycdat.T           = A(:,8);       % cell temperature [K]
    cycdat.timeCha     = A(:,9);       % time spent on charging [s]
    cycdat.AhCha       = A(:,10);      % charged charge [Ah]
    cycdat.WhCha       = A(:,11);      % charged energy [Wh]
    cycdat.timeDis     = A(:,12);      % time spent on discharging [s]
    cycdat.AhDis       = A(:,13);      % discharged charge [Ah]
    cycdat.WhDis       = A(:,14);      % discharged energy [Wh]
    cycdat.timeRest    = A(:,15);      % time spent on rest [s]
    cycdat.FEC         = cycdat.Ahtot / (2*capini); % cumulative full equivalent cycles [-]
    
    clear('A', 'Ai');                  % clear the matrix to free up memory
    
    % Plot the current, voltage and temperature
    figure()
    subplot(2,1,1)
        yyaxis left
            if FECx
                plot(cycdat.FEC,cycdat.I)
                xlabel('full equivalent cycles [-]');
            else
                plot(cycdat.timetot/3600,cycdat.I)
                xlabel('time [hour]')
            end
            ylabel('[A]')
        yyaxis right
            if FECx
                plot(cycdat.FEC,cycdat.V)
                xlabel('full equivalent cycles [-]');
            else
                plot(cycdat.timetot/3600,cycdat.V)
                xlabel('time [hour]')
            end
            ylabel('[V]')
        legend('current','voltage')
        title(strcat('cell current and voltage of', " ", IDs{i}),'Interpreter', 'none') 
    subplot(2,1,2)
        if FECx
            plot(cycdat.FEC,cycdat.T-273)
            xlabel('full equivalent cycles [-]');
        else
            plot(cycdat.timetot/3600,cycdat.T-273)
            xlabel('time [hour]')
        end
        ylabel('[degrees]')
        title(strcat('cell temperature of', " ", IDs{i}), 'Interpreter', 'none')    
        
        clear('cycdat');               % clear the struct with the data to free up memory
end









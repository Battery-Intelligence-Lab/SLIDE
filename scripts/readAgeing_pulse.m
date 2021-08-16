% Script to read and plot the voltage response from the current pulse
% profiles applied during the check-ups in the degradation simulations.
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

fileName = 'DegradationData_CheckupPulse_';     % Name of the file which contains the pulse data

%% Read the pulse curves

nCheckups_max = 100;                            % maximum number of check-ups done 
nCheck = zeros(size(IDs));                      % number of checkups for each ageing regime

% Make a cell array pulse with the data
% the data for degradation regime i is in pulse{i,:} (i.e. one row)
%   The columns contain the data for the consecutive check-ups (i.e. column
%   j has data from check-up j)
% Each cell {i,j} has a structure with the data with the data from check-up 
%   j from degradation regime i). E.g. the field 'I' has the current
% i.e. pulse{i,j}.I gives the current during the pulse test of the
%   j'th check-up of the i'th degradation regime

% A loop to read the files of each cycling regime
for i=1:length(IDs)
    
    nCheck(i) = 0;                              % number of checkups done for this ageing regime
    fol = strcat(pref,'_',ageingID,'_',IDs{i}); % Folder in which to find the file
    nrErrors = 0;                               % number of files not found for this regime
    
    % A loop to read all the check-ups from one cycling regime
    for j=1:nCheckups_max
        
        name = fullfile(pathvar.results_folder, fol, strcat(fileName,num2str(j-1),'.csv'));    % Full name of the file
        
        try
            % Read the file
            A = csvread(name);
            nCheck(i) = j;
        catch
            % No file of this check-up exists, potentially because not all 
            % check-ups were done for this ageing regime
            nrErrors = nrErrors + 1;            % increase the number of files not found
            
            % Store all NaNs
            A = nan(1,15);
            
            % If this was the first check-up, no data was available for
            % this experiment. Print a warning to the user
            if j == 1
                warning(['warning no pulse data for ageing regime ' IDs{i} ' could be found'])
            end
        end
        
        % Store the data
        pulse{i,j}.timetot     = A(:,1);       % total time since the start [s]
        pulse{i,j}.Ahtot       = A(:,2);       % total charge throughput
        pulse{i,j}.Whtot       = A(:,3);       % total energy throughput
        pulse{i,j}.I           = A(:,4);       % current, < 0 for charge [A]
        pulse{i,j}.V           = A(:,5);       % cell voltage [V]
        pulse{i,j}.OCVp        = A(:,6);       % cathode potential [V]
        pulse{i,j}.OCVn        = A(:,7);       % anode potential [V]
        pulse{i,j}.T           = A(:,8);       % cell temperature [K]
        pulse{i,j}.timeCha     = A(:,9);       % time spent on charging [s]
        pulse{i,j}.AhCha       = A(:,10);      % charged charge [Ah]
        pulse{i,j}.WhCha       = A(:,11);      % charged energy [Wh]
        pulse{i,j}.timeDis     = A(:,12);      % time spent on discharging [s]
        pulse{i,j}.AhDis       = A(:,13);      % discharged charge [Ah]
        pulse{i,j}.WhDis       = A(:,14);      % discharged energy [Wh]
        pulse{i,j}.timeRest    = A(:,15);      % time spent on rest [s]
            
        % If we have not found 2 files, assume no further data is
        % available (we search for data of nCheckups_max check-ups but
        % most cells won't have done so many)
        if nrErrors > 2
            break;
        end
        
    end % end the loop for the different check-ups of degradation regime i
end % end the loop for all the degradation regimes


%% Plot the voltages of the pulse curves
% Make one figure per cycling regime
%   On that figure, make one subplot per cycle
%   On that subplot, show the data from all the checkups

figure()

% Find how to organise the subplots
nregimes = length(IDs);
if nregimes <= 15
    nrows = 3;                      % number of rows
else 
    nrows = 4;                      % number of columns
end
ncols = ceil(nregimes/nrows);

for i=1:length(IDs)
    col = winter(nCheck(i));
    subplot(nrows,ncols,i)

    for j=1:nCheck(i)
        plot(pulse{i,j}.Ahtot, pulse{i,j}.V,'color',col(j,:))
        hold on
    end
    xlabel('[Ah]')
    ylabel('[V]')
    tit = strcat(IDs{i});
    title(tit,'Interpreter', 'none')
 
    % try to make a legend giving the full equivalent cycles of when
    % the check-up was done. this will only work if you have previously
    % ran 'readCycleAgeing_BatteryState.m'
    try
        legi = cell(nCheck(i),1);
        if FECx
            xax = state{i}.FEC;
            for t = 1:nCheck(i)
                legi{t} = [num2str(round(xax(t))) ' FEC'];
            end
        else
            xax = state{i}.time;
            for t = 1:nCheck(i)
                legi{t} = [num2str(round(xax(t)/24)) ' days'];
            end
        end
        legend(legi)
    catch
    end
end





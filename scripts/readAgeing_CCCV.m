% Script to read and plot the voltage response from the CC or CCCV cycles
% applied during the check-ups in the degradation simulations.
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

fileName = 'DegradationData_CheckupCycle_';     % Name of the file which contains the CCCV cycles

%% Read the CCCV curves
nCheckups_max = 100;                            % maximum number of check-ups done  
Crates = [0.5 1 2];                             % C rates of the CC phases, MUST BE THE SAME AS THE CYCLES ACTUALLY DONE for the check-up (defined in the struct proc in Degradation.cpp)
nCheck = zeros(size(IDs));                      % number of checkups for each ageing regime

% Make a cell array CCCV with the data
% the data for degradation regime i is in CCCV{i,:} (i.e. one row)
%   The columns contain the data for the consecutive check-ups (i.e. column
%   j has data from check-up j)
% Each cell {i,j} has a structure with the data with the data from check-up 
%   j from degradation regime i). E.g. the field 'I' has the current
%   the data from the different cycles in the CCCV check-up is separated in
%   the different columns.
% i.e. CCCV{i,j}.I(:,k) gives the current during the k'th cycle of the
%   j'th check-up of the i'th degradation regime

% A loop to read the files of each degradation regime
for i=1:length(IDs)
    
    nCheck(i) = 0;                                  % number of checkups done for this ageing regime
    fol = strcat(pref,'_',ageingID,'_',IDs{i});     % Folder in which to find the file
    nrErrors = 0;                                   % number of files not found for this regime
    
    % A loop to read all the check-ups from this one degradation regime
    for j=1:nCheckups_max
        name = fullfile(pathvar.results_folder, fol,strcat(fileName,num2str(j-1),'.csv'));    % Full name of the file
        try
            % Read the file
            A = csvread(name);
            
            % separate between the different cycles
            % We do this by finding when the sign of the current changes
            % (i.e. when we go from charging to discharging or vise versa)
            I = A(:,4);                             % current, < 0 for charge [A]
            f = find( sign(I(1:end-1)) ~= sign(I(2:end)) & I(1:end-1) ~= 0); % find the rows where the sign of the current changes but where the current is not 0
            s = [1 ; f+1];                          % array with the row numbers where the data from one cycle starts
            f = [f ; length(I)];                    % array with the row numberss where the data from one cycle finishes
            maxlength = max(f-s)+1;                 % maximum length of the data from all cycles
            
            % Make the matrices for this check-up
            CCCV{i,j}.timetot     = nan(maxlength,length(s));       % total time since the start [s]
            CCCV{i,j}.Ahtot       = nan(maxlength,length(s));       % total charge throughput
            CCCV{i,j}.Whtot       = nan(maxlength,length(s));       % total energy throughput
            CCCV{i,j}.I           = nan(maxlength,length(s));       % current, < 0 for charge [A]
            CCCV{i,j}.V           = nan(maxlength,length(s));       % cell voltage [V]
            CCCV{i,j}.OCVp        = nan(maxlength,length(s));       % cathode potential [V]
            CCCV{i,j}.OCVn        = nan(maxlength,length(s));       % anode potential [V]
            CCCV{i,j}.T           = nan(maxlength,length(s));       % cell temperature [K]
            CCCV{i,j}.timeCha     = nan(maxlength,length(s));       % time spent on charging [s]
            CCCV{i,j}.AhCha       = nan(maxlength,length(s));       % charged charge [Ah]
            CCCV{i,j}.WhCha       = nan(maxlength,length(s));       % charged energy [Wh]
            CCCV{i,j}.timeDis     = nan(maxlength,length(s));       % time spent on discharging [s]
            CCCV{i,j}.AhDis       = nan(maxlength,length(s));       % discharged charge [Ah]
            CCCV{i,j}.WhDis       = nan(maxlength,length(s));       % discharged energy [Wh]
            CCCV{i,j}.timeRest    = nan(maxlength,length(s));       % time spent on rest [s]
            nCheck(i) = 0;                                          % number of checkups done for this ageing regime
        
            % Store the data from each cycle (k) in this check-up (j) from
            % this degradation regime (i)
            for k = 1:length(s)
                r =  s(k):f(k) ;                             % rows of this cycle
                l = length(r);                               % length of this cycle

                CCCV{i,j}.timetot(1:l,k)     = A(r,1);       % total time since the start [s]
                CCCV{i,j}.Ahtot(1:l,k)       = A(r,2);       % total charge throughput
                CCCV{i,j}.Whtot(1:l,k)       = A(r,3);       % total energy throughput
                CCCV{i,j}.I(1:l,k)           = A(r,4);       % current, < 0 for charge [A]
                CCCV{i,j}.V(1:l,k)           = A(r,5);       % cell voltage [V]
                CCCV{i,j}.OCVp(1:l,k)        = A(r,6);       % cathode potential [V]
                CCCV{i,j}.OCVn(1:l,k)        = A(r,7);       % anode potential [V]
                CCCV{i,j}.T(1:l,k)           = A(r,8);       % cell temperature [K]
                CCCV{i,j}.timeCha(1:l,k)     = A(r,9);       % time spent on charging [s]
                CCCV{i,j}.AhCha(1:l,k)       = A(r,10);      % charged charge [Ah]
                CCCV{i,j}.WhCha(1:l,k)       = A(r,11);      % charged energy [Wh]
                CCCV{i,j}.timeDis(1:l,k)     = A(r,12);      % time spent on discharging [s]
                CCCV{i,j}.AhDis(1:l,k)       = A(r,13);      % discharged charge [Ah]
                CCCV{i,j}.WhDis(1:l,k)       = A(r,14);      % discharged energy [Wh]
                CCCV{i,j}.timeRest(1:l,k)    = A(r,15);      % time spent on rest [s]
                nCheck(i) = j;
            end % end loop for the different cycles of one check-up
            
        catch
            % No file of this check-up exists, potentially because not all 
            % check-ups were done for this ageing regime
            nrErrors = nrErrors + 1;                        % increase the number of files not found
            
            % Store all NaNs
            maxlength = 1;
            ncyc = length(Crates)*2;
            % Make the matrices for this check-up
            CCCV{i,j}.timetot     = nan(maxlength,ncyc);       % total time since the start [s]
            CCCV{i,j}.Ahtot       = nan(maxlength,ncyc);       % total charge throughput
            CCCV{i,j}.Whtot       = nan(maxlength,ncyc);       % total energy throughput
            CCCV{i,j}.I           = nan(maxlength,ncyc);       % current, < 0 for charge [A]
            CCCV{i,j}.V           = nan(maxlength,ncyc);       % cell voltage [V]
            CCCV{i,j}.OCVp        = nan(maxlength,ncyc);       % cathode potential [V]
            CCCV{i,j}.OCVn        = nan(maxlength,ncyc);       % anode potential [V]
            CCCV{i,j}.T           = nan(maxlength,ncyc);       % cell temperature [K]
            CCCV{i,j}.timeCha     = nan(maxlength,ncyc);       % time spent on charging [s]
            CCCV{i,j}.AhCha       = nan(maxlength,ncyc);       % charged charge [Ah]
            CCCV{i,j}.WhCha       = nan(maxlength,ncyc);       % charged energy [Wh]
            CCCV{i,j}.timeDis     = nan(maxlength,ncyc);       % time spent on discharging [s]
            CCCV{i,j}.AhDis       = nan(maxlength,ncyc);       % discharged charge [Ah]
            CCCV{i,j}.WhDis       = nan(maxlength,ncyc);       % discharged energy [Wh]
            CCCV{i,j}.timeRest    = nan(maxlength,ncyc);       % time spent on rest [s]
            
            % If this was the first check-up, no data was available for
            % this experiment. Print a warning to the user
            if j == 1
                warning(['warning no CCCV data for ageing regime ' IDs{i} ' could be found'])
            end
            
            % If we have not found 2 files, assume no further data is
            % available (we search for data of nCheckups_max check-ups but
            % most cells won't have done so many)
            if nrErrors > 2
                break;
            end
        end % end try-catch block to read the data from check-up j of ageing regime i
        
    end % end loop for the different check-ups of ageing regime i
end % end the loop for all the different ageing regimes

% the cumulative variables (such as charge throughput) keep running over
% the different cycles of one check-up. Reset them at the start of every
% cycle, such that they give the cumulative values within this one cycle
for i=1:length(IDs)
    for j=1:nCheck(i)
        for k = 1:(length(Crates)*2)
            CCCV{i,j}.timetot(:,k) = CCCV{i,j}.timetot(:,k) - CCCV{i,j}.timetot(1,k);
            CCCV{i,j}.Ahtot(:,k)   = CCCV{i,j}.Ahtot(:,k) - CCCV{i,j}.Ahtot(1,k);
            CCCV{i,j}.Whtot(:,k)   = CCCV{i,j}.Whtot(:,k) - CCCV{i,j}.Whtot(1,k);
        end
    end
end

%% Plot the voltages of the CCCV curves
% Make one figure per ageing regime
%   On that figure, make one subplot per cycle
%   On that subplot, show the data from all the checkups

for i=1:length(IDs)
    col = winter(nCheck(i));
    figure()
    
    % Loop for all the cycles
    for k = 1:(length(Crates)*2)
        cr = num2str(Crates( ceil(k/2)));       % C rate used for this CCCV cycle
        subplot(length(Crates),2,k)
        
        % loop for all the check-ups
        for j=1:nCheck(i)
            plot(CCCV{i,j}.Ahtot(:,k), CCCV{i,j}.V(:,k),'color',col(j,:))
            hold on
        end
        xlabel('[Ah]')
        ylabel('[V]')
        tit = strcat(IDs{i}, ', cycle at', " ", cr,'C');
        title(tit,'Interpreter', 'none')
    end
        
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





% Script to read and plot the battery states from the check-ups in the 
% degradation simulations. The battery states indicate the remaining 
% capacity, as well as the underlying reasons for the lost capacity (such 
% as the mount of lithium lost)
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

fileName = 'DegradationData_batteryState.csv';      % Name of the file which contains the battery states

%% Read the battery states
% A loop to read the files of each ageing regime
for i=1:length(IDs)
    fol = strcat(pref,'_',ageingID,'_',IDs{i});     % Folder in which to find the file
    na = fullfile(pathvar.results_folder, fol,fileName);                    % Full name of the file
    try
        A = csvread(na);
    catch
        A = nan(50,50);                             % If we can't read the file, print a warning to the user, and fill the data with NaNs
        warning(['warning no file ' na ' could be found'])
    end

    % Store the results in a struct called state
    na = IDs{i};                                    % name
    state{i}.text = na;                             % name of the file
    na(strfind(na,'_')) = ' ';                      % name with space instead of _
    state{i}.name = strcat(ageingID,'-',na);
    state{i}.cap = A(:,5);                          % remaining capacity [Ah]
    capini = state{i}.cap(1);                       % initial capacity [Ah]
    state{i}.cycleNr = A(:,1);                      % cycle number [-]
    state{i}.time = A(:,2);                         % cumulative time since the start [hour]
    state{i}.Ah = A(:,3);                           % cumulative charge throughput since the start [Ah]
    state{i}.FEC = A(:,3) / (2*capini);             % cumulative full equivalent cycles [-]
    state{i}.Wh = A(:,4);                           % cumulative energy throughput since the start [Wh]
    state{i}.cap_rel = state{i}.cap/capini*100;     % relative remaining capacity [%]
    state{i}.T = A(:,16)-273;                       % Cell temperature [degrees]
    state{i}.delta = A(:,17);                       % Thickness of the SEI layer [m]
    state{i}.LLI = A(:,18);                         % Lost lithium [As]
    state{i}.thickp = A(:,19);                      % thickness of the positive electrode [m]
    state{i}.thickn = A(:,20);                      % thickness of the negative electrode [m]
    state{i}.ep = A(:,21);                          % volume fraction of the positive electrode [-]
    state{i}.en = A(:,22);                          % volume fraction of the negative electrode [-]
    state{i}.ap = A(:,23);                          % effective surface area of the positive electrode [m2]
    state{i}.an = A(:,24);                          % effective surface area of the negative electrode [m2]
    state{i}.CS = A(:,25);                          % total area of the surface cracks [m2]
    state{i}.Dp = A(:,26);                          % diffusion constant of the positive electrode at reference temperature [m s-1]
    state{i}.Dn = A(:,27);                          % diffusion constant of the negative electrode at reference temperature [m s-1]
    state{i}.r = A(:,28);                           % resistance of the electrodes [Ohm m2]
    state{i}.deltapl = A(:,29);                     % thickness of the plated lithium layer [m]
    state{i}.R = A(:,30);                           % total resistance of the cell [Ohm]
    state{i}.An = A(:,31);                          % active surface area of the anode (excluding crack surface), i.e. the product of the effective surface area with the electrode volume
end

%% Plot the capacities
col = jet(length(IDs));                             % Use a different colour for each degradation regime
    
figure()
for i=1:length(IDs)
    c = col(i,:);                                   % colour for this regime
    if FECx
        plot(state{i}.FEC,state{i}.cap_rel,'color',c); % use 'full equivalent cyles' as x-axis
    else
        plot(state{i}.time,state{i}.cap_rel,'color',c); % use 'time' as x-axis
    end
    hold on
end

% Make legend, labels & title
legend(IDs,'Interpreter','none');
if FECx
    xlabel('full equivalent cycles [-]');
else
    xlabel('time [hour]')
end
ylabel('capacity [%]');
title('Remaining capacity')

%% Plot the degradation details
% Make a figure with the subplots for every state of the battery related to
% battery degradation
figure()
subplot(4,4,1)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);  
        y = state{i}.delta;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('thickness of the SEI layer')
    ylabel('m')
    ylim([0.9*m 1.1*M]);
subplot(4,4,2)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);   
        y = state{i}.LLI/3600;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c); % /3600 to go from As to Ah
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('lost lithium inventory')
    ylabel('[Ah]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,3)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);    
        y = state{i}.thickp;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('cathode thickness')
    ylabel('[m]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,4)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);    
        y = state{i}.thickn;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('anode thickness')
    ylabel('[m]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,5)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);    
        y = state{i}.ep;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('cathodic volume fraction') % of active material
    ylabel('[-]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,6)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);   
        y = state{i}.en;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('anodic volume fraction') % of active material
    ylabel('[-]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,7)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);    
        y = state{i}.ap;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('cathodic effective surface area')
    ylabel('[m2/m3]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,8)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);    
        y = state{i}.an;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('anodic effective surface area')
    ylabel('[m2/m3]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,9)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);   
        y = state{i}.An;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('anodic active surface area excluding crack surface')
    ylabel('[m2]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,10)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);       
        y = state{i}.CS;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('area of surface cracks at the anode surface')
    ylabel('[m2]')
    ylim([0.9*m 1.1*M]); 
    % Note: the value of this surface can be compared with the surface area
    % of the anode (subplot 9). E.g. if CS/An = 0.1 then the total crack
    % surface area is 10% of the original surface area (as it would be in
    % the absence of cracks)
subplot(4,4,11)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);   
        y = state{i}.Dp;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('cathode diffusion constant')
    ylabel('[m/s]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,12)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);           
        y = state{i}.Dn;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('anode diffusion constant')
    ylabel('[m/s]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,13)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);         
        y = state{i}.R; 
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('total DC resistance')
    ylabel('[Ohm]')
    ylim([0.9*m 1.1*M]);
subplot(4,4,14)
    m = 10^10;                          % minimum value
    M = 0;                              % maximum value
    for i=1:length(IDs)
        c = col(i,:);           
        y = state{i}.deltapl;
        if FECx
            x = state{i}.FEC;
        else
            x = state{i}.time;
        end
        plot(x,y,'color',c);
        hold on
        m = min(m, min(y));
        M = max(M, max(y));
    end
    title('thickness of the plated li layer')
    ylabel('[m]')
    ylim([0.9*m (1.1*M+10^(-9))]);
    


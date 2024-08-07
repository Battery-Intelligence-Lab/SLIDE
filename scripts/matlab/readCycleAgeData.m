% Script to read the file written by the check-up in Procedure.cpp.
% It is one file with the properties of all cells over time.
% It also plots the statistics of the cooling systems

clc
close all
clear

cell_nominalCapacity = 16;          % the nominal capacity of cells. 

%% Read the files with the cell data

% name = 'proctest_p_module_0.000000_checkUp.csv';
pp = {'cool1_capSpread_RSpread_degSpread_contactR_balance_EPFL',...
    'cool2_capSpread_RSpread_degSpread_contactR_balance_EPFL',...
    'cool3_capSpread_RSpread_degSpread_contactR_balance_EPFL',...
    'cool4_capSpread_RSpread_degSpread_contactR_balance_EPFL'};

for kkk = 1:length(pp)
    pre = pp{kkk}
    close all


name = strcat(pre,'_checkUp.csv');

% Read the cell checkups
try
    A = csvread(name,1,0);                         
    % skip the first row, which has the cell IDs since Matlab cannot read
    % strings
catch
    disp(['warning, could not read file ' name])
end
% column 1 = ID of the row
%   1       parameters of the cell-to-cell variation for this cell [cap; resistance; degradation ] 
%   2       total battery charge throughput (same for every cell)
%   3       throughput of each cell [time ; Ah ; Wh]
%   4       cell capacity in Ah
%   9       separator row to separate different data streams
%   10 and above are state variables (10+i is state[i])
%
% Note that the rows with 1 are only written once at the start of the
% procedure. Rows 2,3,4,10+ are written for every check-up. 

% Delete the last column, which is all 0 (the loops in c++ wrote a comma
% after every value, so the last column is empty)
A = A(:,1:end-1);

% Read the cell statistics
try
    C = csvread(strcat(pre,'_cellStats.csv'));       
catch
    disp(['warning, could not read file ' name])
end
%   5       histogram for I (101 values)
%   6       histogram for V (101 values)
%   7       histogram for T (101 values)

% Read the edges of the histogram
B = csvread('histogram_edges.csv');
    % I V T from minimum to maximum

%% Process the cell results
% Per property, make a matrix with column = cell, row = checkup
ID = A(:,1);
[~,ncel] = size(A);
ncel = ncel - 1;            % for the column with the IDs

through.totalAh = A(ID == 2,2:end);
t = A(ID == 3,2:end);
through.time = t(1:3:end,:);
through.Ah = t(2:3:end,:);
through.Wh = t(3:3:end,:);
through.FEC = through.Ah / (2*cell_nominalCapacity);

capacity = A(ID == 4,2:end);
for i=1:ncel
    leg{i} = num2str(i);
end

% states if desired: 
%   20  delta 		the thickness of the SEI layer [m]
%   21  LLI 		the lost lithium [As]
%   22  thickp 		the thickness of the cathode [m]
%   23  thickn		the thickness of the anode [m]
%   24  ep 			the volume fraction of active material in the cathode [-]
%   25  en 			the volume fraction of active material in the anode [-]
%   26  ap 			the effective surface area of the cathode [m2 m-3]
%   27  an			the effective surface area of the anode [m2 m-3]
%   28  CS			the surface area of the cracks at the surface of the negative particle [m2]
%   29  Dp			the diffusion constant at reference temperature of the cathode [m s-1]
%   30  Dn			the diffusion constant at reference temperature of the anode [m s-1]
%   31  rp 			the specific resistance of the cathode [Ohm m2]
%   32  rn 			the specific resistance of the anode
%   33  rcc 		the specific resistance of the separator
%   34  delta_pl 	the thickness of the plated lithium layer [m]
%   35  SoC			the state of charge of the cell
%   36  T 			the cell temperature [K]
%   37  I 			the current of the cell
states.delta = A(ID == 20,2:end);
states.LLI = A(ID == 21,2:end);
states.CS = A(ID == 28,2:end); 
states.delta_pl = A(ID == 34,2:end);
% amount of active material = thickness * effective surface area * elec_surf
    %   elec_surf does not change, so just track relative change of thickness*a
    amp = A(ID == 22,2:end) .* A(ID == 26,2:end);
    amn = A(ID == 23,2:end) .* A(ID == 27,2:end);
states.AMp = amp ./ amp(1,:) * 100;     % [%] of active material left
states.AMn = amn ./ amn(1,:) * 100;     % [%] of active material left
states.thickn =  A(ID == 23,2:end);
states.en =  A(ID == 25,2:end);

% total resistance = Rsei + Rp + Rn + Rcc
    %   = rsei*delta / (thick_n*a_n*elec_surf) + rp / (thick_p*a_p*elec_surf) 
    %       + rn / (thick_n*a_n*elec_surf) + rcc / elec_surf
    % relative change (ignore constant elec_surf)
    rsei = 2037.4; % in constructor of cell
    Rsei = states.delta .* rsei ./ (A(ID == 23,2:end) .* A(ID == 27,2:end));
    Rp = A(ID == 31,2:end) ./ (A(ID == 22,2:end) .* A(ID == 26,2:end));
    Rn = A(ID == 32,2:end) ./ (A(ID == 23,2:end) .* A(ID == 27,2:end));
    Rcc = A(ID == 33,2:end);
    Rtot = Rsei + Rp + Rn + Rcc;
states.Rtot = Rtot ./ Rtot(1,:) * 100;  % relative resistance [%]

ID = C(:,1);
hist.I = C(ID == 5,2:end);
hist.V = C(ID == 6,2:end);
hist.T = C(ID == 7,2:end);
hist.Iedge = B(:,1);
hist.Vedge = B(:,2);
hist.Tedge = B(:,3);

clear('A','C','t','B');

%% Plot the cell results
cols = lines(ncel);

% Plot cell capacities
f1 = figure();
for i=1:ncel
    plot(through.FEC(:,i),capacity(:,i),'.-','color',cols(i,:))
    hold on
end
xlabel('FEC')
ylabel('capacity [Ah]')
grid on
% legend(leg)


% Plot individual states
f2 = figure();
subplot(4,2,1)
    for i=1:ncel
        plot(through.FEC(:,i),states.delta(:,i),'.-','color',cols(i,:))
        hold on
    end
    ylabel('delta')
    grid on
subplot(4,2,2)
    for i=1:ncel
        plot(through.FEC(:,i),states.LLI(:,i),'.-','color',cols(i,:))
        hold on
    end
    ylabel('LLI')
    grid on
subplot(4,2,3)
    for i=1:ncel
        plot(through.FEC(:,i),states.CS(:,i),'.-','color',cols(i,:))
        hold on
    end
    ylabel('crack surface')
    grid on
subplot(4,2,4)
    for i=1:ncel
        plot(through.FEC(:,i),states.delta_pl(:,i),'.-','color',cols(i,:))
        hold on
    end
    ylabel('delta plating')
    grid on
subplot(4,2,5)
    for i=1:ncel
        plot(through.FEC(:,i),states.AMp(:,i),'.-','color',cols(i,:))
        hold on
    end
    ylabel('positive active material [%]')
    grid on
subplot(4,2,6)
    for i=1:ncel
        plot(through.FEC(:,i),states.thickn(:,i),'.-','color',cols(i,:))
%         plot(through.FEC(:,i),states.AMn(:,i),'.-','color',cols(i,:))
        hold on
    end
    ylabel('negative electrode thickness []')
%     ylabel('negative active material [%]')
    grid on
subplot(4,2,7)
    for i=1:ncel
        plot(through.FEC(:,i),states.Rtot(:,i),'.-','color',cols(i,:))
        hold on
    end
    ylabel('total resistance [%]')
    grid on
subplot(4,2,8)
    for i=1:ncel
        plot(through.FEC(:,i),states.en(:,i),'.-','color',cols(i,:))
        hold on
    end
    ylabel('negative volume fraction [-]')
    grid on
    
% 
% % Plot histogram for one random cell
% i = 100;
% figure()
%     subplot(3,1,1)
%     plotHistogram(hist.Iedge,hist.I(:,i))
% %     legend(IDs)
%     xlabel('current')
%     ylabel('frequency')
% %     legend(leg)
% subplot(3,1,2)
%     plotHistogram(hist.Vedge,hist.V(:,i))
% %     legend(IDs)
%     xlabel('voltage')
%     ylabel('frequency')
% subplot(3,1,3)
%     plotHistogram(hist.Tedge-273,hist.T(:,i))
% %     legend(IDs)
%     xlabel('temperature')
%     ylabel('frequency')

% Plot histograms of cell statistics
lin = true; % plot a line per cell. if false, it plots a bar graph
f3 = figure();
subplot(3,1,1)
    plotHistogramMultiple(hist.Iedge,hist.I, lin)
%     legend(IDs)
    xlabel('current')
    ylabel('frequency')
%     legend(leg)
subplot(3,1,2)
    plotHistogramMultiple(hist.Vedge,hist.V, lin)
%     legend(IDs)
    xlabel('voltage')
    ylabel('frequency')
subplot(3,1,3)
    plotHistogramMultiple(hist.Tedge-273,hist.T, lin)
%     legend(IDs)
    xlabel('temperature')
    ylabel('frequency')

%% Read and plot the statistics of the cooling systems

name = strcat(pre,'_checkModules.csv');

try
    A = csvread(name,1,0);                         
    % skip the first row, which has the module IDs since Matlab cannot read strings
catch
    disp(['warning, could not read file ' name])
end
% column 1 = ID of the row
%   1       number of cells connected to this module / coolsystem
%   5       histogram for T (101 values)
%   6       histogram for Q/Ncells (101 values)
%   7       histogram for flowrate/Ncells (101 values)
%   8       histogram for operating power of the fans / Ncells
%   9       histogram for cooling power of the AC system / Ncells (0 for
%               all but top level module)
%   10      histogram for operating power of the AC system / Ncells (0 for
%               all but top level module)
%   99       separator row to separate different data streams

% Delete the last column, which is all 0
A = A(:,1:end-1);
ID = A(:,1);

% Check whether we simulate Battery or a Module
if A(1,2) == A(1,3)           % battery has the same number of cells as the top level module
    bat = true;               % store separately
else
    bat = false;
end

% Sort per number of cells connected
%   cell array, one cell per 'level' of the battery (same number of cells connected)
%   2nd level is also a cell array with one value per cell
%   3rd level is a struct with the statistics of that cell
nc = A(1,2:end);                % number of cells of each module
nci = unique(nc);               % unique number of cells
if bat
    nci = [nci nci(end)];       % additional level for battery with same number of cells
end
Thist = cell(length(nci),1);    % cell array with one cell per level
for i=1:length(nci)       
    
    % Logical index indicating which module is at this level
    if ~ bat || i < length(nci)-1   % no battery or not the top level
        ind = nc == nci(i);         % modules with this number of cells
    elseif i == length(nci)-1       % top level module
        ind = nc == nci(i);
        ind(1) = false;             % remove the Battery
    else                            % Battery level
        ind = nc == nci(i);
        ind(2) = false;             % remove the top level module
    end
    
    % Combine into a struct
    ncs = sum(ind);             % number of modules at this level
    Thist{i} = cell(ncs,1);
    k = 1;
    for j=find(ind)             % loop for all modules at this level
        Thist{i}{k}.Ncells = nci(i);   % number of cells connected to this level -> edges should be multiplied by this value
        Thist{i}{k}.T = A(ID == 5,j+1); %+1 for 1st column with ID numbers
        Thist{i}{k}.Q = A(ID == 6,j+1);
        Thist{i}{k}.flr = A(ID == 7,j+1); 
        Thist{i}{k}.E = A(ID == 8,j+1); 
        Thist{i}{k}.Qac = A(ID == 9,j+1); 
        Thist{i}{k}.Eac = A(ID == 10,j+1); 
        k = k+1;
    end
end

% Read the edges of the histogram
B = csvread('histogram_thermal_edges.csv');
Tedge.T = B(:,1);
Tedge.Q = B(:,2);
Tedge.flr = B(:,3);
Tedge.E = B(:,4);
Tedge.Qac = B(:,5);
Tedge.Eac = B(:,6);

% plot
lin = true; % plot a line per cell. if false, it plots a bar graph
k = 1;
f4 = figure();
for i=1:length(nci)
    
    % Combine all statistics of this level to one matrix
    Tbins = nan(length(Thist{i}{1}.T),length(Thist{i}));
    Qbins = nan(length(Thist{i}{1}.T),length(Thist{i}));
    flrbins = nan(length(Thist{i}{1}.T),length(Thist{i}));
    Ebins = nan(length(Thist{i}{1}.T),length(Thist{i}));
    QACbins = nan(length(Thist{i}{1}.T),length(Thist{i}));
    EACbins = nan(length(Thist{i}{1}.T),length(Thist{i}));
    for j = 1:length(Thist{i})
        Tbins(:,j) = Thist{i}{j}.T;
        Qbins(:,j) = Thist{i}{j}.Q;
        flrbins(:,j) = Thist{i}{j}.flr;
        Ebins(:,j) = Thist{i}{j}.E;
        QACbins(:,j) = Thist{i}{j}.Qac;
        EACbins(:,j) = Thist{i}{j}.Eac;
    end
    
    % get edges of this level (for total module, ie not per cell)
    T_edge = Tedge.T - 273;  % T in degrees
    Q_edge = Tedge.Q * Thist{i}{1}.Ncells; % total cooling power [W]
    flr_edge = Tedge.flr * Thist{i}{1}.Ncells; % total flow rate [W]
    E_edge = Tedge.E * Thist{i}{1}.Ncells; % total operating power [W]
    Qac_edge = Tedge.Qac * Thist{i}{1}.Ncells; % cooling power of the AC system [W]
    Eac_edge = Tedge.Eac * Thist{i}{1}.Ncells; % operating power of the AC system [W]
    
    % Calculate the mean operating power
    dx = abs(E_edge(2) - E_edge(1));
    xmin = E_edge(1) - 10*dx;
    xmax = E_edge(end)+10*dx;
    x = [xmin ; E_edge ; xmax];
    xmid = x(1:end-1) + 1/2*diff(x);
    for j = 1:length(Thist{i})
        Emean(k) = sum(xmid .* Ebins(:,j) / sum(Ebins(:,j)));
        k = k+1;
    end    
    
    % Plot this level
    subplot(length(nci),4,(i-1)*4+1)
        plotHistogramMultiple(T_edge,Tbins, lin)
        xlabel('temperature')
        ylabel('frequency')
        title(strcat('T of modules with-',num2str(Thist{i}{1}.Ncells),'-cells'))
    subplot(length(nci),4,(i-1)*4+2)
        plotHistogramMultiple(Q_edge,Qbins, lin)
        % If there are HVAC bins, draw them separately
        if sum(sum(QACbins)) > 0
            hold on
            plotHistogramMultiple(Qac_edge, QACbins, lin)
            legend('Q from children','Q of AC system')
        end
        xlabel('cooling power [W]')
        ylabel('frequency')
        title(strcat('heat evacuated in modules with-',num2str(Thist{i}{1}.Ncells),'-cells'))
    subplot(length(nci),4,(i-1)*4+3)
        plotHistogramMultiple(flr_edge,flrbins, lin)
        xlabel('flow rate [m3 s-1]')
        ylabel('frequency')
        title(strcat('flow rate of modules with-',num2str(Thist{i}{1}.Ncells),'-cells'))
    subplot(length(nci),4,(i-1)*4+4)
        plotHistogramMultiple(E_edge,Ebins, lin)            % E bin is only operating power of the fan
        % If there are HVAC bins, draw them separately
        if sum(sum(EACbins)) > 0
            hold on
            plotHistogramMultiple(Eac_edge, EACbins, lin) 
            legend('E for fan','E for AC system')
        end
        xlabel('operating power [W]')
        ylabel('frequency')
        title(strcat('operating power of modules with-',num2str(Thist{i}{1}.Ncells),'-cells'))
end
% 
% % Note on how to calculate the operating power for thermal management system
% %   power to operate fan is stored for every module in histogram E
% %   power to operate HVAC system is calculated below in 3 ways (last one is most correct) 
% %       COP * heat exctraced from 2nd level modules using histogram
% %       COP * heat exctraced from 2nd level modules using C++ variable
% %       COP * heat extracted by the HVAC system from the top level module
% %           using histogram
% %   Additionally, the throughput file contains the energy required during
% %       every half cycle. This is the total energy (operating fans +
% %       operating HVAC where HVAC is using 3rd method)
% % NOTE2: statistics are based on data snapshots every ndata seconds
% % (currently 100) while the stored total values (i.e. 2nd method and
% % throughput) is continuously updated. 
% % Therefore, the most accurate one is what is stored in throughput.
% % Second most accurate is 2.2 (COP * heat exctraced from 2nd level modules
%     % using C++ variable) since it uses the total C++ variable, but it uses the
%     % heat evacuated from racks, not the battery
% % 2.1 and 2.3 use statistics so are rough (Especially for on/off control)
% 
% % Calculate the power for the HVAC system
% COP = 3; % value from Schimpe's paper
% Q_edge = Tedge.Q * Thist{end}{1}.Ncells;
% dx = abs(Q_edge(2) - Q_edge(1));
% xmin = Q_edge(1) - 10*dx;
% xmax = Q_edge(end)+10*dx;
% x = [xmin ; Q_edge ; xmax];
% xmid = x(1:end-1) + 1/2*diff(x);
% Qcool = sum(xmid .* Thist{end}{1}.Q / sum(Thist{end}{1}.Q));
% Phvac = COP * Qcool;
% 
% % Extract the value from the written file
% tottime = A(3,2)/3600;      % 2nd column = top level module, [hour]
% Qevac = A(4,2)/3600;   % heat evaculated from children of top level module, [Wh]
% Qabs = A(5,2)/3600;    % heat absorbed in the top level module [Wh]
% Phvac2 = COP * Qevac / tottime;
% 
% % Read the heat exctraced by the HVAC system from the top level module
% edge = Tedge.Qac*Thist{end}{1}.Ncells;
% bins = Thist{end}{1}.Qac;
% dx = abs(edge(2) - edge(1));
% xmin = edge(1) - 10*dx;
% xmax = edge(end)+10*dx;
% x = [xmin ; edge ; xmax];
% xmid = x(1:end-1) + 1/2*diff(x);
% Qhvac = sum(xmid .* bins / sum(bins));
% Phvac3 = COP * Qhvac;
% 
% % Display the mean power to operate all cooling systems
% disp([' the mean total power to operate all fans in the rack is ' num2str(sum(Emean/1000)) ' kW'])
% disp(['the mean power for the HVAC system is ' num2str(Phvac/1000) ' kW or ' num2str(Phvac2/1000) ' kW or ' num2str(Phvac3/1000) ' kW']) % slight difference dependin on how it is calculated (here or in C++)
% disp([' so total ancillary load for thermal management system is ' num2str(Phvac2/1000+sum(Emean/1000)) ' kW'])
% % the throughput contains the energy for every full cycle. Given that a
% % cycle takes about 2h, the numbers should be about twice the figure given
% % here


%% analyse the efficiency and charge throughput

% Read the file
name = strcat(pre,'_throughput.csv');
try
    A = csvread(name,1,0);                         
    % skip the first row, which has the column headings since Matlab cannot read strings                      
catch
    disp(['warning, could not read file ' name])
end
% column 1 = ID of the row
%   1       CC charge (<0)
%   2       CV charge (< 0)
%   3       CC discharge ( >0)
%   4       CV discharge (>0)
% column 2  charge throuhgput of all the cells
% column 3  energy throughput of all the cells
% column 4  energy requipred to run the cooling system during this half cycle [Wh]
% column 5  energy losses in the power electronic converter of the battery [Wh] 

ID = A(:,1);
Ah = A(:,2);
Wh = A(:,3); % energy throughput of the battery compartment (including losses in contact R, cells and useful stored energy)
Ecool = A(:,4);
convloss = A(:,5);
x = 1:length(ID);
rebalance = 0:10*2:length(x);     % cycles we rebalance
ymax = max(abs(Wh))*1.1;
ymin = min(abs(Wh))*0.9;
    % *4 cause 4 actions in every cycle
    % note: sometimes an action is skipped because an error happened
    % so the x-axis goes slightly wrong in that case

f5 = figure();
yyaxis left
for i=1:4
    plot(x(ID == i), abs(Wh(ID == i)),'.')
    hold on
end
yyaxis right
for i=1:4
    plot(x(ID == i), Ecool(ID == i) + convloss(ID == i),'.')
    hold on    
end
l = legend('CC charge','CC discharge','ancillary losses');
l.AutoUpdate = 'Off';
grid on
xlabel ('action = cycle/2')
ylabel('Wh')
yyaxis left
hold on
% for i=1:length(rebalance)
%     plot([rebalance(i) rebalance(i)],[ymin ymax],'-r')
%     hold on
% end

% NOTE: energy for cooling as fraction of charge throughput should be about
% twice as the printed number of the mean HVAC power as fraction of 1C power
%   here: energy for cooling vs throughput of one cycle
%   there: mean power to operate thermal management system (vs power of
%       total battery at 1C, about 130 kW)
%   and it is double cause the ancillary load is both during charge and
%       discharge

% Calculate energy efficiency per cycle
cha = abs(Wh(ID == 1));% + Wh (ID == 2));
dis = abs(Wh(ID == 3));% + Wh (ID == 4));
N = min(length(cha),length(dis));                   % in case you finished mid-way a cycle
cha = cha(1:N);
dis = dis(1:N);
cells = cha - dis;                                  % losses in the cells
Ecool = Ecool(1:2*N);
convloss = convloss(1:2*N);
ID = ID(1:2*N);
cool = (Ecool(ID == 1) + Ecool(ID == 3));           % operating power of the coolsystem during charge and discharge
conv = (convloss(ID == 1) + convloss(ID == 3));     % losses in the converter during charge and discharge                             
loss = cells + cool + conv;

% Efficiency = out / in
%   (Wh_dis - conv_dis - therm_dis) / 
%       (Wh_cha + conv_cha + therm_cha)
%   since Wh is energy to/from the battery compartment
idis = ID == 3;
icha = ID == 1;
IN = abs(Wh(icha)) + abs(convloss(icha)) + abs(Ecool(icha));
OUT = abs(Wh(idis)) - abs(convloss(idis)) - abs(Ecool(idis)); 
eta2 = OUT ./ IN * 100;

% Relative to total input power, which is energy charged + converter loses
% + ancillary power during charge
chatot = abs(Wh(ID == 1)) + convloss(ID == 1) + Ecool(ID == 1);
chatot = chatot(1:N);

% L = [cool conv cells];
L = [cool./chatot conv./chatot cells./chatot]  * 100;

eta = 100-loss ./ chatot * 100;    
f6 = figure();
subplot(2,1,1)
    plot(eta(2:end),'.')
    ylabel('[%]')
    xlabel('cycle')
    title('efficiency')
%     ylim([60 110])
hold on
plot(eta2(2:end),'.')
    grid on
subplot(2,1,2)
    bar(L(2:end,:),'stacked')            % stacked bar graph with all losses
    ylabel('[%]')
    xlabel('cycle')
    legend('cool system power','converter losses','cell losses')
    title('breakdown of the losses')
    grid on
%     ylim ([0 3*10^5])
%     ylim([-5 30])
% NOTE FOR simulations with no thermal model: except for single cell, the
% thermal system was still cooling the losses of the PE converter, hence
% the operating losses for the thermal management system are not 0.
% They don't affect anything else so simply hard-code that they should be 0
warning('if no thermal model, set coollosses to 0 if needed')

% Note that losses during charge are significantly smaller than during
% discharge
    % During charge, the voltage of the battery is higher. This causes
    % fewer losses in the converter (smaller V ratio to convert) while at
    % the same time increasing the 'base' by which they are divided to get
    % the %. So you get a smaller amount of losses divided by a larger
    % amount of power
    

%% Save figures
saveas(f1,strcat(pre,'_cap.png'))
saveas(f2,strcat(pre,'_cellstates.png'))
saveas(f3,strcat(pre,'_cellstats.png'))
saveas(f4,strcat(pre,'_thermstats.png'))
saveas(f5,strcat(pre,'_usable.png'))
saveas(f6,strcat(pre,'_efficiency.png'))


end




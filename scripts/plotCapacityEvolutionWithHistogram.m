close all
clear
clc

COOLS = [1];

    
%% read Cells
cell_nominalCapacity = 16;          % the nominal capacity of cells. 

% Struct with cell arrays with values over time for plotting
capcell.x = cell(5,1);              % x values for mean and std
capcell.mean_rel = cell(5,1);
capcell.std_rel = cell(5,1);
capcell.min_rel = cell(5,1);
capcell.max_rel = cell(5,1);
capcel.x2 = cell(5,1);              % x values for usable and efficiency
capcell.usable_rel = cell(5,1);
capcell.eta = cell(5,1);
Thist.cells = cell(5,1);

for kk = COOLS
    
    % Read the file
    pp = num2str(kk);
    pre = strcat('cool',pp,'_capSpread_RSpread_degSpread_contactR_balance_EPFL');
    name = strcat(pre,'_checkUp.csv');

    try
        A = csvread(name,1,0);  
    catch
        A = NaN(50,50);
        disp(['warning, could not read file ' name])
    end
    % Delete the last column
    A = A(:,1:end-1);
    B = csvread('histogram_edges.csv');

    % Process the file
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

    hist.I = A(ID == 5,2:end);
    hist.V = A(ID == 6,2:end);
    hist.T = A(ID == 7,2:end);
    hist.Iedge = B(:,1);
    hist.Vedge = B(:,2);
    hist.Tedge = B(:,3);
    Thist.cells{kk} = hist.T;
    Thist.edges =  hist.Tedge;

    % Store capacity values over time
    capcell.mean_rel{kk} = mean(capacity,2)/cell_nominalCapacity*100; % mean of all cells
    capcell.std_rel{kk} = std(capacity,0,2)/cell_nominalCapacity*100; % mean of all cells
    capcell.min_rel{kk} = min(capacity,[],2)/cell_nominalCapacity*100; % mean of all cells
    capcell.max_rel{kk} = max(capacity,[],2)/cell_nominalCapacity*100; % mean of all cells
    capcell.x{kk} = mean(through.FEC,2);% ignore difference in utilisation between cells
    capcell.capend{kk} = capacity(end,:);

    % Usable (1C charge) capacity and 1C efficiency
    % Read the file
    name = strcat(pre,'_throughput.csv');
    try
        A = csvread(name,1,0);                                  
    catch
        disp(['warning, could not read file ' name])
    end
    ID = A(:,1);
    Ah = A(:,2);
    Wh = A(:,3);    
    Ecool = A(:,4);
    convloss = A(:,5);
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
    L = [cool./cha conv./cha cells./cha]  * 100;
    eta1 = 100-loss ./ cha * 100;  
    cap1c = -cha/1000;
    

    % set values after balancing to NaN
    bal = eta1 > 93.5;
    cap1c(bal) = NaN;
    eta1(bal) = NaN;

    % Values over time
    capcell.usable_rel{kk} = cap1c/130*100;
    capcell.eta{kk} = eta1;
    capcell.x2{kk} = 1:length(eta1);
end


%% Try to plot capacity with distribution of end capacity


for kk = COOLS

% Calculate the histogram
capend = capcell.capend{kk};
figure()
h = histogram(capend);
x = h.BinEdges / cell_nominalCapacity*100;
y = h.BinCounts;
y = y / sum(y) * 100; % frequency

% plot with swapped axis
% figure()
% for i=1:length(y)
%     xx = [x(i) x(i+1) x(i+1) x(i) x(i)];
%     yy = [0 0 y(i) y(i) 0];
%     fill(yy,xx,'b')
%     hold on
% end

f = 20; % fontsize
f2 = 15;%font size for histogram
leg = {'on','on/off local','on/off global','proportional local','proportional global'};
col = lines(5);
fig = figure();
 % plot line with mean
        % plot line with mean
        for i=kk
            plot(capcell.x{i},capcell.mean_rel{i},'color',col(i,:),'linewidth',3)
            hold on
        end
        % shade area mean +- std
        for i=kk
            ymax = capcell.mean_rel{i} + capcell.std_rel{i};
            ymin = capcell.mean_rel{i} - capcell.std_rel{i};
            ystd = [ymax ; flipud(ymin)];
            xstd = [capcell.x{i} ; flipud(capcell.x{i})];
            fill(xstd,ystd,col(i,:),'facealpha',0.2,'linestyle','none')
            hold on
        end
        % plot min and max
        for i=kk
            plot(capcell.x{i},capcell.min_rel{i},':','color',col(i,:),'linewidth',0.5)
            plot(capcell.x{i},capcell.max_rel{i},':','color',col(i,:),'linewidth',0.5)
            hold on
        end
        % shade area between min and max
        for i=kk
            ymax = capcell.max_rel{i};
            ymin = capcell.min_rel{i};
            ystd = [ymax ; flipud(ymin)];
            xstd = [capcell.x{i} ; flipud(capcell.x{i})];
            fill(xstd,ystd,col(i,:),'facealpha',0.1,'linestyle','none')
            hold on
        end
%         title(strcat('cell capacity for ',{' '},leg{kk}))
        xlabel('FEC')
        ylabel('[%]')
        grid on
        xlim([0 15000])
        ylim([0 110])
        set(gca,'FontSize',f)

% Add new axis on this plot
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes [left bottom width height]
xl = ax1.XLim;          
    % min and max -> xlim(1) is at ax1_pos(1)
    % xlim(2) is at ax1_pos(1) + ax1_pos(3)
xend = max(capcell.x{i});
    % last x-value on the graph
    % location is ax1_pos(1) + ax1_pos(3)*xend/xlim(1)
ax2_x1 = ax1_pos(1) + ax1_pos(3)*xend/xl(end) + 0.01;
    % left point where 2nd axis should start [add 1% space]
ax2_x2 = ax1_pos(1) + ax1_pos(3) - ax2_x1;
    % width of 2nd axis if you want to fill all the way to the end
ax2_x2 = ax1_pos(3) / 5; % width to just fill 1/5 of the original plot
ax2_pos = [ax2_x1  ax1_pos(2) ax2_x2 ax1_pos(4)];
ax2 = axes('Position',ax2_pos,...
    'XAxisLocation','top',...
    'YaxisLocation','right',...
    'Color','g'); 
for i=1:length(y)
    xx = [x(i) x(i+1) x(i+1) x(i) x(i)];
    yy = [0 0 y(i) y(i) 0];
    fill(yy,xx,col(kk,:))
    hold on
end
ylim([0 110])
xlim([0 20])
grid on
xlabel('[%] of cells')
ax2.YTick = [];% remove this if you want to check the axes align correctly
set(ax2,'YAxisLocation','right','XAxisLocation','top');%,'FontSize',f2)


    saveas(fig,strcat('A system ',num2str(kk),' with histogram.png'))
    
end





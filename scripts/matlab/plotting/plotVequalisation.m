%
% plotting Vequalisation results for 1 mOhm contact resistance.
% paperCode.cpp
%
% Created on: 04 Jan 2023
%  Author(s): Jorn Reniers, Volkan Kumtepeli
%

clear all; close all; clc;

% Load the data:

Ncell = 5; % Number of cells simulated.

results = cell(Ncell,1);

folder = "../results";

for i=1:5
    file_name = "paper_Vequalisation_0.001000_pmod_cell"+ num2str(i)+"_cellData.csv";
    path = fullfile(folder, file_name);
    results{i} = readmatrix(path,'Range',7);
end

%% Plot currents:
legends = split(sprintf('Cell-%d\n',1:Ncell));
legends = legends(1:Ncell);
figure;

for i=1:5
    result = results{i};
    plot(result(:,5), result(:,1));
    hold on;
end

legend(legends(:));
grid on;
title('Current [A]');

%% Plot voltage:
figure;

for i=1:5
    result = results{i};
    plot(result(:,5), result(:,2));
    hold on;
end

legend(legends(:));
grid on;
title('Voltage [V]');

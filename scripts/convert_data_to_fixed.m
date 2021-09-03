% This file takes non-constantly sampled data to make fixed-sampled data
% for faster interpolation! 

% NOT IMPLEMENTED YET.

data_folder = fullfile('..', 'data');
file_name = 'Kokam_OCV_C.csv'; %Kokam_OCV_NMC
data_path = fullfile(data_folder, file_name);

data = readmatrix(data_path);


%%

dts = diff(data(:,1));

min_dt = min(dts)/2;
n = (data(end,1) - data(1,1)) / min_dt;

new_x = linspace(data(1,1), data(end,1), length(data(:,1))); %ceil(n)
new_y = interp1(data(:,1), data(:,2), new_x);


close all;

plot(data(:,1), data(:,2),'*-'); hold on;
plot(new_x, new_y, 'r');

%%

new_file_name = [file_name(1:end-4), '_fixed', file_name(end-3:end)];
new_data_path = fullfile(data_folder, new_file_name);

writematrix([new_x(:), new_y(:)], new_data_path);
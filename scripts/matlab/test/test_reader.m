clear variables; close all; clc;
tic;
tab = slide_to_table('mydata.slide');
toc
tic;
tabcsv = readtable('mydata.csv');
toc
error = norm(tab{:,:} - tabcsv{:,:}, 1);
fprintf('Error is: %4.4f\n', error);
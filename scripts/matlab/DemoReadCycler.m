% Test file for demo.
clear variables; close all; clc; 
filename = 'Cycler1_Cell_SPM_cellData.csv';

states = readmatrix(fullfile('../results',filename),"NumHeaderLines",3);

I = states(:,1);
V = states(:,2);
SOC = 100*states(:,3); 
T = states(:,4) - 273.15; 
t = states(:,5)/3600; 
Ah = states(:,6);
Wh = states(:,7);

figure; 

subplot(3,1,1);
plot(t,V); grid on;
xlabel('Time (h)');
ylabel('Voltage (V)');


subplot(3,1,2);
plot(t,I); grid on;
xlabel('Time (h)');
ylabel('Current (A)');


subplot(3,1,3);
plot(t,SOC); grid on;
xlabel('Time (h)');
ylabel('SOC (%)');





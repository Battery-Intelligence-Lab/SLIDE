% This file is for simulating ECM cases
clear variables; close all; clc;

Tminutes = 5;

ECM_single_default_pulse = sim('Cell_ECM_single_default_pulse.slx');

ECM_single_default_pulse.I = squeeze(ECM_single_default_pulse.I);
ECM_single_default_pulse.SOC = squeeze(ECM_single_default_pulse.SOC);
ECM_single_default_pulse.V = squeeze(ECM_single_default_pulse.V);


ECM_single_default_pulse.SLIDE = readmatrix('../../results/Cell_ECM_single_default_pulse_Cell_ECM_cellData.csv','NumHeaderLines',3);

% Plotting:

figure;
subplot(3,1,1);
plot(ECM_single_default_pulse.tout, ECM_single_default_pulse.V,'LineWidth',1.5); hold on;
plot(ECM_single_default_pulse.SLIDE(:,5), ECM_single_default_pulse.SLIDE(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Terminal voltage [V]')
grid on;

subplot(3,1,2);
plot(ECM_single_default_pulse.tout, ECM_single_default_pulse.SOC,'LineWidth',1.5); hold on;
plot(ECM_single_default_pulse.SLIDE(:,5), ECM_single_default_pulse.SLIDE(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(3,1,3);
plot(ECM_single_default_pulse.tout, ECM_single_default_pulse.I,'LineWidth',1.5); hold on;
plot(ECM_single_default_pulse.SLIDE(:,5), ECM_single_default_pulse.SLIDE(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

sgtitle('ECM single default pulse')





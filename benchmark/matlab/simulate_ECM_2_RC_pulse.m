% This file is for simulating ECM cases
clear variables; close all; clc;

Tminutes = 5;

pulse = sim('Cell_ECM_2_RC_single_default_pulse.slx');

pulse.I = squeeze(pulse.I);
pulse.SOC = squeeze(pulse.SOC);
pulse.V = squeeze(pulse.V);
%%

pulse.SLIDE = readmatrix('../../results/Cell_ECM_2_RC_single_default_pulse_Cell_ECM_cellData.csv','NumHeaderLines',3);

% Plotting:

figure;
subplot(3,1,1);
plot(pulse.tout, pulse.V,'LineWidth',1.5); hold on;
plot(pulse.SLIDE(:,5), pulse.SLIDE(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Terminal voltage [V]')
grid on;

subplot(3,1,2);
plot(pulse.tout, pulse.SOC,'LineWidth',1.5); hold on;
plot(pulse.SLIDE(:,5), pulse.SLIDE(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(3,1,3);
plot(pulse.tout, pulse.I,'LineWidth',1.5); hold on;
plot(pulse.SLIDE(:,5), pulse.SLIDE(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

sgtitle('ECM single 2 RC pair default pulse')





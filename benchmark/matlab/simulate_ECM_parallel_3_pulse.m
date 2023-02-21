% This file is for simulating ECM cases
clear variables; close all; clc;

Tminutes = 5;

pulse = sim('Cell_ECM_parallel_3_default_pulse.slx');

pulse.I1 = squeeze(pulse.I1);
pulse.I2 = squeeze(pulse.I2);
pulse.I3 = squeeze(pulse.I3);

pulse.I = squeeze(pulse.I);
pulse.SOC = squeeze(pulse.SOC);
pulse.V = squeeze(pulse.V);

%%
pulse.SLIDE_1 = readmatrix('../../results/Cell_ECM_parallel_3_default_pulse_Par3_1_cellData.csv','NumHeaderLines',3);
pulse.SLIDE_2 = readmatrix('../../results/Cell_ECM_parallel_3_default_pulse_Par3_2_cellData.csv','NumHeaderLines',3);
pulse.SLIDE_3 = readmatrix('../../results/Cell_ECM_parallel_3_default_pulse_Par3_3_cellData.csv','NumHeaderLines',3);



% Plotting:

figure;
subplot(5,1,1);
plot(pulse.tout, pulse.V,'LineWidth',1.5); hold on;
plot(pulse.SLIDE_1(:,5), pulse.SLIDE_1(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Terminal voltage [V]')
grid on;

subplot(5,1,2);
plot(pulse.tout, pulse.SOC,'LineWidth',1.5); hold on;
plot(pulse.SLIDE_1(:,5), pulse.SLIDE_1(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(5,1,3);
plot(pulse.tout, pulse.I1,'LineWidth',1.5); hold on;
plot(pulse.SLIDE_1(:,5), pulse.SLIDE_1(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

subplot(5,1,4);
plot(pulse.tout, pulse.I2,'LineWidth',1.5); hold on;
plot(pulse.SLIDE_2(:,5), pulse.SLIDE_2(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

subplot(5,1,5);
plot(pulse.tout, pulse.I3,'LineWidth',1.5); hold on;
plot(pulse.SLIDE_3(:,5), pulse.SLIDE_3(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

sgtitle('ECM parallel Cell-1 pulse')






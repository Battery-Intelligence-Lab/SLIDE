% This file is for simulating Bucket cases
clear variables; close all; clc;

Tminutes = 5;

Bucket_single_default_pulse = sim('Cell_Bucket_single_default_pulse.slx');

Bucket_single_default_pulse.I = squeeze(Bucket_single_default_pulse.I);
Bucket_single_default_pulse.SOC = squeeze(Bucket_single_default_pulse.SOC);
Bucket_single_default_pulse.V = squeeze(Bucket_single_default_pulse.V);


Bucket_single_default_pulse.SLIDE = readmatrix('../../results/Cell_Bucket_single_default_pulse_Cell_ECM_cellData.csv','NumHeaderLines',3);

% Plotting:

figure;
subplot(3,1,1);
plot(Bucket_single_default_pulse.tout, Bucket_single_default_pulse.V,'LineWidth',1.5); hold on;
plot(Bucket_single_default_pulse.SLIDE(:,5), Bucket_single_default_pulse.SLIDE(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Terminal voltage [V]')
grid on;

subplot(3,1,2);
plot(Bucket_single_default_pulse.tout, Bucket_single_default_pulse.SOC,'LineWidth',1.5); hold on;
plot(Bucket_single_default_pulse.SLIDE(:,5), Bucket_single_default_pulse.SLIDE(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(3,1,3);
plot(Bucket_single_default_pulse.tout, Bucket_single_default_pulse.I,'LineWidth',1.5); hold on;
plot(Bucket_single_default_pulse.SLIDE(:,5), Bucket_single_default_pulse.SLIDE(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

sgtitle('Bucket single default pulse')





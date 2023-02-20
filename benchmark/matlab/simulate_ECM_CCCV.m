% This file is for simulating ECM cases
clear variables; close all; clc;

Tminutes = 5;
tic; 
ECM_single_default_CCCV = sim('Cell_ECM_single_default_CCCV.slx');
toc

ECM_single_default_CCCV.I = squeeze(ECM_single_default_CCCV.I);
ECM_single_default_CCCV.SOC = squeeze(ECM_single_default_CCCV.SOC);
ECM_single_default_CCCV.V = squeeze(ECM_single_default_CCCV.V);



%% Plotting:
ECM_single_default_CCCV.SLIDE = readmatrix('../../results/Cell_ECM_single_default_CCCV_Cell_ECM_cellData.csv','NumHeaderLines',3);

figure;
subplot(3,1,1);
plot(ECM_single_default_CCCV.tout, ECM_single_default_CCCV.V,'LineWidth',1.5); hold on;
plot(ECM_single_default_CCCV.SLIDE(:,5), ECM_single_default_CCCV.SLIDE(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Terminal voltage [V]')
grid on;

subplot(3,1,2);
plot(ECM_single_default_CCCV.tout, ECM_single_default_CCCV.SOC,'LineWidth',1.5); hold on;
plot(ECM_single_default_CCCV.SLIDE(:,5), ECM_single_default_CCCV.SLIDE(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(3,1,3);
plot(ECM_single_default_CCCV.tout, ECM_single_default_CCCV.I,'LineWidth',1.5); hold on;
plot(ECM_single_default_CCCV.SLIDE(:,5), ECM_single_default_CCCV.SLIDE(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

sgtitle('ECM single 1 RC pair default CCCV')

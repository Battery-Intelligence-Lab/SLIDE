% This file is for simulating ECM cases
clear variables; close all; clc;

Tminutes = 5;
tic; 
CCCV = sim('Cell_ECM_2_RC_single_default_CCCV.slx');
toc

CCCV.I = squeeze(CCCV.I);
CCCV.SOC = squeeze(CCCV.SOC);
CCCV.V = squeeze(CCCV.V);



%% Plotting:
CCCV.SLIDE = readmatrix('../../results/Cell_ECM_2_RC_single_default_CCCV_Cell_ECM_cellData.csv','NumHeaderLines',3);

figure;
subplot(3,1,1);
plot(CCCV.tout, CCCV.V,'LineWidth',1.5); hold on;
plot(CCCV.SLIDE(:,5), CCCV.SLIDE(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Terminal voltage [V]')
grid on;

subplot(3,1,2);
plot(CCCV.tout, CCCV.SOC,'LineWidth',1.5); hold on;
plot(CCCV.SLIDE(:,5), CCCV.SLIDE(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(3,1,3);
plot(CCCV.tout, CCCV.I,'LineWidth',1.5); hold on;
plot(CCCV.SLIDE(:,5), CCCV.SLIDE(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

sgtitle('ECM single 2 RC pair default CCCV')

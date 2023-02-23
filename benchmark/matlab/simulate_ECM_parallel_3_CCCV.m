% This file is for simulating ECM cases
clear variables; close all; clc;

Tminutes = 5;
tic; 
CCCV = sim('Cell_ECM_parallel_3_default_CCCV.slx');
toc

CCCV.I1 = squeeze(CCCV.I1);
CCCV.I2 = squeeze(CCCV.I2);
CCCV.I3 = squeeze(CCCV.I3);

CCCV.I = squeeze(CCCV.I);
CCCV.SOC = squeeze(CCCV.SOC);
CCCV.V = squeeze(CCCV.V);

%%
CCCV.SLIDE_1 = readmatrix('../../results/Cell_ECM_parallel_3_default_CCCV_Par3_1_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_2 = readmatrix('../../results/Cell_ECM_parallel_3_default_CCCV_Par3_2_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_3 = readmatrix('../../results/Cell_ECM_parallel_3_default_CCCV_Par3_3_cellData.csv','NumHeaderLines',3);



% Plotting:

figure;
subplot(5,1,1);
plot(CCCV.tout, CCCV.V,'LineWidth',1.5); hold on;
plot(CCCV.SLIDE_1(:,5), CCCV.SLIDE_1(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Terminal voltage [V]')
grid on;

subplot(5,1,2);
plot(CCCV.tout, CCCV.SOC,'LineWidth',1.5); hold on;
plot(CCCV.SLIDE_1(:,5), CCCV.SLIDE_1(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(5,1,3);
plot(CCCV.tout, CCCV.I1,'LineWidth',1.5); hold on;
plot(CCCV.SLIDE_1(:,5), CCCV.SLIDE_1(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;
title('Cell-1')

subplot(5,1,4);
plot(CCCV.tout, CCCV.I2,'LineWidth',1.5); hold on;
plot(CCCV.SLIDE_2(:,5), CCCV.SLIDE_2(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;
title('Cell-2')

subplot(5,1,5);
plot(CCCV.tout, CCCV.I3,'LineWidth',1.5); hold on;
plot(CCCV.SLIDE_3(:,5), CCCV.SLIDE_3(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;
title('Cell-3')

sgtitle('ECM parallel 3 cells CCCV')






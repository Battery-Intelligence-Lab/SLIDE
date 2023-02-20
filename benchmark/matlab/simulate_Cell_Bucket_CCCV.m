% This file is for simulating Bucket cases
clear variables; close all; clc;

Tminutes = 5;
tic; 
Bucket_single_default_CCCV = sim('Cell_Bucket_single_default_CCCV.slx');
toc

Bucket_single_default_CCCV.I = squeeze(Bucket_single_default_CCCV.I);
Bucket_single_default_CCCV.SOC = squeeze(Bucket_single_default_CCCV.SOC);
Bucket_single_default_CCCV.V = squeeze(Bucket_single_default_CCCV.V);



%% Plotting:
Bucket_single_default_CCCV.SLIDE = readmatrix('../../results/Cell_Bucket_single_default_CCCV_Cell_ECM_cellData.csv','NumHeaderLines',3);

figure;
subplot(3,1,1);
plot(Bucket_single_default_CCCV.tout, Bucket_single_default_CCCV.V,'LineWidth',1.5); hold on;
plot(Bucket_single_default_CCCV.SLIDE(:,5), Bucket_single_default_CCCV.SLIDE(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Terminal voltage [V]')
grid on;

subplot(3,1,2);
plot(Bucket_single_default_CCCV.tout, Bucket_single_default_CCCV.SOC,'LineWidth',1.5); hold on;
plot(Bucket_single_default_CCCV.SLIDE(:,5), Bucket_single_default_CCCV.SLIDE(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(3,1,3);
plot(Bucket_single_default_CCCV.tout, Bucket_single_default_CCCV.I,'LineWidth',1.5); hold on;
plot(Bucket_single_default_CCCV.SLIDE(:,5), Bucket_single_default_CCCV.SLIDE(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

sgtitle('Bucket single 0 RC pair default CCCV')

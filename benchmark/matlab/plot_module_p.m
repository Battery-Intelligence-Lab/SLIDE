% This file plots SLIDE vs Simulink outputs for comparison
% Author(s): Volkan Kumtepeli
%      Date: 2023.03.01
function [] = plot_module_p(simOut, test)

figure;
subplot(5,1,1);
plot(simOut.tout, simOut.v1,'LineWidth',1.5); hold on;
plot(simOut.SLIDE_1(:,5), simOut.SLIDE_1(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Cell-1 voltage [V]')
grid on;

subplot(5,1,2);
plot(simOut.tout, simOut.SOC,'LineWidth',1.5); hold on;
plot(simOut.SLIDE_1(:,5), simOut.SLIDE_1(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(5,1,3);
plot(simOut.tout, simOut.I1,'LineWidth',1.5); hold on;
plot(simOut.SLIDE_1(:,5), simOut.SLIDE_1(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;
title('Cell-1')

subplot(5,1,4);
plot(simOut.tout, simOut.I2,'LineWidth',1.5); hold on;
plot(simOut.SLIDE_2(:,5), simOut.SLIDE_2(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;
title('Cell-2')

subplot(5,1,5);
plot(simOut.tout, simOut.I3,'LineWidth',1.5); hold on;
plot(simOut.SLIDE_3(:,5), simOut.SLIDE_3(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;
title('Cell-3')

sgtitle(test.Name)
end
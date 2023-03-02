% This file plots SLIDE vs Simulink simOutputs for comparison
% Plots single cell results. 
% Author(s): Volkan Kumtepeli
%      Date: 2023.03.02
function [] = plot_single(simOut, test)

figure;
subplot(3,1,1);
plot(simOut.tout, simOut.V,'LineWidth',1.5); hold on;
plot(simOut.SLIDE(:,5), simOut.SLIDE(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Terminal voltage [V]')
grid on;

subplot(3,1,2);
plot(simOut.tout, simOut.SOC,'LineWidth',1.5); hold on;
plot(simOut.SLIDE(:,5), simOut.SLIDE(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(3,1,3);
plot(simOut.tout, simOut.I,'LineWidth',1.5); hold on;
plot(simOut.SLIDE(:,5), simOut.SLIDE(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;

sgtitle(test.Name)

end
% This file plots SLIDE vs Simulink outputs for comparison
% Author(s): Volkan Kumtepeli
%      Date: 2023.03.01
function out = plot_variables(out)

figure;
subplot(5,1,1);
plot(out.tout, out.v1,'LineWidth',1.5); hold on;
plot(out.SLIDE_1(:,5), out.SLIDE_1(:,2),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Cell-1 voltage [V]')
grid on;

subplot(5,1,2);
plot(out.tout, out.SOC,'LineWidth',1.5); hold on;
plot(out.SLIDE_1(:,5), out.SLIDE_1(:,3),'--','LineWidth',1.5);
ylabel('SOC [%]');
legend('MATLAB', 'SLIDE');
grid on;

subplot(5,1,3);
plot(out.tout, out.I1,'LineWidth',1.5); hold on;
plot(out.SLIDE_1(:,5), out.SLIDE_1(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;
title('Cell-1')

subplot(5,1,4);
plot(out.tout, out.I2,'LineWidth',1.5); hold on;
plot(out.SLIDE_2(:,5), out.SLIDE_2(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;
title('Cell-2')

subplot(5,1,5);
plot(out.tout, out.I3,'LineWidth',1.5); hold on;
plot(out.SLIDE_3(:,5), out.SLIDE_3(:,1),'--','LineWidth',1.5);
legend('MATLAB', 'SLIDE');
ylabel('Current [A]')
grid on;
title('Cell-3')

sgtitle('ECM parallel 3 cells out')
end
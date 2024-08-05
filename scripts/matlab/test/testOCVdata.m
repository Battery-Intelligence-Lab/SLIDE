% This is a test function for OCV 

clear all; close all; clc;


table_data = readtable('Kokam_NMC_16Ah_OCV.xlsx');
data = table2array(table_data);

%%
figure; 
plot(data(:,7), data(:,2),'LineWidth',2)
xlim([-1,101]);
grid on;
hold on;

plot(data(:,12), data(:,2),'--','LineWidth',2);
legend('By surface Li', 'By Z');

xlabel('SOC [%]');
ylabel('Voltage [V]')

%%

V_0_100 = interp1(data(:,7), data(:,2), (0:10:100)');
SOC_0_100 = (0:10:100)';

for i=1:length(SOC_0_100)
    fprintf('%3.0f%%   %4.4f [V]\n', SOC_0_100(i),V_0_100(i));
end

zp_0_100 = interp1(data(:,7), data(:,8), (0:10:100)');
zn_0_100 = interp1(data(:,7), data(:,9), (0:10:100)');


fp_0_100 = interp1(data(:,7), data(:,3), (0:10:100)');
fn_0_100 = interp1(data(:,7), data(:,4), (0:10:100)');

fprintf('%4.6f, ',fp_0_100);
fprintf('\n');
fprintf('%4.6f, ',fn_0_100);


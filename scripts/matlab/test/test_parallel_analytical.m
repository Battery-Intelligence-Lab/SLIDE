close all; clc;

folder = '';


Npar = 10; 

results = cell(Npar,1);

tic;
for i = 1:Npar
    results{i} = slide_to_table("Cycler1_parECM_cell" + (i-1) + ".slide");
end
toc

%%

figure; 
for i=1:Npar    
    c = results{i};
    plot(c{:, "time [s]"},c{:, "Voltage [V]"});
    hold on;
end

grid on;


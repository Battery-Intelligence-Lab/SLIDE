clc;

folder = '';


Npar = 10; 

results = cell(Npar,1);

tic; %Cycler1_parECM_cell0 %parECM_boost_parECM_boost_cell
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
%%
figure; 
for i=1:Npar    
    c = results{i};
    plot(c{:, "time [s]"},c{:, "Current [A]"});
    hold on;
end

grid on;

%%
figure; 
for i=1:Npar    
    c = results{i};
    plot(c{:, "time [s]"},c{:, "SOC [-]"});
    hold on;
end

grid on;
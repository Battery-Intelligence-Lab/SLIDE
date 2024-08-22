clc;

folder = '';


Npar = 500; 

results = cell(Npar,1);

tic; %Cycler1_parECM_cell0 %parECM_boost_parECM_boost_cell
for i = 1:Npar
    results{i} = slide_to_table("Cycler1_parECM_cell" + (i-1) + ".slide");
end
toc

%%
width = 750;
length = 410;
line_width = 2; 
font_size = 18; 
title_size = 20;    
color1 = '#009083';

% Plot for Voltage vs. Time
figure('Position', [100, 100, width, length]); % Set the window size
hold on;
for i = 1:Npar    
    c = results{i};
    if i == 1
        h1 = plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 2
        h2 = plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 3
        h3 = plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width);
    elseif i == Npar
        hn = plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)], 'color', color1);
    else
        plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width);
    end
end
xlabel('Time [s]', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('Voltage [V]', 'FontSize', font_size, 'Interpreter', 'latex');
title('Voltages at Each Branch Over Time', 'FontSize', title_size, 'Interpreter', 'latex');
legend([h1 h2 h3 hn], {'Cell 1', 'Cell 2', '\vdots', ['Cell ', num2str(Npar)]}, 'interpreter','latex', 'FontSize', font_size, 'Location', 'best');
grid on;
% g = gca;
% set(g, 'fontsize', font_size);

% Plot for Current vs. Time
figure('Position', [100, 100, width, length]); % Set the window size
hold on;
for i = 1:Npar    
    c = results{i};
    if i == 1
        h1 = plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 2
        h2 = plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 3
        h3 = plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width);
    elseif i == Npar
        hn = plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)],'color', color1);
    else
        plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width);
    end
end
xlabel('Time [s]', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('Current [A]', 'FontSize', font_size, 'Interpreter', 'latex');
title('Currents at Each Branch Over Time', 'FontSize', title_size, 'Interpreter', 'latex');
legend([h1 h2 h3 hn], {'Cell 1', 'Cell 2', '\vdots', ['Cell ', num2str(Npar)]}, 'interpreter','latex', 'FontSize', font_size, 'Location', 'best');
grid on;

% Plot for SOC vs. Time
figure('Position', [100, 100, width, length]); % Set the window size
hold on;
for i = 1:Npar    
    c = results{i};
    if i == 1
        h1 = plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 2
        h2 = plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 3
        h3 = plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width);
    elseif i == Npar
        hn = plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)], 'color', color1);
    else
        plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width);
    end
end
xlabel('Time [s]', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('State of Charge (SOC)', 'FontSize', font_size, 'Interpreter', 'latex');
title('SOC at Each Branch Over Time', 'FontSize', title_size, 'Interpreter', 'latex');
legend([h1 h2 h3 hn], {'Cell 1', 'Cell 2', '\vdots', ['Cell ', num2str(Npar)]}, 'interpreter','latex', 'FontSize', font_size, 'Location', 'best');
grid on;

%%
% Create a figure with subplots
figure('Position', [10, 10, 800, 3*390]); % Adjust the window size to fit three subplots
font_size = 14; 
title_size = 16;  
legend_size = 12;
% Subplot 1: Voltage vs. Time
subplot(3, 1, 1); % Create the first subplot in a 3x1 grid
sgtitle(['Implementation for ', num2str(Npar), ' Cells \& 200 Cycles'], 'FontSize', 18, 'Interpreter', 'latex')
hold on;
for i = 1:Npar    
    c = results{i};
    if i == 1
        h1 = plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 2
        h2 = plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 3
        h3 = plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width);
    elseif i == Npar
        hn = plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)], 'color', color1);
    else
        plot(c{:, "time [s]"}, c{:, "Voltage [V]"}, 'LineWidth', line_width);
    end
end
xlabel('Time [s]', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('Voltage [V]', 'FontSize', font_size, 'Interpreter', 'latex');
title('Voltages at Each Branch Over Time', 'FontSize', title_size, 'Interpreter', 'latex');
legend([h1 h2 h3 hn], {'Cell 1', 'Cell 2', '\vdots', ['Cell ', num2str(Npar)]}, 'interpreter','latex', 'FontSize', legend_size, 'Location', 'best');
grid on;

% Subplot 2: Current vs. Time
subplot(3, 1, 2); % Create the second subplot in a 3x1 grid
hold on;
for i = 1:Npar    
    c = results{i};
    if i == 1
        h1 = plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 2
        h2 = plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 3
        h3 = plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width);
    elseif i == Npar
        hn = plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)], 'color', color1);
    else
        plot(c{:, "time [s]"}, c{:, "Current [A]"}, 'LineWidth', line_width);
    end
end
xlabel('Time [s]', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('Current [A]', 'FontSize', font_size, 'Interpreter', 'latex');
title('Currents at Each Branch Over Time', 'FontSize', title_size, 'Interpreter', 'latex');
legend([h1 h2 h3 hn], {'Cell 1', 'Cell 2', '\vdots', ['Cell ', num2str(Npar)]}, 'interpreter','latex', 'FontSize', legend_size, 'Location', 'best');
grid on;

% Subplot 3: SOC vs. Time
subplot(3, 1, 3); % Create the third subplot in a 3x1 grid
hold on;
for i = 1:Npar    
    c = results{i};
    if i == 1
        h1 = plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 2
        h2 = plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)]);
    elseif i == 3
        h3 = plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width);
    elseif i == Npar
        hn = plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width, 'DisplayName', ['Cell ', num2str(i)], 'color', color1);
    else
        plot(c{:, "time [s]"}, c{:, "SOC [-]"}, 'LineWidth', line_width);
    end
end
xlabel('Time [s]', 'FontSize', font_size, 'Interpreter', 'latex');
ylabel('State of Charge (SOC)', 'FontSize', font_size, 'Interpreter', 'latex');
title('SOC at Each Branch Over Time', 'FontSize', title_size, 'Interpreter', 'latex');
legend([h1 h2 h3 hn], {'Cell 1', 'Cell 2', '\vdots', ['Cell ', num2str(Npar)]}, 'interpreter','latex', 'FontSize', legend_size, 'Location', 'best');
grid on;
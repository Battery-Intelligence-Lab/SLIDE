% This file is for simulating ECM cases
clear variables; close all; clc;

% Default settings:
testNow = defaultSettings();

%% Case: Parallel 3 ECM Cell without contact resistances: 
testNow.Name = "Parallel 3 ECM Cell CCCV w/o contact resistances";
testNow.Rc = 0*[0.5, 1, 0.7]*1e-3; % % [Ohm] contact resistance. 
testNow.type = "module_p"; 

CCCV = run_test(testNow);

% Loading:
CCCV.SLIDE_1 = readmatrix('../../results/Cell_ECM_parallel_3_default_CCCV_Par3_1_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_2 = readmatrix('../../results/Cell_ECM_parallel_3_default_CCCV_Par3_2_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_3 = readmatrix('../../results/Cell_ECM_parallel_3_default_CCCV_Par3_3_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(CCCV, testNow);

%% Case: Parallel 3 ECM Cell with contact resistances: 
testNow.Name = "Parallel 3 ECM Cell CCCV with contact resistances";
testNow.Rc = [0.5, 1, 0.7]*1e-3; % [Ohm] contact resistance. 
testNow.Procedure = 1; % 0-> Pulse, 1-> CCCV

CCCV = run_test(testNow);

% Loading:
CCCV.SLIDE_1 = readmatrix('../../results/Cell_ECM_parallel_3_withRcontact_CCCV_Par3_1_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_2 = readmatrix('../../results/Cell_ECM_parallel_3_withRcontact_CCCV_Par3_2_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_3 = readmatrix('../../results/Cell_ECM_parallel_3_withRcontact_CCCV_Par3_3_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(CCCV, testNow);

%%
testNow = defaultSettings();
testNow.Name = "Parallel 3 ECM Cell pulse w/o contact resistances";
testNow.type = "module_p"; 
testNow.Rc = 0*testNow.Rc;
testNow.R0 = [1,2,2]*1e-3; 
testNow.Procedure = 0;
testNow.Tend = 10*60*3;

pulse = run_test(testNow);

% Loading:
pulse.SLIDE_1 = readmatrix('../../results/Cell_ECM_parallel_3_default_pulse_Par3_1_cellData.csv','NumHeaderLines',3);
pulse.SLIDE_2 = readmatrix('../../results/Cell_ECM_parallel_3_default_pulse_Par3_2_cellData.csv','NumHeaderLines',3);
pulse.SLIDE_3 = readmatrix('../../results/Cell_ECM_parallel_3_default_pulse_Par3_3_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(pulse, testNow);
% This file is for simulating ECM cases
clear variables; close all; clc;

% Default settings:
testNow = defaultSettings();

%% Case: Parallel 3 ECM Cell without contact resistances: 
testNow.Rc = 0*[0.5, 1, 0.7]*1e-3; % % [Ohm] contact resistance. 

tic; 
CCCV = sim('Cell_ECM_parallel_3_with_Rcontact_CCCV.slx');
toc

CCCV = squeeze_variables(CCCV);

% Loading:
CCCV.SLIDE_1 = readmatrix('../../results/Cell_ECM_parallel_3_default_CCCV_Par3_1_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_2 = readmatrix('../../results/Cell_ECM_parallel_3_default_CCCV_Par3_2_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_3 = readmatrix('../../results/Cell_ECM_parallel_3_default_CCCV_Par3_3_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(CCCV);

%% Case: Parallel 3 ECM Cell with contact resistances: 
testNow.Rc = [0.5, 1, 0.7]*1e-3; %#ok<NASGU> % [Ohm] contact resistance. 
testNow.Procedure = 1; % 0-> Pulse, 1-> CCCV

tic; 
CCCV = sim('Cell_ECM_parallel_3_with_Rcontact_CCCV.slx');
toc

CCCV = squeeze_variables(CCCV);

% Loading:
CCCV.SLIDE_1 = readmatrix('../../results/Cell_ECM_parallel_3_withRcontact_CCCV_Par3_1_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_2 = readmatrix('../../results/Cell_ECM_parallel_3_withRcontact_CCCV_Par3_2_cellData.csv','NumHeaderLines',3);
CCCV.SLIDE_3 = readmatrix('../../results/Cell_ECM_parallel_3_withRcontact_CCCV_Par3_3_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(CCCV);

%%
testNow = defaultSettings();
testNow.Rc = 0*testNow.Rc;
testNow.R0 = [1,2,2]*1e-3; 
testNow.Procedure = 0;
testNow.Tend = 10*60*4;
pulse = sim('Cell_ECM_parallel_3_with_Rcontact_CCCV.slx');
pulse = squeeze_variables(pulse);

% Loading:
pulse.SLIDE_1 = readmatrix('../../results/Cell_ECM_parallel_3_default_pulse_Par3_1_cellData.csv','NumHeaderLines',3);
pulse.SLIDE_2 = readmatrix('../../results/Cell_ECM_parallel_3_default_pulse_Par3_2_cellData.csv','NumHeaderLines',3);
pulse.SLIDE_3 = readmatrix('../../results/Cell_ECM_parallel_3_default_pulse_Par3_3_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(pulse);


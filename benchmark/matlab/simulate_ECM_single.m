% This file is for simulating ECM cases
clear variables; close all; clc;

% Default settings:
testNow = defaultSettings();
testNow.type = "single"; 

%% Single ECM (2 RC pair) pulse: 
testNow.Name = "Single ECM (2 RC pair) pulse";
testNow.Procedure = 0; % 0 -> pulse, 1-> CCCV;
testNow.Tend = 30*60; % [seconds]. 
testNow.R0(1) = 2e-3;
simOut = run_test(testNow);

% Loading:
simOut.SLIDE = readmatrix('../../results/Cell_ECM_2_RC_single_default_pulse_Cell_ECM_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(simOut, testNow);

%% Single ECM (2 RC pair) CCCV: 
testNow.Name = "Single ECM (2 RC pair) CCCV";
testNow.Procedure = 1; % 0 -> pulse, 1-> CCCV;
testNow.Tend = 42724; % [seconds]. 
simOut = run_test(testNow);

% Loading:
simOut.SLIDE = readmatrix('../../results/Cell_ECM_2_RC_single_default_CCCV_Cell_ECM_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(simOut, testNow);

%% Single ECM pulse:
testNow.Name = "Single ECM (1 RC pair) pulse";
testNow.R2(1) = 1e-9; % Make RC pair's RC zero.
testNow.Procedure = 0; % 0 -> pulse
testNow.Tend = 30*60; % 30 minutes. 

ECM_single = run_test(testNow);

ECM_single.SLIDE = readmatrix('../../results/Cell_ECM_single_default_pulse_Cell_ECM_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(ECM_single, testNow);

%% Single ECM CCCV:
testNow.Name = "Single ECM (1 RC pair) CCCV";
testNow.Procedure = 1; % 1 -> CCCV;
testNow.Tend = 42724; % [seconds]
simOut = run_test(testNow);

% Loading:
simOut.SLIDE = readmatrix('../../results/Cell_ECM_single_default_CCCV_Cell_ECM_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(simOut, testNow);

%% Single ECM-Bucket (0 RC pair) pulse: 
testNow.Name = "Single ECM-Bucket (0 RC pair) pulse";
testNow.Procedure = 0; % 0 -> pulse, 1-> CCCV;
testNow.Tend = 30*60; % 30 minutes. 
testNow.R1(1) = 1e-9; % Make RC pair's RC zero.
simOut = run_test(testNow);

% Loading:
simOut.SLIDE = readmatrix('../../results/Cell_Bucket_single_default_pulse_Cell_ECM_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(simOut, testNow);

%% Single ECM-Bucket (0 RC pair) CCCV: 
testNow.Name = "Single ECM-Bucket (0 RC pair) CCCV";
testNow.Procedure = 1; % 0 -> pulse, 1-> CCCV;
testNow.Tend = 17450; % [seconds]. 
simOut = run_test(testNow);

% Loading:
simOut.SLIDE = readmatrix('../../results/Cell_Bucket_single_default_CCCV_Cell_ECM_cellData.csv','NumHeaderLines',3);

% Plotting:
plot_variables(simOut, testNow);





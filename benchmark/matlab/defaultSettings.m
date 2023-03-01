function default = defaultSettings()

%% Case: Parallel 3 ECM Cell with contact resistances: 
default.Tend = 19*60*60; % [seconds] 
default.Rc = [0.5, 1, 0.7]*1e-3; %#ok<NASGU> % [Ohm] contact resistance. 
default.R0 = [1, 3, 2]*1e-3; % [Ohm] series resistance. 
default.R1 = [15.8, 15.8, 15.8]*1e-3; %[Ohm] 1st RC pair's R. 
default.Tau1 = 38e3*default.R1;
default.Procedure = 1; % 0-> Pulse, 1-> CCCV
default.PulseWidth = 2*60*5; 
default.PulseAmplitude = 16;
end
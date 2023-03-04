% This function runs the Simulink file corresponding to the case
% interested. 
%
% Author(s): Volkan Kumtepeli
%      Date: 2023.03.01
function out = run_test(testNow)
if(strcmpi(testNow.type,'single'))
    simIn  = Simulink.SimulationInput('Cell_ECM_single');

elseif(strcmpi(testNow.type,'module_p'))
    simIn  = Simulink.SimulationInput('Cell_ECM_parallel');

elseif(strcmpi(testNow.type,'module_s'))
    simIn  = Simulink.SimulationInput('Cell_ECM_series');
else
    error(strcat('Test case ', testNow.type, ' is NOT found!'));
end

simIn = setVariable(simIn,'testNow',testNow);

tic; 
out = sim(simIn);
toc

out = squeeze_variables(out);

end
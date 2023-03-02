% This file plots SLIDE vs Simulink outputs for comparison
% Author(s): Volkan Kumtepeli
%      Date: 2023.03.01
function [] = plot_variables(simOut, test)

if(strcmpi(test.type,'single'))
    plot_single(simOut, test);
elseif(strcmpi(test.type,'module_p'))
    plot_module_p(simOut, test);
elseif(strcmpi(test.type,'module_s'))
    plot_module_p(simOut, test);
else
    error(strcat('Test case ', test.type, ' is NOT found!'));
end

end
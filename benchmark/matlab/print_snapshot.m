% This function prints the snapshot values for snapshot testing.
%
% Author(s): Volkan Kumtepeli
%      Date: 2023.03.02
function [] = print_snapshot(simOut, test, folder)

if(nargin<3)
    folder = '.';
end

results_matrix = [simOut.tout, simOut.I, simOut.V, simOut.I1, simOut.I2, simOut.I3, simOut.v1, simOut.v2, simOut.v3];

path = fullfile(folder, "test.csv");

if(strcmpi(test.type,'module_p'))
    writematrix(results_matrix(1:40:end,:), path);
else
    writematrix(results_matrix(1:40:end,1:3), path);

end


end
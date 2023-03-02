% This file reduces the dimensionality of output signals. 
% Author(s): Volkan Kumtepeli
%      Date: 2023.03.01
function out = squeeze_variables(out)

out.I1 = squeeze(out.I1);
out.I2 = squeeze(out.I2);
out.I3 = squeeze(out.I3);

out.v1 = squeeze(out.v1);
out.v2 = squeeze(out.v2);
out.v3 = squeeze(out.v3);

out.I = squeeze(out.I);
out.SOC = squeeze(out.SOC);
out.V = squeeze(out.V);

end
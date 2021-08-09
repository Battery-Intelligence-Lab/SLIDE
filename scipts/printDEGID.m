
function str = printDEGID(sei_id, sei_por, CS_id, CS_diff, LAM_id, pl_id)
% Function to get the string from the degradatio models. It works similar
% to the function in C++: the user specifies which degradation models were
% used, and this function makes the string used to identify said
% combination of degradation models
%
% IN
% sei_id    array with the identifications of the SEI models
% sei_por   integer deciding whether we reduce the volume fraction due to SEI growth
% CS_id     array with the identifications of the crack growth models
% CS_diff   integer deciding whether we reduce the negative diffusion constant due to surface cracks
% LAM_id    array with the identifications of the LAM models
% pl_id     array with the identifications of the lithium plating models
% 
% OUT
% str       string-representation of the degradation settings
%           this string will be part of the names of the subfolder in which
%           the simulation results will be written by C++
%
%
% Copyright (c) 2019, The Chancellor, Masters and Scholars of the University 
% of Oxford, VITO nv, and the 'Slide' Developers.
% See the licence file LICENCE.txt for more information.
%%

str = '';   % start with an empty string

% Append the SEI models, separated by a hyphen
for i=1:length(sei_id)
    str = strcat(str, num2str(sei_id(i)),'-');    
end

% Append SEI_porosity, and terminate with an underscore to separate the SEI
% models from the crack growth
str = strcat(str, num2str(sei_por),'_');  

% Append the crack growth models, separated by a hyphen
for i=1:length(CS_id)
    str = strcat(str, num2str(CS_id(i)),'-');    
end

% Append CS_diff, and terminate with an underscore to separate the crack
% grawth models from the LAM models
str = strcat(str, num2str(CS_diff),'_');  

% Append the LAM models, separated by a hyphen
for i=1:length(LAM_id)
    str = strcat(str, num2str(LAM_id(i)));   
    if i < length(LAM_id)
        str = strcat(str,'-');  
    end
end

% Add an underscore to separate and ppend the lithium plating model
str = strcat(str,'_', num2str(pl_id));  














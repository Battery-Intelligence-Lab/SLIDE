function nodes = get_nodes( data,N)
%GET_NODES Computes the Chebyshev nodes and the mapping function to
%calculate the physical coordinate xp from the computation coordinates xc.
%
% INPUTS
% data          Structure containg the model parameters, see get_modelData
% N             The number N+1 of Chebyshev nodes used to discretize the
%               diffusion PDE in each particle.
%
% OUTPUTS
% nodes         Structure containing the Chebyshev nodes xr used to
%               discretize the diffusion equation in each particle and 
%               the function-mapping from the Chebyshev computation nodes
%               xc in [-1,1] to the physical nodel xp in [0,Rs]
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.

M = 2*N;
% Computational coordinates (Chebyshev nodes)
nodes.xm = chebdif(M+1,0);  %  Full domain [-R,R]->[-1,1]
nodes.xr = nodes.xm(1:N+1); %  Half domain [ 0,R]->[ 0,1]

% Coordinate mapping functions: Mapping from the computational coordinates 
% xc to the physical coordinates xp
nodes.xc2xp = @xc2xp;
    function xp = xc2xp(xc,domain)
        if  strcmp(domain,'r1') == 1
            xp = data.Rs1*xc;
        elseif  strcmp(domain,'r3') == 1
            xp = data.Rs3*xc;
        else
            error('ERROR: Domain must be r1 for anode or r3 for cathode');
        end
    end

end


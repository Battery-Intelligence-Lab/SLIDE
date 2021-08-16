% Script that calculates the model parameters for the Chebyshev
% discretisation
%
%
% Copyright (c) 2019, The Chancellor, Masters and Scholars of the University 
% of Oxford, VITO nv, and the 'Slide' Developers.
% See the licence file LICENCE.txt for more information.

clc
close all
clear

%% input parameters

% Adding required folder / subfolder on the MATLAB path
addpath(genpath('MatlabSetup'));

% Input parameters, must be the same as in the C++ model
nch = 5;            % Number of positive inner chebyshev nodes. 
Rp = 8.5*10^(-6);   % radius of the positive particle
Rn = 1.25*10^(-5);  % radius of the negative particle
  
N = nch + 1;        % number of positive chebyshev points including the surface node
                    % but excluding the centre node (and the negative half)
M = 2*N+1;         % total number of nodes, including the negative half and centre nodes

%% Get the location of the nodes
xm = chebdif(M,0);
x = xm(2:N);    % only the inner nodes

% Write the locations to 16 significant digits
dlmwrite('Cheb_Nodes.csv',x,'delimiter', ',', 'precision', 16); 

%% Get the matrices needed for the model

% Get matrices
matrices_spm    = get_model(N, xm, Rp, Rn);  % Create the SPM model
matrices_spm.Q(1,:) = 0;                               % due to numerical errors it is not exactly 0, but it is the integral from the centre to the centre, so it should be 0
Ccentre = (matrices_spm.DM(1,1:N,1) + matrices_spm.DM(1,N+2:M,1)*matrices_spm.P )';  % matrix to get the concentration at the centre node
VN = inv(matrices_spm.Vn);
VP = inv(matrices_spm.Vp);

% One of the eigenvectors represents a uniform concentration. The
% corresponding eigenvalue has to be 0 (because the time derivative of a 
% constant concentration has to be 0). 
% Due to numerical errors, the value is not exactly 0. So we need to force
% it to be 0. (The C++ code takes advantage of the fact it is exactly 0).
d = diag(matrices_spm.Ap);          % get the eigenvalues
d = d ./ max(abs(d));               % get relative values
zer = find(abs(d) < 10^(-10));      % the location of the 0 eigenvalue
matrices_spm.Ap(zer,zer) = 0;       % Set the values to exactly 0;
matrices_spm.An(zer,zer) = 0;
% 
%% Write the parameters
 csvwrite('Cheb_input.csv',[nch ; Rp ; Rn ; zer-1]); % write the input parameters for use in C++
%     % write zer-1 because C++ indices start at 0 while Matlab indices start
%     % at 1
 dlmwrite('Cheb_Ap.csv',diag(matrices_spm.Ap),'delimiter', ',', 'precision', 16);  
 dlmwrite('Cheb_An.csv',diag(matrices_spm.An),'delimiter', ',', 'precision', 16);    
 dlmwrite('Cheb_Bp.csv',matrices_spm.Bp,'delimiter', ',', 'precision', 16);   
 dlmwrite('Cheb_Bn.csv',matrices_spm.Bn,'delimiter', ',', 'precision', 16);     
 dlmwrite('Cheb_Cp.csv',matrices_spm.Cp,'delimiter', ',', 'precision', 16);   
 dlmwrite('Cheb_Cn.csv',matrices_spm.Cn,'delimiter', ',', 'precision', 16);     
 dlmwrite('Cheb_Cc.csv',Ccentre,'delimiter', ',', 'precision', 16);   
 dlmwrite('Cheb_Dp.csv',matrices_spm.Dp,'delimiter', ',', 'precision', 16);   
 dlmwrite('Cheb_Dn.csv',matrices_spm.Dn,'delimiter', ',', 'precision', 16);   
 dlmwrite('Cheb_Q.csv',matrices_spm.Q,'delimiter', ',', 'precision', 16);      
 
 dlmwrite('Cheb_Vn.csv',VN,'delimiter', ',', 'precision', 16);       
 dlmwrite('Cheb_Vp.csv',VP,'delimiter', ',', 'precision', 16);   

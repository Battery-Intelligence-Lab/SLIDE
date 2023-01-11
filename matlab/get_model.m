function model = get_model( N, xm, Rp, Rn )
%GET_MODEL This function computes the matrices A, B, C and D of the 
% particel diffusion state space model of the single particle model for
% each electrode:
%    du/dt  = A*u + B*j
%    cs     = C*u + D*j
% with:
%   u       the transformed concentration.
%   j       volumetric reaction rate.
%   cs      particle concentration profile (including at the surface)
% See the word document explaining the battery model and the Chebyshev
% discretisation for an explanation of the matrices and an explanation of
% where the formulas come from. 
% 2 main differences with that pdf: 
% 0<r<R while this code uses -R < r < R such that the indices are a bit
    % different (as well as a factor '2' and a minus sign from the different
    % transformation of variables-function)
% in the pdf, index 1 is the centre and index N+1 is the surface while here
    % we flip the matrices (the first row is the surface while the last row
    % is the centre)
%
%
% IN
% N             The number N+1 of Chebyshev nodes used to discretize the
%               diffusion PDE in each particle.
% xm            location of the chebyshev nodes
% Rp            Radius of the positive particle
% Rn            Radius of the negative particle
%
% OUT
% model         A structure constaining the matrices A, B, C and D for each
%               electrode, and also the Chebyshev differentiation matrices
%               DNr.
%               Also contains the Chebyshev integration matrix Q
%
%
% Copyright (c) 2019, The Chancellor, Masters and Scholars of the University 
% of Oxford, VITO nv, and the 'Slide' Developers.
% See the licence file LICENCE.txt for more information.

data.R_n = Rn;
data.R_p = Rp;

xp = xm*Rp;
xn = xm*Rn;

    M = 2*N;
    % Computing the Chebyshev differentiation matrices
    [~,DM] = chebdif(M+1,2);
    model.DM = DM;
    DM1 = squeeze(DM(1,:,1));
    DM2 = squeeze(DM(1:N,:,2));

    P   = fliplr(eye(N)); 
    model.P = P;  % Permutation matrix, backward-identity (flipped the indices)
    
%% first diffusion matrices
    % Modified differentation matrices accounting for the solution symmetry
    DN2 = DM2(:,1:N) - DM2(:,N+2:M+1)*P;
    DN1 = DM1(1,1:N) - DM1(1,N+2:M+1)*P;


    A = DN2(2:N,2:N) + (DN2(2:N,1)*DN1(1,2:N))/(1-DN1(1,1));
    B = DN2(2:N,1)/(1-DN1(1,1));
    C = DN1(1,2:N)/(1-DN1(1,1));
    D = 1/(1-DN1(1,1));

    % Anode diffusion state-space model
    A1 = (1/data.R_n^2)*A;
    B1 = B;
    C1 = vertcat( ...
                C/data.R_n , ...
                diag(1./xn(2:N)) );
    D1 = vertcat(...
                data.R_n*D , ...
                zeros(N-1,1) );

    % Cathode diffusion state-space model
    A3 = (1/data.R_p^2)*A;
    B3 = B;
    C3 = vertcat( ...
                C/data.R_p , ...
                diag(1./xp(2:N)) );
    D3 = vertcat( ...
                data.R_p*D , ...
                zeros(N-1,1) );

%% Transformation to eigenvector-basis z (instead of u)
    % Note: eigenvalue (i) = 0 = eigenvector for uniform concentration
    %           -> no diffusion going on (dz(i)/dt = 0 * z(i))
    %         if start with uniform concentration, only z0(i) != 0
    [Vn,Dn] = eig(A1);
    model.Ln = Dn;
    model.Vn = Vn;
    model.An = Dn;       % diagonal matrix with eigenvalues
    model.Bn = Vn\B1;    % V^(-1)*B
    model.Cn = C1*Vn;    % row1 = surface (node N+1), last row = close to centre (node 2)
    model.Dn = D1;

    [Vp,Dp] = eig(A3);
    model.Lp = Dp;
    model.Vp = Vp;
    model.Ap = Dp;       % diagonal matrix with eigenvalues
    model.Bp = Vp\B3;    % V^(-1)*B
    model.Cp = C3*Vp;    % row1 = surface (node N+1), last row = close to centre (node 2)
    model.Dp = D3;

%% Chebyshev integration matrix
    model.Q = cumsummat(M+1);

end


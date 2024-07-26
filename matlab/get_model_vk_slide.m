function model = get_model_vk_slide(nch)
%GET_MODEL This function computes the matrices A, B, C and D of the 
% particel diffusion state space model of the single particle model for
% each electrode:
%    du/dt  = A*u + B*j
%    cs     = C*u + D*j
% with:
%   u=rc    state of the diffusion model.
%   j       volumetric reaction rate.
%   cs      particle concentration profile (including at the surface)
%
% INPUTS
% nodes         Structure containing the Chebyshev nodes, see get_nodes
% N             The number N+1 of Chebyshev nodes used to discretize the
%               diffusion PDE in each particle.
%
% OUTPUTS
% model         A structure constaining the matrices A, B, C and D for each
%               electrode, and also the Chebyshev differentiation matrices
%               DNr and the Clenshaw-Curtis quadrature weights wn.
%
%
% Copyright (c) 2016, The Chancellor, Masters and Scholars of the University 
% of Oxford, and the 'Spectral li-ion SPM' Developers.
% See the licence file LICENCE.txt for more information.


% M+1: number of Chebyshev nodes in each particle on the domain [-R,R]
% mapped onto [-1,1], with N = M/2 defined by the user.


%-------------
% Solid particles of the electrodes
Rn = 12.5e-6; % Anode solid particles' radius [m]
Rp =  8.5e-6; % Cathode solid particles' radius [m]

N = nch + 1; 
M = 2*N;
model.N = N;
model.M = M;

Ncheb = M+1; 
Mcheb = 2;

dtheta = pi/(Ncheb-1);

% Computational coordinates (Chebyshev nodes)
xm = sin((Ncheb-1:-2:1-Ncheb)'*dtheta/2);  %  Full domain [-R,R]->[-1,1]
xr = xm(2:N); %  Half domain [ 0,R]->[ 0,1]

xp = xr*Rp;
xn = xr*Rn;


% Computing the Chebyshev differentiation matrices
D_vk = eye(Ncheb);
for ell = 1:Mcheb
    for i=1:Ncheb
        row_sum = 0;
        for j = 1:Ncheb
            if(i==j)
                continue;
            end

            DX_vk = cos(dtheta*(i-1)) - cos(dtheta*(j-1));

            C_vk = 1;

            if(i==1 || i==Ncheb)
                C_vk = 2*C_vk;
            end
    
            if(j==1 || j==Ncheb)
                C_vk = C_vk/2;
            end
    
            if(mod(i+j,2)==1)
                C_vk = -C_vk;
            end            


            D_vk(i,j) = ell*(C_vk*D_vk(i,i) - D_vk(i,j))/DX_vk;
            row_sum  =  row_sum - D_vk(i,j);
        end
        D_vk(i,i) = row_sum;
    end
    DM(:,:,ell) = D_vk;
end




model.DM = DM;
DM1 = squeeze(DM(1,:,1));
DM2 = squeeze(DM(1:N,:,2));

P   = fliplr(eye(N)); 
model.P = P;  % Permutation matrix, backward-identity (flipped the indices)

% Modified differentation matrices accounting for the solution symmetry
DN2 = DM2(:,1:N) - DM2(:,N+2:M+1)*P;
DN1 = DM1(1,1:N) - DM1(1,N+2:M+1)*P;


A = DN2(2:N,2:N) + (DN2(2:N,1)*DN1(1,2:N))/(1-DN1(1,1));
B = DN2(2:N,1)/(1-DN1(1,1));
C = DN1(1,2:N)/(1-DN1(1,1));
D = 1/(1-DN1(1,1));

% Anode diffusion state-space model
A1 = (1/Rn^2)*A;
B1 = B;
C1 = vertcat( C/Rn , ...
                    diag(1./xn ) ); 
% Coordinate mapping functions: Mapping from the computational coordinates 
% xc to the physical coordinates xp
D1 = vertcat(...
            Rn*D , ...
            zeros(N-1,1) );

% Cathode diffusion state-space model
A3 = (1/Rp^2)*A;
B3 = B;
C3 = vertcat( ...
            C/Rp, ...
            diag(1./xp) );
D3 = vertcat( ...
           	Rp*D , ...
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

%% Corrections?
    model.Q(1,:) = 0;        % due to numerical errors it is not exactly 0, but it is the integral from the centre to the centre, so it should be 0
    
    model.Cc = (model.DM(1,1:N,1) + model.DM(1,N+2:Ncheb,1)*model.P )';  % matrix to get the concentration at the centre node
    model.Vn = inv(model.Vn);
    model.Vp = inv(model.Vp);
    
    % One of the eigenvectors represents a uniform concentration. The
    % corresponding eigenvalue has to be 0 (because the time derivative of a 
    % constant concentration has to be 0). 
    % Due to numerical errors, the value is not exactly 0. So we need to force
    % it to be 0. (The C++ code takes advantage of the fact it is exactly 0).

  
    d = diag(model.Ap);                 % get the eigenvalues
    d = d ./ max(abs(d));               % get relative values
    zer = find(abs(d) < 1e-10);      %  the location of the 0 eigenvalue
    model.Ap(zer,zer) = 0;          % Set the values to exactly 0;
    model.An(zer,zer) = 0;
    

    model.zer = zer;
    model.nodes = xr;
end


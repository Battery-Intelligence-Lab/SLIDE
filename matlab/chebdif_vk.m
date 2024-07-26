function [x, DM] = chebdif_vk(N, M)

%  The function [x, DM] =  chebdif(N,M) computes the differentiation
%  matrices D1, D2, ..., DM on Chebyshev nodes.
%
%  Input:
%  N:        Size of differentiation matrix.
%  M:        Number of derivatives required (integer).
%  Note:     0 < M <= N-1.
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
%  The code implements two strategies for enhanced
%  accuracy suggested by W. Don and S. Solomonoff in
%  SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric
%  identities to avoid the computation of differences
%  x(k)-x(j) and (b) the use of the "flipping trick"
%  which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
%  Note added May 2003:  It may, in fact, be slightly better not to
%  implement the strategies (a) and (b).   Please consult the following
%  paper for details:   "Spectral Differencing with a Twist", by
%  R. Baltensperger and M.R. Trummer, to appear in SIAM J. Sci. Comp.

%  J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by
%  JACW, May 2003.
%  reference:
%  A MATLAB differentiation matrix suite, ACM TOMS, DOI 10.1145/365723.365727
%
%  Available at:
%  see https://uk.mathworks.com/matlabcentral/fileexchange/29-dmsuite
%  http://appliedmaths.sun.ac.za/~weideman/research/differ.html

I = eye(N);                          % Identity matrix.
L = logical(I);                      % Logical identity matrix.

n1 = floor(N/2); n2  = ceil(N/2);     % Indices used for flipping trick.

k = (0:N-1)';                        % Compute theta vector.
th = k*pi/(N-1);

x = sin(pi*(N-1:-2:1-N)'/(2*(N-1))); % Compute Chebyshev points.

DX = 2*sin((th' + th)/2).*sin((th' - th)/2);   % Trigonometric identity.

DX_vk = zeros(N,N);
C_vk = zeros(N,N);

for i=1:N
    for j = 1:N
        if(i==j)
            DX_vk(i,j) = 1;
        elseif(i <= n1)
            DX_vk(i,j) = 2*sin((th(j) + th(i))/2).*sin((th(j) - th(i))/2);
        else
            DX_vk(i,j) = -2*sin((th(N-j+1) + th(N-i+1))/2).*sin((th(N-j+1) - th(N-i+1))/2);
        end
    end
end

for i=1:N
    for j = 1:N

        C_vk(i,j) = 1;

        if(i==1 || i==N)
            C_vk(i,j) = 2*C_vk(i,j);
        end

        if(j==1 || j==N)
            C_vk(i,j) = C_vk(i,j)/2;
        end

        if(mod(i+j,2)==1)
            C_vk(i,j) = -C_vk(i,j);
        end
    end
end

D_vk = eye(N);
for ell = 1:M
    for i=1:N
        row_sum = 0;
        for j = 1:N
            if(i==j)
                continue;
            end
            D_vk(i,j) = ell*(C_vk(i,j)*D_vk(i,i) - D_vk(i,j))/DX_vk(i,j);
            row_sum  =  row_sum - D_vk(i,j);
        end
        D_vk(i,i) = row_sum;
    end
    DM_vk(:,:,ell) = D_vk;
end





DX = [DX(1:n1,:); -rot90(DX(1:n2,:),2)];   % Flipping trick.
DX(L) = 1;                       % Put 1's on the main diagonal of DX.


fprintf('Norm = %14.16f\n',  norm(DX - DX_vk));



C = toeplitz((-1).^k);               % C is the matrix with
C(1,:) = C(1,:)*2;
C(N,:) = C(N,:)*2;     % entries c(k)/c(j)

C(:,1) = C(:,1)/2;
C(:,N) = C(:,N)/2;
fprintf('Norm C_vk = %14.16f\n',  norm(C - C_vk));



Z = 1./DX;                           % Z contains entries 1/(x(k)-x(j))
Z(L) = zeros(N,1);                      % with zeros on the diagonal.

D = eye(N);                          % D contains diff. matrices.
for ell = 1:M
    D = ell*Z.*(C.*diag(D) - D); % Off-diagonals
    D(L) = -sum(D,2);                            % Correct main diagonal of D
    DM(:,:,ell) = D;                                   % Store current D in DM
end




fprintf('Norm DM_vk = %14.16f\n',  norm(DM_vk(:) - DM(:)));






end




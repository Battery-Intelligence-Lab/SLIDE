function [x, DM] = chebdif_vk2(N, M)

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


k = (0:N-1)';                        % Compute theta vector.
th = k*pi/(N-1);

x = sin(pi*(N-1:-2:1-N)'/(2*(N-1))); % Compute Chebyshev points.

D_vk = eye(N);
for ell = 1:M
    for i=1:N
        row_sum = 0;
        for j = 1:N
            if(i==j)
                continue;
            end

            if(i <= floor(N/2)) % Indices used for flipping trick.
                DX_vk = 2*sin((th(j) + th(i))/2).*sin((th(j) - th(i))/2);
            else
                DX_vk = -2*sin((th(N-j+1) + th(N-i+1))/2).*sin((th(N-j+1) - th(N-i+1))/2);
            end

            C_vk = 1;

            if(i==1 || i==N)
                C_vk = 2*C_vk;
            end
    
            if(j==1 || j==N)
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

end




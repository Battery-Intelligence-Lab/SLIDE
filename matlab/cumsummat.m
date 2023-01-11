function Q = cumsummat(N)
% CUMSUMMAT  Chebyshev integration matrix.
% Q = CUMSUMMAT(N) is the matrix that maps function values at N Chebyshev
% points to values of the integral of the interpolating polynomial at
% those points, with the convention that the first value is zero. 
% 
% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.
%
% available at:
% https://uk.mathworks.com/matlabcentral/fileexchange/23972-chebfun-v4-old-version--please-download-current-version-instead?focused=5592790&tab=function
% 
% license:
% Copyright (c) 2015, The Chancellor, Masters and Scholars of the University of Oxford, and the Chebfun Developers
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%

N = N-1;

% persistent cache    % stores computed values for fast return
% if isempty(cache), cache = {}; end    % first call

% if length(cache) >= N && ~isempty(cache{N})
%   Q = cache{N};
%   return
% else
%   cache{N} = [];
% end

% Matrix mapping coeffs -> values.
T = cp2cdm(N);

% Matrix mapping values -> coeffs.
Tinv = cd2cpm(N);

% Matrix mapping coeffs -> integral coeffs. Note that the highest order
% term is truncated. 
k = 1:N;
k2 = 2*(k-1);  k2(1) = 1;  % avoid divide by zero
B = diag(1./(2*k),-1) - diag(1./k2,1);
v = ones(N,1); v(2:2:end) = -1;
B(1,:) = sum( diag(v)*B(2:N+1,:), 1 );
B(:,1) = 2*B(:,1); 

Q = T*B*Tinv;               
% cache{N} = Q;

end

function T = cp2cdm(N)
% Values of Cheb. polys at Cheb nodes, x(n)=-cos(pi*n/N).
theta = pi*(N:-1:0)'/N;
T = cos( theta*(0:N) );
end

function C = cd2cpm(N)
% Three steps: Double the data around the circle, apply the DFT matrix,
% and then take half the result with 0.5 factor at the ends.
theta = (pi/N)*(0:2*N-1)';
F = exp( -1i*theta*(0:2*N-1) );  % DFT matrix
rows = 1:N+1;  % output upper half only
% Impose symmetries on data and coeffs.
C = real( [ F(rows,N+1) F(rows,N:-1:2)+F(rows,N+2:2*N) F(rows,1) ] );
C = C/N;  C([1 N+1],:) = 0.5*C([1 N+1],:);
end
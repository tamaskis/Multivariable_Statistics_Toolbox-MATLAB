%==========================================================================
%
% randmvu  Generate a random sample from a multivariate uniform 
% distribution.
%
%   x = randmvu(a,b)
%   x = randmvu(a,b,N)
%
% See also randmvn.
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-09-28
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   a       - (n×1 double) lower bound
%   b       - (n×n double) upper bound
%   N       - (OPTIONAL) (1×1 double) sample size (defaults to 1)
%
% -------
% OUTPUT:
% -------
%   x       - (n×N double) random sample of size N from X ~ U(a,b)
%
%==========================================================================
function x = randmvu(a,b,N)
    
    % defaults sample size to 1 if not input
    if (nargin < 3) || isempty(N)
        N = 1;
    end
    
    % dimension of random vector
    n = length(a);
    
    % preallocates array to store samples of X
    x = zeros(n,N);
    
    % generates random samples
    for k = 1:n
        x(k,:) = a(k)+(b(k)-a(k))*rand(1,N);
    end
    
end
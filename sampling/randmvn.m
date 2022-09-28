%==========================================================================
%
% randmvn Generate random samples from a multivariate normal distribution.
%
%   x = randmvn(mu,Sigma)
%   x = randmvn(mu,Sigma,N)
%
% See also randmvu.
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-09-27
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   mu      - (n×1 double) mean, μ
%   Sigma   - (n×n double) covariance, Σ
%   N       - (OPTIONAL) (1×1 double) sample size (defaults to 1)
%
% -------
% OUTPUT:
% -------
%   x       - (n×N double) N samples from X ~ N(μ,Σ)
%
% -----
% NOTE:
% -----
%   --> The covariance (Σ) must be positive semidefinite.
%
%==========================================================================
function x = randmvn(mu,Sigma,N)
    
    % defaults sample size to 1 if not input
    if (nargin < 3) || isempty(N)
        N = 1;
    end
    
    % dimension of random vector
    n = length(mu);
    
    % square root of covariance matrix (Σ¹ᐟ²) via Cholesky decomposition
    Sigma_sqrt = chol(Sigma)';
    
    % N random samples from standard normal distribution
    z = randn(n,N);
    
    % n×N matrix storing mean in each column
    mu = repmat(mu,1,N);
    
    % converts to Gaussian distribution with given parameters
    x = Sigma_sqrt*z+mu;
    
end
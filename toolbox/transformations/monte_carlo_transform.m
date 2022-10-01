%==========================================================================
%
% monte_carlo_transform  Monte Carlo transformation for passing a Gaussian 
% through a nonlinearity.
%
%   [mu_y,Sigma_yy] = monte_carlo_transform(mu_x,Sigma_xx,f)
%   [mu_y,Sigma_yy] = monte_carlo_transform(mu_x,Sigma_xx,f,N)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-09-27
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% TOOLBOX DOCUMENTATION:
% https://tamaskis.github.io/Multivariable_Statistics_Toolbox-MATLAB/
%
% TECHNICAL DOCUMENTATION:
% https://tamaskis.github.io/files/Multivariable_Statistics_Tools.pdf
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   mu_x        - (n×1 double) mean of X
%   Sigma_xx    - (n×n double) covariance of X
%   f           - (1×1 function_handle) multivariate, vector-valued 
%                 function, Y = f(X) (f : ℝⁿ → ℝᵐ)
%   N           - (OPTIONAL) (1×1 double) sample size (defaults to 1000)
%
% -------
% OUTPUT:
% -------
%   mu_y        - (m×1 double) mean of Y
%   Sigma_yy    - (m×m double) covariance of Y
%
% -----
% NOTE:
% -----
%   --> This implementation assumes that X is a Gaussian random vector.
%
%==========================================================================
function [mu_y,Sigma_yy] = monte_carlo_transform(mu_x,Sigma_xx,f,N)
    
    % defaults sample size to 1000 if not specified
    if nargin < 4
        N = 1000;
    end
    
    % generate random samples of X
    x = randmvn(mu_x,Sigma_xx,N);
    
    % dimension of Y
    m = length(f(mu_x));
    
    % preallocates matrix to store transformations of X
    y = zeros(m,N);
    
    % transforms X to Y
    for k = 1:N
        y(:,k) = f(x(:,k));
    end
    
    % approximates statistics of Y using sample statistics
    mu_y = sample_mean(y);
    Sigma_yy = sample_covariance(y);
    
end
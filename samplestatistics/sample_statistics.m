%==========================================================================
%
% sample_statistics  Mean and covariance of a sample.
%
%   [mu,Sigma] = sample_statistics(x)
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
%   x       - (n×N double) N samples of X
%
% -------
% OUTPUT:
% -------
%   mu      - (n×1 double) mean of X, μ
%   Sigma   - (n×n double) covariance of X, Σ
%
%==========================================================================
function [mu,Sigma] = sample_statistics(x)
    
    % dimension of X
    n = size(x,1);
    
    % sample size
    N = size(x,2);
    
    % mean
    mu = zeros(n,1);
    for i = 1:N
        mu = mu+x(:,i);
    end
    mu = mu/N;
    
    % covariance
    Sigma = zeros(n,n);
    for i = 1:N
        Sigma = Sigma+(x(:,i)-mu)*(x(:,i)-mu).';
    end
    Sigma = Sigma/(N-1);
    
end
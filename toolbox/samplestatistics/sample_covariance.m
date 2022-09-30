%==========================================================================
%
% sample_covariance  Covariance of a sample.
%
%   Qxx = sample_covariance(x)
%   Qxx = sample_covariance(x,w)
%   Qxx = sample_covariance(__,bessel)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-09-29
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
%   x       - (n×N double) sample of X (sample size = N)
%   w       - (OPTIONAL) (N×1 double) weight vector (defaults to vector of
%             1/N's)
%   bessel  - (OPTIONAL) (n×1 double) true if Bessel's correction should be
%             included, false otherwise (defaults to true)
%
% -------
% OUTPUT:
% -------
%   Qxx     - (n×n double) sample covariance
%
% -----
% NOTE:
% -----
%   --> The weight vector, w, does not need to sum to 1. If the weight 
%       vector does not sum to 1, it is automatically normalized to do so.
%
%==========================================================================
function Qxx = sample_covariance(x,w)
    
    % sample size
    N = size(x,2);
    
    % default weight vector to 1/N in each element if not input
    if (nargin == 1) || isempty(w)
        w = (1/N)*ones(N,1);
    end
    
    % determine the β coefficient
    if bessel
        beta = N/(N-1);
    else
        beta = 1;
    end
    
    % vector of ones
    ones_vec = ones(N,1);
    
    % auxiliary matrices
    A = x-((x*w)/(w.'*ones_vec))*ones_vec.';
    W = diag(w);
    
    % sample covariance
    Qxx = (beta*A*W*A.')/(w.'*ones_vec);
    
end
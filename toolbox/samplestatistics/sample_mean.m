%==========================================================================
%
% sample_mean  Mean of a sample.
%
%   xbar = sample_mean(x)
%   xbar = sample_mean(x,w)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-10-01
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
%
% -------
% OUTPUT:
% -------
%   xbar    - (n×1 double) sample mean
%
% -----
% NOTE:
% -----
%   --> The weight vector, w, does not need to sum to 1. If the weight 
%       vector does not sum to 1, it is automatically normalized to do so.
%
%==========================================================================
function xbar = sample_mean(x,w)
    
    % sample size
    N = size(x,2);
    
    % default weight vector to 1/N in each element if not input
    if (nargin == 1) || isempty(w)
        w = (1/N)*ones(N,1);
    end
    
    % sample mean
    xbar = x*w/(w.'*ones(N,1));
    
end
%==========================================================================
%
% sample_cross_covariance  Cross covariance of two samples.
%
%   Qxy = sample_covariance(x,y)
%   Qxy = sample_covariance(x,y,w)
%   Qxy = sample_covariance(x,y,wm,wc)
%   Qxy = sample_covariance(__,bessel)
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
%   x       - (n×N double) sample of X ∈ ℝⁿ (sample size = N)
%   y       - (m×N double) sample of Y ∈ ℝᵐ (sample size = N)
%   wm      - (OPTIONAL) (N×1 double) mean weight vector (defaults to 
%             vector of 1/N's)
%   wc      - (OPTIONAL) (N×1 double) covariance weight vector (defaults to 
%             wm)
%   bessel  - (OPTIONAL) (n×1 double) true if Bessel's correction should be
%             included, false otherwise (defaults to true)
%
% -------
% OUTPUT:
% -------
%   Qxy     - (n×m double) sample cross covariance
%
% -----
% NOTE:
% -----
%   --> The mean weight vector, wm, does not need to sum to 1; if it
%       doesn't already, it is normalized to do so.
%   --> The covariance weight vector, wc, does not need to sum to 1; if it
%       doesn't already, it is normalized to do so.
%
%==========================================================================
function Qxy = sample_cross_covariance(x,y,wm,wc,bessel)
    
    % include Bessel's correction by default unless otherwise specified
    if (nargin < 3) || isempty(bessel)
        bessel = true;
    end
    
    % sample size
    N = size(x,2);
    
    % default mean weight vector to vector of 1/N's if not input
    if (nargin < 3) || isempty(wm)
        wm = (1/N)*ones(N,1);
    end
    
    % default covariance weight vector to mean weight vector if not input
    if (nargin < 4) || isempty(wc)
        wc = wm;
    end
    
    % determine the β coefficient
    if bessel
        beta = N/(N-1);
    else
        beta = 1;
    end
    
    % vector of ones
    ones_vec = ones(N,1);
    
    % mean difference matrices
    Mx = x-((x*wm)/(wm.'*ones_vec))*ones_vec.';
    My = y-((y*wm)/(wm.'*ones_vec))*ones_vec.';
    
    % weight matrix
    W = diag(wc);
    
    % sample cross covariance
    Qxy = (beta*Mx*W*My.')/(wc.'*ones_vec);
    
end
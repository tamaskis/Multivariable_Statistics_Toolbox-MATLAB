%==========================================================================
%
% sample_covariance  Covariance of a sample.
%
%   Qxx = sample_covariance(x)
%   Qxx = sample_covariance(x,w)
%   Qxx = sample_covariance(x,wm,wc)
%   Qxx = sample_covariance(__,bessel)
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
%   Qxx     - (n×n double) sample covariance
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
function Qxx = sample_covariance(x,wm,wc,bessel)
    
    % include Bessel's correction by default unless otherwise specified
    if (nargin < 3) || isempty(bessel)
        bessel = true;
    end
    
    % sample size
    N = size(x,2);
    
    % default mean weight vector to vector of 1/N's if not input
    if (nargin < 2) || isempty(wm)
        wm = (1/N)*ones(N,1);
    end
    
    % default covariance weight vector to mean weight vector if not input
    if (nargin < 3) || isempty(wc)
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
    
    % mean difference matrix
    M = x-((x*wm)/(wm.'*ones_vec))*ones_vec.';
    
    % weight matrix
    W = diag(wc);
    
    % sample covariance
    Qxx = (beta*M*W*M.')/(wc.'*ones_vec);
    
end
%==========================================================================
%
% samples_full_factorial  Full factorial sampling plan.
%
%   X = samples_full_factorial(a,b,m)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-09-27
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% REFERENCES:
%   [1] Kochenderfer and Wheeler, "Algorithms for Optimization"
%       (pp. 235-236)
%   [2] https://www.mathworks.com/matlabcentral/answers/263932-unknown-number-of-output-variables
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   a       - (n×1 double) lower bound
%   b       - (n×1 double) upper bound
%   m       - (n×1 double) number of samples in each direction
%
% -------
% OUTPUT:
% -------
%   x       - (n×N double) N samples of x
%
% -----
% NOTE:
% -----
%   --> M = Πᵢmᵢ (i.e. product of elements of m)
%
%==========================================================================
function X = samples_full_factorial(a,b,m)
    
    % dimension of x
    n = length(a);
    
    % 1-dimensional grids for each design variable
    x_1d = cell(1,n);
    for i = 1:n
        x_1d{i} = a(i):((b(i)-a(i))/(m(i)-1)):b(i);
    end
    
    % n-dimensional grids for each design variable
    x_nd = cell(1,n);
    [x_nd{:}] = ndgrid(x_1d{:});
    
    % defines an n×N matrix storing N samples of x ∈ ℝⁿ
    X = zeros(n,length(x_nd{i}(:)));
    for i = 1:n
        X(i,:) = x_nd{i}(:);
    end
    
end
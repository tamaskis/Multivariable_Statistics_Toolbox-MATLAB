%==========================================================================
%
% eval_surrogate  Evaluates a surrogate model defined using radial basis 
% functions.
%
%   f = eval_surrogate(theta,x,c,psi)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-09-27
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% REFERENCES:
%   [1] Kochenderfer and Wheeler, "Algorithms for Optimization"
%       (pp. 255, 259, 262-263)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   theta   - (m×1 double) surrogate model coefficients
%   x       - (n×1 double) evaluation point
%   c       - (n×m double) centers for radial basis functions
%   psi     - (1×1 function_handle) kernel for radial basis functions
%
% -------
% OUTPUT:
% -------
%   f       - (1×1 double) evaluation of surrogate model
%
%==========================================================================
function f = eval_surrogate(theta,x,c,psi)
    
    % evaluates radial basis functions
    b = radial_basis(x,c,psi);
    
    % evaluates surrogate model
    f = (theta.')*b;
    
end
%==========================================================================
%
% radial_basis  Evaluates m radial basis functions.
%
%   b = radial_basis(x,c,psi)
%
% Copyright © 2022 Tamas Kis
% Last Update: 2022-09-27
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% REFERENCES:
%   [1] Kochenderfer and Wheeler, "Algorithms for Optimization"
%       (pp. 255, 259, 262-263)
%   [2] https://en.wikipedia.org/wiki/Radial_basis_function
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   x       - (n×1 double) evaluation point
%   c       - (n×N double) centers of radial basis functions
%   psi     - (1×1 function_handle) kernel
%
% -------
% OUTPUT:
% -------
%   b       - (N×1 double) evaluation of radial basis functions at x
%
%==========================================================================
function b = radial_basis(x,c,psi)
    
    % number of radial basis functions
    N = size(c,2);
    
    % preallocates vector b to store evaluations of radial basis functions
    b = zeros(N,1);
    
    % evaluates radial basis functions
    for i = 1:N
        b(i) = psi(sqrt((x-c(:,i)).'*(x-c(:,i))));
    end
    
end
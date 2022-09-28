%==========================================================================
%
% get_surrogate  Develops a surrogate model using radial basis functions.
%
%   f_hat = get_surrogate(f,x)
%   f_hat = get_surrogate(f,x,type)
%   f_hat = get_surrogate(f,x,type,sigma)
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
%   f       - (1×1 function_handle) scalar-valued, multivarite function, 
%             f(x) (f : ℝⁿ → ℝ)
%   X       - (n×N double) N samples of x
%   type    - (OPTIONAL) (char) 'linear', 'cubic', 'thin plate spline', 
%             'Gaussian', 'multiquadratic', or 'inverse multiquadratic'
%             (defaults to 'cubic')
%   sigma   - (OPTIONAL) (1×1 double) hyperparameter (needed for 
%             'Gaussian', 'multiquadratic', or 'inverse multiquadratic'
%             kernels) (defaults to 1)
%
% -------
% OUTPUT:
% -------
%   f_hat   - (1×1 function_handle) surrogate model, f_hat(x) 
%             (f_hat : ℝⁿ → ℝ)
%
%==========================================================================
function f_hat = get_surrogate(f,X,type,sigma)
    
    % defaults "type" to 'cubic'
    if (nargin < 3) || isempty(type)
        type = 'cubic';
    end
    
    % defaults "sigma" to 1
    if (nargin < 4) || isempty(sigma)
        sigma = 1;
    end
    
    % kernel for radial basis functions
    psi = kernel(type,sigma);
    
    % sample size
    N = size(X,2);
    
    % sample points serve as centers for radial basis functions
    c = X;
    
    % preallocates design matrix
    B = zeros(N,N);
    
    % populates design matrix
    for i = 1:N
        B(i,:) = radial_basis(X(:,i),c,psi).';
    end
    
    % evaluates function at each sample point
    y = zeros(N,1);
    for i = 1:N
        y(i) = f(X(:,i));
    end
    
    % solves for surrogate model coefficients
    theta = B\y;
    
    % surrogate model
    f_hat = @(x) (theta.')*radial_basis(x,c,psi);
    
end
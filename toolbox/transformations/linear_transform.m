%==========================================================================
%
% linear_transform  Linear transformation of mean and covariance.
%
%   [mu_y,Sigma_yy] = linear_transform(mu_x,Sigma_xx,A)
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
%   A           - (m×n double) matrix defining linear transformation Y = AX
%
% -------
% OUTPUT:
% -------
%   mu_y        - (m×1 double) mean of Y
%   Sigma_yy    - (m×m double) covariance of Y
%
%==========================================================================
function [mu_y,Sigma_yy] = linear_transform(mu_x,Sigma_xx,A)
    
    % mean of Y
    mu_y = A*mu_x;
    
    % covariance of Y
    Sigma_yy = A*Sigma_xx*A.';
    
end
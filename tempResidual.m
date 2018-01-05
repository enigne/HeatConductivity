% Construct F(x) for the least-square problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'K'             - heat conductivity for the model;
%   'HeatSolver'    - function handler for the heat equation;
%   'z'             - the z-coordinates of the numerical model;
%   't'             - the time of the numerical model;
%   'f_data'        - measured data projected to the computational nodes;
%   'heatParam'     - other coefficients for the heat equation:      
% The return values:
%   'F'             - minimize f(x) = F(x)^T*F(x), F(x) is the mismatch
%                     between the model and the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = tempResidual(K, HeatSolver, z, t, f_data, heatParam)
    T_sol = HeatSolver(t, z, K, heatParam);
   
    err = T_sol - f_data;   
    mask = (f_data<-2);
    
    compareErr = err(mask);
    F = compareErr(:);
end

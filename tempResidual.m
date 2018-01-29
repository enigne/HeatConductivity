% Construct F(x) for the least-square problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'K'             - heat conductivity for the model;
%   'HeatSolver'    - function handler for the heat equation;
%   't'             - the time of the numerical model;
%   'z'             - the z-coordinates of the numerical model;
%   'f_data'        - measured data projected to the computational nodes;
%   'w'             - the weights get from STD;
%   'heatParam'     - other coefficients for the heat equation:      
% The return values:
%   'F'             - minimize f(x) = F(x)^T*F(x), F(x) is the mismatch
%                     between the model and the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = tempResidual(K, HeatSolver, t, z, T_data, t_data, z_data, w, heatParam)
    % Project data to the computational domain
%     f_data = project2D(T_data, t_data, z_data, t, z);

    T_sol = HeatSolver(t, z, K, heatParam);
   
    T_sol_int = project2D(T_sol, t, z, t_data, z_data);
    err = w .* (T_sol_int - T_data);   
    
    mask = (((T_data<-2) & (~isnan(err))));
    
    compareErr = err(mask);
    F = compareErr(:);
end

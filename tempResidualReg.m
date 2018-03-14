% Construct F(x) for the least-square problem with a regularization term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'K'             - heat conductivity for the model;
%   'HeatSolver'    - function handler for the heat equation;
%   't'             - the time of the numerical model;
%   'z'             - the z-coordinates of the numerical model;
%   'f_data'        - measured data projected to the computational nodes;
%   'w'             - the weights get from STD;
%   'heatParam'     - other coefficients for the heat equation,      
%   'gamma'         - regularization parameter.     
% The return values:
%   'F'             - minimize f(x) = F(x)^T*F(x), F(x) is the mismatch
%                     between the model and the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = tempResidualReg(x, HeatSolver, t, z, T_data, t_data, z_data, w, heatParam, gamma)
    % solve heat equation with given K and rho in x
    T_sol = HeatSolver(t, z, x, heatParam);
   
    % Project solutions to the same mesh as the measurements
    T_sol_int = project2D(T_sol, t, z, t_data, z_data);
    err = w .* (T_sol_int - T_data);   
    
    % create mask for temperature lower than maskConst
    mask = (((T_data < heatParam.maskConst) & (~isnan(err))));
    
    % Apply mask
    compareErr = err(mask);
    n = length(compareErr);
    
    % Construct regularization for rho
    Nk = length(heatParam.zK);
    rho.rho = x(Nk+1:end) * heatParam.rhoScale;
    rho.z = heatParam.zK;
    rhoP = density(rho, heatParam.rho.z);
    reg = rhoP(:) - heatParam.rho.rho(:);
    
    % Add regularization
    F = [1./sqrt(n) .*compareErr(:); gamma*reg];
end

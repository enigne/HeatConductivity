% To compute the sensitivity matrix for k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'z_data'        - the z-coordinates of the measured data;
%   't_data'        - the time of the data;
%   'T_data'        - measured data(temperature);
%   'dZfine'        - refinement factor for the computational nodes;
%   'zK'            - the z-coordinates of K;
%   'K0'            - the deepth dependent heat conductivity;
%   'dK'            - delta k for the numerical solution of the gradient;
%   'rho'           - density;
%   'C'             - heat capacity;
%   'interpOption'  - intepolation method (linear by default).
% The return values:
%   'dTdz'          - The derivatives of T with respect to z at t-z plan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A] = computeKSensitivity(z_data, t_data, T_data, dZfine, zK, K0, dK, rho, C, interpOption)
    Nz = length(z_data);
    Nzfine = dZfine * (Nz - 1) + 1;
    xInd = [1:dZfine:Nzfine]';

    % Time discretization same size as measurment
    Nt = length(t_data);
    
    % Set initial and boundary conditions
    [Tbc, T0, z, t, dz, dt] = setIBCs(z_data, t_data, Nzfine, Nt, T_data, interpOption);
    
    % Set Parameters for solving
    heatParam = setHeatParam(dt, Nt, dz, Nzfine, rho, C, T0, Tbc.Up, Tbc.Down, zK);
    
    %% Solve for the optimal solution
    T0_sol = solveHeat(t, z, K0, heatParam);
    
    %%
    Nk = length(K0);
    
    dKvec = dK * speye(Nk);
    A = zeros(Nz * Nt, Nk);
    
    for i = 1 : Nk
        T_sol_temp = solveHeat(t, z, K0+dKvec(:,i), heatParam);
        dT_temp = (T0_sol - T_sol_temp) /dK;
        dT_temp = dT_temp(xInd, :);
        A(:, i) = dT_temp(:);
    end
end
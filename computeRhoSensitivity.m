% To compute the sensitivity matrix for rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'z_data'        - the z-coordinates of the measured data;
%   't_data'        - the time of the data;
%   'T_data'        - measured data(temperature);
%   'dZfine'        - refinement factor for the computational nodes;
%   'zK'            - the z-coordinates of K;
%   'K0'            - the depth dependent heat conductivity;
%   'dRho'          - delta rho for the numerical solution of the gradient;
%   'rho'           - density;
%   'zRho'          - the depth of rho;
%   'C'             - heat capacity;
%   'interpOption'  - intepolation method (linear by default).
% The return values:
%   'dTdz'          - The derivatives of T with respect to z at t-z plan.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A] = computeRhoSensitivity(z_data, t_data, T_data, dZfine, zK, K0, dRho, rho, zRho, C, interpOption)
    Nz = length(z_data);
    Nzfine = dZfine * (Nz - 1) + 1;
    xInd = [1:dZfine:Nzfine]';

    % Time discretization same size as measurment
    Nt = length(t_data);
    
    % Set initial and boundary conditions
    [Tbc, T0, z, t, dz, dt] = setIBCs(z_data, t_data, Nzfine, Nt, T_data, interpOption);
    
    % Set Parameters for solving
    heatParam = setHeatParam(dt, Nt, dz, Nzfine, rho, C, T0, Tbc.Up, Tbc.Down, zK, zRho);
    
    %% Solve for the optimal solution
    T0_sol = solveHeat(t, z, K0, heatParam);
    
    %%
    NRho = length(rho);
    
    dRhovec = dRho * speye(NRho);
    A = zeros(Nz * Nt, NRho);
    
    heatParam_temp = heatParam;
    
    for i = 1 : NRho
        heatParam_temp.rho = heatParam.rho + dRhovec(:, i);
        T_sol_temp = solveHeat(t, z, K0, heatParam_temp);
        dT_temp = (T0_sol - T_sol_temp) / dRho;
        dT_temp = dT_temp(xInd, :);
        A(:, i) = dT_temp(:);
    end
end
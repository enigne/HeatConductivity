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
%   'rho'           - density data including deepth(rho.z) and the density(rho.rho);
%   'C'             - heat capacity;
%   'interpOption'  - intepolation method (linear by default).
% The return values:
%   'A'             - The snesivitity matrix of rho.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A] = computeRhoSensitivity(z_data, t_data, T_data, dZfine, zK, K0, dRho, rho, C, mask, interpOption)
    Nz = length(z_data);
    Nzfine = dZfine * (Nz - 1) + 1;
    xInd = [1:dZfine:Nzfine]';
    
    % Set initial and boundary conditions
    [Tbc, T0, z, t, dz, Nt, dt] = setIBCs(z_data, t_data, Nzfine, T_data, interpOption);
    
    % Set Parameters for solving
    heatParam = setHeatParam(dt, Nt, dz, Nzfine, rho, C, T0, Tbc.Up, Tbc.Down, zK);
    
    %% Solve for the optimal solution
    T0_sol = solveHeat(t, z, K0, heatParam);
    
    %%
    NRho = length(rho.rho);
    Nt_data = length(t_data);

    dRhovec = dRho * speye(NRho);
    A = zeros(Nz * Nt_data, NRho);
    
    heatParam_temp = heatParam;
    
    for i = 1 : NRho
        heatParam_temp.rho.rho = heatParam.rho.rho + dRhovec(:, i);
        T_sol_temp = solveHeat(t, z, K0, heatParam_temp);
        dT_temp = (T_sol_temp - T0_sol) / dRho;

        dT_temp = project2D(dT_temp, t, z, t_data, z_data);
        A(:, i) = dT_temp(:);
    end
    A(mask, :) = 0;
end
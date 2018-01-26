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
%   'A'             - The sensitivity matrix of K.
%   'z'             - The z-coordinates for the computation
%   't'             - The time spots for the computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, z, t] = computeKSensitivity(z_data, t_data, T_data, dZfine, zK, K0, dK, rho, C, mask, interpOption)
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
    Nk = length(K0);
    Nt_data = length(t_data);
    
    dKvec = dK * speye(Nk);
    A = zeros(Nz * Nt_data, Nk);
    
    for i = 1 : Nk
        T_sol_temp = solveHeat(t, z, K0+dKvec(:,i), heatParam);
        dT_temp = ( T_sol_temp - T0_sol ) /dK;
            
        dT_temp = project2D(dT_temp, t, z, t_data, z_data);

        A(:, i) = dT_temp(:);
    end

    % Set the area with mask to 0
    A(mask, :) = 0;
    
    % Scale the unit in time from days to seconds
    A = A .* (24*3600);
end
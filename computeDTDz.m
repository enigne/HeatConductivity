% To compute dT/dz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'z_data'        - the z-coordinates of the measured data;
%   't_data'        - the time of the data;
%   'T_data'        - measured data(temperature);
%   'dZfine'        - refinement factor for the computational nodes;
%   'zK'            - the z-coordinates of K;
%   'K0'            - the deepth dependent heat conductivity;
%   'rho'           - density;
%   'C'             - heat capacity;
%   'interpOption'  - intepolation method (linear by default).
% The return values:
%   'dTdz'          - The derivatives of T with respect to z at t-z plan.
%   'matDTDz'       - The derivatives in the matrix form: dT=matDTdz*dz.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dTdz, matDTDz] = computeDTDz(z_data, t_data, T_data, dZfine, zK, K0, rho, C, interpOption)
    % Number of spatial discretization
    Nz = dZfine * (length(z_data)-1) + 1;

    % Time discretization same size as measurment
    Nt = length(t_data);

    % Set initial and boundary conditions
    [Tbc, T0, z, t, dz, dt] = setIBCs(z_data, t_data, Nz, Nt, T_data, interpOption);

    % Set Parameters for solving
    heatParam = setHeatParam(dt, Nt, dz, Nz, rho, C, T0, Tbc.Up, Tbc.Down, zK);

    %% Compute dT/dz
    % Solve heat equation on a finear mesh
    T0_sol = solveHeat(t, z, K0, heatParam);

    % Shift the indices for finite differences
    xInd = 1:dZfine:Nz;
    xIndP1 = xInd + 1;
    xIndP1(end) = Nz;

    Tcenter = T0_sol(xInd, :);
    Tup = T0_sol(xIndP1, :);

    % one-side finite differences
    dTdz = 1.0 / dz .* (Tup - Tcenter);
    
    % convert to matrix-vecotr multiplication form
    nDz = length(z_data);
    eyeDz = speye(nDz);
    matDTDz = sparse(repmat(eyeDz, Nt, 1)) .* dTdz(:);
end
% To compute dT/dz only from the measurments
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
%   'z'             - The z-coordinates for the computation
%   't'             - The time spots for the computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dTdz, matDTDz, z, t] = computeDTDzFromData(z_data, t_data, T_data, dZfine, mask, interpOption)
    % Number of spatial discretization
    Nz = dZfine * (length(z_data)-1) + 1;

    % Set initial and boundary conditions
    [~, ~, z, t, dz, ~, ~] = setIBCs(z_data, t_data, Nz, T_data, interpOption);

    % project T_data on the finer mesh
    T_data_finer = project2D(T_data, t_data, z_data, t, z);

    % Shift the indices for finite differences
    xInd = 1:dZfine:Nz;
    xIndP1 = xInd + 1;
    xIndP1(end) = Nz;

    Tcenter = T_data_finer(xInd, :);
    Tup = T_data_finer(xIndP1, :);

    % one-side finite differences
    dTdz = 1.0 / dz .* (Tup - Tcenter);
    
    % Project to the data points 
    dTdz = project2D(dTdz, t, z(xInd), t_data, z_data);

    % Set the area with mask to 0
    dTdz(mask) = 0;
    
    % convert to matrix-vecotr multiplication form
    nDz = length(z_data);
    nTdata = length(t_data);
    eyeDz = speye(nDz);
    
    matDTDz = sparse(repmat(eyeDz, nTdata, 1)) .* dTdz(:);        
end
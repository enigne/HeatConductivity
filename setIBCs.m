% Set initial and boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'z_data'        - the z-coordinates of the measured data;
%   't_data'        - the time of the data;
%   'Nz'            - number of spatial discretization;
%   'T_data'        - measured data(temperature);
%   'interpOption'  - intepolation method (linear by default);
%   'zRange'        - cut-off range for z-coordinates
% The return values:
%   'Tbc'           - Dirichlet boundary condition for the given Nt grid;
%   'T0'            - intial condition for the given Nz grid;
%   'z'             - grid point of the spatial discretization;
%   't'             - grid point of the time discretization;       
%   'dz'            - step size of the spatial discretization;
%   'Nt'            - number of time discretization;
%   'dt'            - step size of the time discretization;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tbc, T0, z, t, dz, Nt, dt] = setIBCs(z_data, t_data, Nz, T_data, interpOption, zRange)
    if nargin < 6
        % Spatial grid
        z = linspace(min(z_data), max(z_data), Nz)';
    else
        % take the minimal range between z_data and zRange
        zmin = max(min(z_data), zRange(1));
        zmax = min(max(z_data), zRange(2));
        z = linspace(zmin, zmax, Nz)';
    end
    dz = abs(z(2) - z(1));

    % Interpolate initial condition
    T0 = interp1(z_data, T_data(:,1), z, interpOption);

    % Time discretization 
    dt = min(diff(t_data));
    t = min(t_data) : dt : max(t_data);
    Nt = length(t);

    % Interpolate Boundary conditions
    Tbc.Up = interp1(t_data, T_data(1,:), t, interpOption);
    Tbc.Down = interp1(t_data, T_data(end,:), t, interpOption);
end
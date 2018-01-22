% Solve the least-square problem to optimize K for the Heat equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'data'      - the whole data set as a structure, with:
%               	'data.t'        - the time vector of the measurements;
%                   'data.z'        - the depth of the measurements;
%                   'data.T'        - the temperature measurements;
%                   'data.T_corr'   - the temperature after correction;
%                   'data.z_corr'   - the depth after correction;
%                   'data.T_i'      - the interpolated temperature;
%                   'data.z_i'      - the depth for the interpolation;
%                   'data.z_a'      - the averaged depth;
%                   'data.T_a'      - the averaged interpolated temperature;
%                   'data.T_sd'     - the standard deviation of the
%                                     temperature over the 9 holes
%   'dataIndex'	- index of the holes to be used, 0 means to use the average
%                 date
%   'zK'    	- z-coordinate of the K parameter
%   'K0'     	- the initial guess of K;
%   'Nz'        - number of grid for the computation;
%   'rho'       - density of the ice;
%   'noise'     - artificial noise.
% The return values:
%   'K_opt'     - the optimal solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K_opt = inverseK(data, dataIndex, zK, K0, Nz, rho, noise)
    % Check the input variables
    if nargin < 7
        noise = 0;
        if nargin < 6
            rho = 900;
        end
    end
    %% Initialize

    % Settings
    interpOption = 'linear';
    
    % Measurements
    if (dataIndex > 0) && (dataIndex < length(data.T))
        % Time vector (row)
        t_data = data.t';
        % Position vector (column)
        z_data = data.z_i{dataIndex};
        % Temperature
        T_data = data.T_i{dataIndex};
    else
        % Take the average value
        % Time vector (row)
        t_data = data.t';
        % Position vector (column)
        z_data = data.z_a;
        % Temperature
        T_data = data.T_a;        
    end
    
    % Add noise 
    if noise > 0 
        noiseT = noise*(1 - 2*rand(size(T_data)));
        T_data = T_data + noiseT;
    end
    
    % Physical parameters
    C = 152.5 + 7.122 * (273.15 - 10);

    %% Solve Heat equation
  
    % cut the data according to the range of K
    [T_data, z_data] = cutData(T_data, z_data, [zK(1),zK(end)]);

    % Set initial and boundary conditions
    [Tbc, T0, z, t, dz, Nt, dt] = setIBCs(z_data, t_data, Nz, T_data, interpOption);

    % Set Parameters for solving
    heatParam = setHeatParam(dt, Nt, dz, Nz, rho, C, T0, Tbc.Up, Tbc.Down, zK);

    % Project data to the computational domain
    f_data = project2D(T_data, t_data, z_data, t, z);

    %% Optimisation
    % Create the objective function
    objF = @(K) tempResidual(K, @solveHeat, z, t, f_data, heatParam);
    % Set options
    options = optimoptions('lsqnonlin','Display','iter','typicalX', K0,'TolFun', 1e-8, 'TolX', 1e-8);
    % Solve
    [K_opt,resnorm,residual,exitflag,output] = lsqnonlin(objF, K0, [], [],options);
end
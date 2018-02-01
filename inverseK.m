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
%   'noise'     - artificial noise matrix.
% The return values:
%   'K_opt'     - the optimal solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K_opt = inverseK(data, dataIndex, zK, K0, Nz, rho, w, noise)
    % Check the input variables
    if nargin < 8
        noise = [];
        if nargin < 7
            w = 1;
            if nargin < 6
                rho = 900;
            end
        end
    end
    %% Initialize

    % Settings
    interpOption = 'linear';
    
    % Load Measurements
    [t_data, z_data, T_data] = loadData(data, dataIndex);
    
    % Add noise 
    if ~isempty(noise)
        noiseT = zeros(size(T_data));
        noiseT(noise.zInd, noise.tInd) = noise.noise;
        T_data = T_data + noiseT;
    end
    
    % Physical parameters
    C = 152.5 + 7.122 * (273.15 - 10);

    %% Solve Heat equation
  
    % cut the data according to the range of K
    [T_data, z_data, ~] = cutData(T_data, z_data, [zK(1),zK(end)]);

    % Set initial and boundary conditions
    [Tbc, T0, z, t, dz, Nt, dt] = setIBCs(z_data, t_data, Nz, T_data, interpOption);

    % Set Parameters for solving
    heatParam = setHeatParam(dt, Nt, dz, Nz, rho, C, T0, Tbc.Up, Tbc.Down, zK);

    %% Optimisation
%     T_data(2:end-1,2:end-1)= T_data(2:end-1,2:end-1) + 0.1;

    % Create the objective function
    objF = @(K) tempResidual(K, @solveHeat, t, z, T_data, t_data, z_data, w, heatParam);
    % Set options
    options = optimoptions('lsqnonlin','Display','iter', 'algorithm', 'levenberg-marquardt', ...
        'typicalX', K0,'TolFun', 1e-10, 'TolX', 1e-10);
    % Solve
    [K_opt,resnorm,residual,exitflag,output] = lsqnonlin(objF, K0, [], [],options);
end
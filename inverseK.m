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
%   'zK'    	- z-coordinate of the K parameter;
%   'x0'     	- the initial guess of the unknown;
%   'Nz'        - number of grid for the computation;
%   'rho'       - density of the ice;
%   'noise'     - artificial noise matrix (not in used);
%   'timePeriod'- only part of t_data are taken into account. [0,1] indicates
%                 the full data piece;
%   'includeRho'- flag to optimize rho with a regularization term.
% The return values:
%   'x_opt'     - the optimal solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_opt = inverseK(data, dataIndex, zK, x0, Nz, rho, w, noise, timePeriod, includeRho, gamma)
    % Check the input variables
    if nargin < 11
        gamma = 1e-2;
        if nargin < 10
            includeRho = 0;
            if nargin < 9
                timePeriod = [0, 1];
                if nargin < 8
                    noise = [];
                    if nargin < 7
                        w = 1;
                        if nargin < 6
                            rho = 900;
                        end
                    end
                end
            end
        end
    end
    %% Initialize

    % Settings
    interpOption = 'linear';
    
    % Load Measurements
    [t_data, z_data, T_data] = loadData(data, dataIndex, timePeriod);
    
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
    heatParam = setHeatParam(dt, Nt, dz, Nz, rho, C, T0, Tbc, zK);

    %% Optimisation
    % Create the objective function
    if includeRho
        % optimize with Rho
        objF = @(x) tempResidualReg(x, @solveHeat, t, z, T_data, t_data, z_data, w, heatParam, gamma);
    else
        %without Rho
        objF = @(x) tempResidual(x, @solveHeat, t, z, T_data, t_data, z_data, w, heatParam);
    end
    
    % Set options
    options = optimoptions('lsqnonlin','Display','iter', ...
        'typicalX', x0,'TolFun', 1e-10, 'TolX', 1e-10, 'MaxFunEvals', 10000);
    %'algorithm', 'levenberg-marquardt', 
    
    % Solve
    [x_opt,resnorm,residual,exitflag,output] = lsqnonlin(objF, x0, zeros(size(x0)), 3*ones(size(x0)), options);
end
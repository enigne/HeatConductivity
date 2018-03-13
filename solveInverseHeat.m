% Test with solving inverse heat equation to find the optimal K and use the
% solution to solve forwardly, and generate the solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'yearIndex'	- index of the year to be used, from 2012-2015
%   'dataIndex'	- index of the holes to be used, 0 means to use the average
%                 date
%   'zK'    	- z-coordinate of the K parameter
%   'timePeriod'- only part of t_data are taken into account. [0,1] indicates
%                 the full data piece;
%   'perturbT'  - artificial noise matrix.
% The return values:
%   'K_opt'     - the optimal solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K_opt_out, t_data] = solveInverseHeat(yearIndex, dataIndex, zK, timePeriod, perturbT)
    %% Check the input
    if nargin < 5
        % perturbation test on the data
        perturbT = [];
        if nargin < 4
            % use part of the time series of the whole data
            timePeriod = [0, 1];
            if nargin < 3
                % solve K on zK only
                zKmasked = linspace(1, 8, 5)';
                if nargin < 2
                    % use the averaged data
                    dataIndex = 0;
                end
            end
        end
    end
    %% Initialize
    % Settings
    interpOption = 'linear';

    %% Load data
    try
        load('LF_4_aver.mat');
    catch
        error('Please check the original data file LF_4_aver.mat');
    end
    
    try
        load('densityData.mat');
    catch
        error('densityData.mat not found. Try to run preprocessRho.m first.');
    end
    
    try
        load('summary.mat');
    catch
        error('summary.mat not found. Try to run averageAndVariance.m first.');
    end

    % Assign data
    data = LF{yearIndex}.T;
    rho = rhoData{yearIndex};

    % Load Measurements
    [t_data, z_data, T_data, indt] = loadData(data, dataIndex, timePeriod);
    
    % Physical parameters
    C = 152.5 + 7.122 * (273.15 - 10);

    %% Initialize
    [umaskedZ, ~] = findUnmaskedZ(T_data, z_data);
    
    % Cut zK and K0 to the size according to umaskedZ from z_data
    [zKmasked, K0, zKMflag] = initK(zK, umaskedZ);

    % cut the data according to the range of K
    [T_data, z_data, indCutZ] = cutData(T_data, z_data, [zKmasked(1),zKmasked(end)]);
   
    % Take weights into account
    if ((yearIndex == 1) || (yearIndex == 3))
        T_S = dataS{yearIndex}.T_S(indCutZ, indt);
    elseif (yearIndex == 2)
        T_S = data.T_sd(indCutZ, indt);
        T_S = T_S.^2;
    else 
        T_S = ones(size(T_data));
    end
    
    w = 1./ (T_S.^(0.5));
    
    %% Optimize for the averaged value
    Nz = length(z_data);
    K_opt = inverseK(data, dataIndex, zKmasked, K0, Nz, rho, w, perturbT, timePeriod);

    %% Solve Heat equation

    % Set initial and boundary conditions
    [Tbc, T0, z, t, dz, Nt, dt] = setIBCs(z_data, t_data, Nz, T_data, interpOption);

    % Set Parameters for solving
    heatParam = setHeatParam(dt, Nt, dz, Nz, rho, C, T0, Tbc, zKmasked);

    [T_sol] = solveHeat(t, z, K_opt, heatParam);

    %% Create output vector
    K_opt_out = [zK(:), zK(:)];
    K_opt_out(zKMflag, 2) = K_opt;
	K_opt_out(~zKMflag, 2) = nan;
    
    %% Visualize measurement
    % Scale t_data to days
%     t_data = scaleTimeUnit(t_data);
    T_sol_ondata = project2D(T_sol, t, z, t_data, z_data);

    % mesh for measurment
    [X_data, Y_data] = meshgrid(t_data, z_data);

    figure
    subplot(2, 1, 1)
    surf(X_data, Y_data, T_data)
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    caxis([-20, -2]);
    xlabel('t')
    ylabel('z')
    title(['Measurements at 201', num2str(yearIndex+1)]);
    grid off
    

    % Scale t to days
%     t = scaleTimeUnit(t);
    % mesh for numerical ex
%     [X, Y] = meshgrid(t, z);
    subplot(2,1,2)
%     surf(X,Y,T_sol)
    surf(X_data, Y_data, abs(T_data-T_sol_ondata))
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    caxis([0,0.1]);
    xlabel('t')
    ylabel('z')
    title(['Errors at 201', num2str(yearIndex+1)]);
    grid off

end
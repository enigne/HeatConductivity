function K_opt = testheat(yearIndex)
    %% Initialize
    % Settings
    interpOption = 'linear';

    %% Load data
    load('LF_4_aver.mat');
    load('densityData.mat');
    load('summary.mat');

    % Assign data
    data = LF{yearIndex}.T;
    rho = rhoData{yearIndex};

    % Load Measurements
    [t_data, z_data, T_data] = loadData(data);
    
    % Physical parameters
    C = 152.5 + 7.122 * (273.15 - 10);

    %% Initialize
    Nk = 5;
    zK = linspace(1, 8, Nk)';
    K0 = 1e5*ones(Nk, 1);

    % cut the data according to the range of K
    [T_data, ~, indCutZ] = cutData(T_data, z_data, [zK(1),zK(end)]);
   
    % Take weights into account
    if ((yearIndex == 1) || (yearIndex == 3))
        T_S = dataS{yearIndex}.T_S(indCutZ, :);
    else 
        T_S = ones(size(T_data));
    end
    w = 1./ (T_S.^(0.5));

    z_data = z_data(indCutZ);
    
    %% Optimize for the averaged value
    Nz = length(z_data);
    K_opt = inverseK(data, 0, zK, K0, Nz, rho, w);

    %% Solve Heat equation

    % Set initial and boundary conditions
    [Tbc, T0, z, t, dz, Nt, dt] = setIBCs(z_data, t_data, Nz, T_data, interpOption);

    % Set Parameters for solving
    heatParam = setHeatParam(dt, Nt, dz, Nz, rho, C, T0, Tbc.Up, Tbc.Down, zK);

    [T_sol] = solveHeat(t, z, K_opt, heatParam);

    %% Visualize measurement
    % mesh for measurment
    [X_data, Y_data] = meshgrid(t_data, z_data);

    figure
    subplot(2, 1, 1)
    surf(X_data, Y_data, T_S)
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
%     caxis([-20, -2]);
    xlabel('t')
    ylabel('z')
    title(['Measurements at 201', num2str(yearIndex+1)]);
    grid off
    

    [X, Y] = meshgrid(t, z);
    subplot(2,1,2)
    surf(X,Y,T_sol)
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    caxis([-20,-2]);
    xlabel('t')
    ylabel('z')
    title(['Optimal solution at 201', num2str(yearIndex+1)]);
    grid off

end
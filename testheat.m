function K_opt = testheat(yearIndex)
    %% Initialize
    % Settings
    interpOption = 'linear';

    %% Load data
    load('LF_4_aver.mat');
    load('densityData.mat');

    %%
    data = LF{yearIndex}.T;
    rho = rhoData{yearIndex};

    dIndex = 0;

    % Measurements
    % Time vector (row)
    t_data = data.t';
    % Position vector (column)
    z_data = data.z_a;
    % z_data = data.z_corr{dIndex};
    % Temperature
    T_data = data.T_a;
    % T_data = data.T_corr{dIndex};

    % % Add noise
    % noise = 0.0*(1-2*rand(size(T_data)));
    % T_data = T_data + noise;

    % Physical parameters
    C = 152.5 + 7.122 * (273.15 - 10);


    %% Initialize
    Nk = 5;
    zK = linspace(0, 8, Nk)';
    K0 = 0.4e5*ones(Nk, 1);

    % cut the data according to the range of K
    [T_data, z_data] = cutData(T_data, z_data, [zK(1),zK(end)]);

    %% Optimize for the averaged value
    Nz = length(z_data);
    K_opt = inverseK(data, 0, zK, K0, Nz, rho);

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
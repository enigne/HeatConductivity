function [K_opt] = inverseK(data, dataIndex, zK)
%% Initialize

% Settings
interpOption = 'linear';

%
% % Load the data set
% load('LF_4_aver.mat');
% data = LF{1,1}.T;
% dataIndex = 4;

% Measurements
% Time vector (row)
t_data = data.t';
% Position vector (column)
z_data = data.z_corr{dataIndex};
% Temperature
T_data = data.T_corr{dataIndex};

% Physical parameters
rho = 900;
C = 152.5 + 7.122 * (273.15 - 10);

% mesh for measurment
% [X_data, Y_data] = meshgrid(t_data, z_data);

%% Visualize measurement
% figure
% subplot(2, 1, 1)
% surf(X_data, Y_data, T_data)
% view(2)
% shading interp;
% colorbar
% colormap(jet)
% axis tight
% caxis([-20, -2]);
% grid off

%% Solve Heat equation
% Spactial Discretization
Nx = 201;
z = linspace(min(z_data), max(z_data),Nx)';

% Time discretization same size as measurment
Nt = length(t_data);
t = linspace(min(t_data), max(t_data),Nt);
dt = abs(t(2) - t(1));

% Mesh for solving PDE
% [X, Y] = meshgrid(t, z);

% Initial Heat conductivity
Nk = length(zK);
% zK = linspace(0, 12, Nk)';
Kp = 0.4e6*ones(Nk, 1);

% Interpolate initial condition
T0 = interp1(z_data, T_data(:,1), z, interpOption);

% Interpolate Boundary conditions
TbcUp = interp1(t_data, T_data(1,:), t, interpOption);
TbcDown = interp1(t_data, T_data(end,:), t, interpOption);

heatParam = setHeatParam(dt, rho, C, T0, TbcUp, TbcDown, zK);

% Project data to the computational domain
f_data = projectY2D(T_data, t_data, z_data, t, z);

%% Optimisation

f = @(K) objF(K, @solveHeat, heatParam, z, t, f_data);

% options = optimoptions('lsqnonlin','Display','iter','typicalX', Kp,'TolFun',1e-10,'TolX',1e-10);
options = optimoptions('lsqnonlin','Display','iter','typicalX', Kp);
[K_opt,resnorm,residual,exitflag,output] = lsqnonlin(f, Kp, [], [],options);

%% Plot optimal solution
% [T_sol] = solveHeat(t, z, K_opt, heatParam);
%
% subplot(2,1,2)
% surf(X,Y,T_sol)
% view(2)
% shading interp;
% colorbar
% colormap(jet)
% axis tight
% caxis([-20,-2]);
% grid off
end
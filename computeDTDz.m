clear
close all;
%%

load('invK.mat');
K0 = K_opt_ave(:,end);


%% Initialize

% Settings
interpOption = 'linear';

% Load the data set
load('LF_4_aver.mat');
data = LF{1,1}.T;

% Measurements
% Time vector (row)
t_data = data.t';
% Position vector (column)
z_data = data.z_a;
% Temperature
T_data = data.T_a;
    
% Physical parameters
rho = 900;
C = 152.5 + 7.122 * (273.15 - 10);

%% Solve Heat equation
% Spatial Discretization
Nfine = 200;
Nx = Nfine * (length(z_data)-1) + 1;
z = linspace(min(z_data), max(z_data),Nx)';

% Time discretization same size as measurment
Nt = length(t_data);
t = linspace(min(t_data), max(t_data),Nt);
dt = abs(t(2) - t(1));

% Initial Heat conductivity
Nk = 5;
zK = linspace(0, 12, Nk)';

% Interpolate initial condition
T0 = interp1(z_data, T_data(:,1), z, interpOption);

% Interpolate Boundary conditions
TbcUp = interp1(t_data, T_data(1,:), t, interpOption);
TbcDown = interp1(t_data, T_data(end,:), t, interpOption);

heatParam = setHeatParam(dt, rho, C, T0, TbcUp, TbcDown, zK);

% Project data to the computational domain
f_data = projectY2D(T_data, t_data, z_data, t, z);

%% Plot optimal solution
T0_sol = solveHeat(t, z, K0, heatParam);

%% Compute dT/dz
xInd = [1:Nfine:Nx];
xIndP1 = xInd + 1;
xIndP1(end) = Nx;
xIndM1 = xInd - 1;
xIndM1(1) = 1;

dz = abs(z(2)-z(1));
Tcenter = T0_sol(xInd, :);
Tup = T0_sol(xIndP1, :);
Tdown = T0_sol(xIndM1, :);

% one-side finite differences
dTdz = (Tup - Tcenter) ./ dz;

% %% Visualize measurement
% % mesh for measurment
% [X_data, Y_data] = meshgrid(t_data, z_data);
% 
% figure
% subplot(2, 1, 1)
% surf(X_data, Y_data, Tcenter)
% view(2)
% shading interp;
% colorbar
% colormap(jet)
% axis tight
% caxis([-20, -2]);
% grid off
% 
% 
% subplot(2, 1, 2)
% surf(X_data, Y_data, dTdz)
% view(2)
% shading interp;
% colorbar
% colormap(jet)
% axis tight
% caxis([-10, 10]);
% grid off
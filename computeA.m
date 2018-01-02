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
grid off

%% Solve Heat equation
% Spactial Discretization
Nx = 601;
z = linspace(min(z_data), max(z_data),Nx)';

% Time discretization same size as measurment
Nt = length(t_data);
t = linspace(min(t_data), max(t_data),Nt);
dt = abs(t(2) - t(1));

% Mesh for solving PDE
[X, Y] = meshgrid(t, z);

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

subplot(2,1,2)
surf(X,Y,T0_sol)
view(2)
shading interp;
colorbar
colormap(jet)
axis tight
caxis([-20,-2]);
grid off


%%
dK = 0.1;
dKvec = dK * speye(Nk);
A = zeros(Nx*Nt, Nk);

for i = 1 : Nk
    T_sol_temp = solveHeat(t, z, K0+dKvec(:,i), heatParam);
    dT_temp = (T0_sol - T_sol_temp) /dK;
    A(:, i) = dT_temp(:);
end


%%
B = inv(A'*A);
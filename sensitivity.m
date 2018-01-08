clear
close all;
%% Initialize
% Settings
interpOption = 'linear';

% Load optimal K
% load('invK_1e8_Nz241.mat');
load('invK_1e8_Nz601.mat');
K0 = K_opt_ave(:,end);
% Initial Heat conductivity
Nk = length(K0);
zK = linspace(0, 12, Nk)';

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

%% Compute dTdz
dZfine = 2;
[dTdz, matDTDz] = computeDTDz(z_data, t_data, T_data, dZfine, zK, K0, rho, C, interpOption);

%% Compute A  
dZfine = 5;
dK = 0.001;
AK = computeA(z_data, t_data, T_data, dZfine, zK, K0, dK, rho, C, interpOption);
B = inv(AK'*AK);

%% Compute rho
dZfine = 5;
dRho = 0.001;
NRho = 5;
rho = 900 * ones(NRho,1);
zRho = linspace(0, 12, NRho)';

ARho = computeRhoSensitivity(z_data, t_data, T_data, dZfine, zK, K0, dRho, rho, zRho, C, interpOption);
%%
D = B * AK' * ARho

% %%
% figure
% C = B*A'*matDTDz;
% Nz = length(z_data);
% R = tril(ones(Nz));
% D = C*R;
% subplot(2,1,1)
% plot(C');
% subplot(2,1,2)
% plot(D');

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
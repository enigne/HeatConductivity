clear
close all;
%% Initialize
% Settings
interpOption = 'linear';

% Load optimal K
% load('invK_1e8_Nz241.mat');
% load('invK_1e8_Nz601.mat');
load('invK_0_8.mat');
K0 = K_opt_ave(:,end);
% Initial Heat conductivity
Nk = length(K0);
zK = linspace(0, 8, Nk)';

% Load the data set
load('LF_4_aver.mat');
data = LF{1,1}.T;

% Load average and mean
load('ML_data.mat');

% Measurements
% Time vector (row)
t_data = data.t';
% Position vector (column)
z_data = data.z_a;
% Temperature
T_data = data.T_a;

% cut the data according to the range of K
[T_data, z_data] = cutData(T_data, z_data, [zK(1),zK(end)]);
[T_E, ~] = cutData(T_E, z_data, [zK(1),zK(end)]);
[T_S, ~] = cutData(T_S, z_data, [zK(1),zK(end)]);

% mask for T < -2
mask = find( T_data < 20 );

% Physical parameters
rho = 900;
C = 152.5 + 7.122 * (273.15 - 10);

%% Compute dTdz
dZfine = 5;
[dTdz, matDTDz] = computeDTDz(z_data, t_data, T_data, dZfine, zK, K0, rho, C, mask, interpOption);

%% Compute A  
dZfine = 10;
dK = 0.001;
AK = computeKSensitivity(z_data, t_data, T_data, dZfine, zK, K0, dK, rho, C, mask, interpOption);
B = inv(AK'*AK);
%%
figure
zALeg = {};
for i = 1: Nk
    zALeg{i} = ['K', num2str(i)];
end
Az = B * AK' * matDTDz;

Nz = length(z_data);
R = tril(ones(Nz));
D = Az*R;
subplot(2,1,1)
plot(z_data, Az');
xlim([0, max(z_data)])
% legend(zALeg);
xlabel('z');
ylabel('Az');
subplot(2,1,2)
plot(z_data, D');
xlim([0, max(z_data)])
xlabel('z');
ylabel('Az*R');
legend(zALeg);


%% Compute rho
dZfine = 5;
dRho = 0.001;
NRho = 5;
rho = 900 * ones(NRho,1);
zRho = linspace(0, 8, NRho)';

ARho = computeRhoSensitivity(z_data, t_data, T_data, dZfine, zK, K0, dRho, rho, zRho, C, mask, interpOption);
%%
D = - B * AK' * ARho

%% compute weight
w = 1./ (T_S.^(0.5));
W = w(mask);
weightedAK = AK' * spvardiag(W);
weightB = inv(weightedAK * AK);

% % %% Visualize measurement
% mesh for measurment
[X_data, Y_data] = meshgrid(t_data, z_data);

figure3
indicator = 1.* (T_data < -2);
subplot(2, 1, 1)
surf(X_data, Y_data, w.*indicator)
view(2)
shading interp;
colorbar
colormap(jet)
axis tight
% caxis([-20, -2]);
grid off


subplot(2, 1, 2)

%%
figure
surf(X_data, Y_data,  reshape(AK(:,4),size(X_data)))
view(2)
shading interp;
colorbar
colormap(jet)
axis tight
% caxis([-1, 1]);
grid off
clear
close all;
%% Initialize
% Settings
interpOption = 'linear';
yearIndex = 1;

%% Load data
load('LF_4_aver.mat');
load('densityData.mat');
load('invK_realRho.mat');

data = LF{yearIndex}.T;
rho = rhoData{yearIndex};
K0 = z_05_8_ondata.K_opt{yearIndex};

% Initial Heat conductivity
Nk = length(K0);

% Measurements
[t_data, z_data, T_data] = loadData(data);

% z-coordinates of K
zK = linspace(0.5, 8, Nk)';

% Cut the data according to the range of K
[T_data, z_data] = cutData(T_data, z_data, [zK(1),zK(end)]);

% Cut the time series with Nan in T_data
[T_data, t_data, noNanInd] = cutNan(T_data, t_data);

% Load average and mean
load('summary.mat');
[T_E, ~] = cutData(dataS{yearIndex}.T_E, z_data, [zK(1),zK(end)]);
[T_S, ~] = cutData(dataS{yearIndex}.T_S, z_data, [zK(1),zK(end)]);

T_E = T_E(: , noNanInd);
T_S = T_S(: , noNanInd);

% mask for T >= -2 put 0
mask = find( T_data >= -2);

% Physical parameters
C = 152.5 + 7.122 * (273.15 - 10);

%% Compute dTdz
dZfine = 5;
[dTdz, matDTDz, ~, t] = computeDTDz(z_data, t_data, T_data, dZfine, zK, K0, rho, C, mask, interpOption);

%% Compute A  
dZfine = 5;
dK = 0.001;
A = computeKSensitivity(z_data, t_data, T_data, dZfine, zK, K0, dK, rho, C, mask, interpOption);
B = A'*A;
AK = B \ A';

%% compute weight
w = 1./ (T_S.^(0.5));
W = w(:);
W(mask) = 0;
weightedA = A' * spvardiag(W);
weightedB = weightedA * A;
weightedAK = weightedB \ weightedA;
%% Plot AK and A
figure
[X_data, Y_data] = meshgrid(t_data, z_data);
for i = 1 : 6
    if i > 1
        p_data = reshape(AK(i-1,:), size(X_data));
%         p_data = reshape(A(:, i-1), size(X_data));
        subTitle = ['K',num2str(i-1)];
    else
        p_data =  T_data;
        subTitle = ['Measurements in 201',num2str(yearIndex+1)];
    end
    
    p_data(mask) = 0;
    
    subplot(3, 2, i)
    surf(X_data, Y_data, p_data);
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    title(subTitle);
    if i == 1
        caxis([-20, -2]);
    else
        caxis([-1, 1]);

    end
    grid off
end
%% Plot Weighted AK and A
figure
for i = 1 : 6
    if i > 1
%         p_data = reshape(weightedA(i-1,:), size(X_data));
        p_data = reshape(weightedAK(i-1,:), size(X_data));
        subTitle = ['K',num2str(i-1)];
    else
        p_data = reshape(W, size(X_data));
        subTitle = ['Weights in 201',num2str(yearIndex+1)];
    end
    
    p_data(mask) = 0;
    
    subplot(3, 2, i)
    surf(X_data, Y_data, p_data);
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    title(subTitle);
    if i == 1
        caxis([0, 15]);
    else
        caxis([-6, 6]);
    end
    grid off
end

%% Plot Az and Az*R
figure
zALeg = {};
for i = 1: Nk
    zALeg{i} = ['K', num2str(i)];
end
Az = B \ A' * matDTDz;

Nz = length(z_data);
R = tril(ones(Nz));
D = Az*R;
subplot(2,1,1)
plot(z_data, Az');
xlim([min(z_data), max(z_data)])
% legend(zALeg);
xlabel('z');
ylabel('Az');
subplot(2,1,2)
plot(z_data, D');
xlim([min(z_data), max(z_data)])
xlabel('z');
ylabel('Az*R');
legend(zALeg);


%% Compute rho
dZfine = 5;
dRho = 0.001;
ARho = computeRhoSensitivity(z_data, t_data, T_data, dZfine, zK, K0, dRho, rho, C, mask, interpOption);


%%
% D = - B * AK' * ARho

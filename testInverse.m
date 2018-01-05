clear
close all

%% Load data
load('LF_4_aver.mat');
data = LF{1,1}.T;
dataIndex = [1:9];
nInd = length(dataIndex);

%% Initialize
Nk = 5;
zK = linspace(0, 12, Nk)';
K0 = 0.4e6*ones(Nk, 1);

%% Optimize for the averaged value with different resolution
N_opt_ave = [31, 61, 121, 241, 361, 481, 601];
K_opt_ave = zeros(Nk, length(N_opt_ave));

for i = 1 : length(N_opt_ave)
    K0 = inverseK(data, 0, zK, K0, N_opt_ave(i));
    K_opt_ave(:,i) = K0;
end

%% Optimize for each hole
Nz = 241;
K_opt_ind = zeros(Nk, nInd);
for i = 1:nInd
    K0 = inverseK(data, dataIndex(i), zK, K0, Nz);
    K_opt_ind(:,i) = K0;
end

%% Optimize for the averaged value with noises
noise = [0.1, 0.2, 0.5, 1, 2];
Nz = 241;
K_opt_noiseT = zeros(Nk, length(noise));

for i = 1 : length(noise)
    K0 = inverseK(data, 0, zK, K0, Nz, noise(i));
    K_opt_noiseT(:,i) = K0;
end




%%
% figure
% for i = 1:nInd
%    plot(zK, K_opt_ind(:,i));
%    hold on
% end

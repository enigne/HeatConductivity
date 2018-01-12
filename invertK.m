clear
close all

%% Load data
load('LF_4_aver.mat');
data = LF{1,1}.T;
dataIndex = [1:9];
nInd = length(dataIndex);

%% Initialize
Nk = 5;
zK = linspace(0, 8, Nk)';
K0 = 0.4e6*ones(Nk, 1);

%% Optimize for the averaged value with different resolution
N_opt_ave = [241];
K_opt_ave = zeros(Nk, length(N_opt_ave));

for i = 1 : length(N_opt_ave)
    K0 = inverseK(data, 4, zK, K0, N_opt_ave(i));
    K_opt_ave(:,i) = K0;
end

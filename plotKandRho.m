% Script for plot optimal K and rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% options
year = 2;
dataIndex = 1;
saveFig = 0;
saveK_Rho = 1;

%% Load data
% Predefined parameters
NK = 8;
NRho = NK;
gamma = [ 1e-3, 1e-1, 1e1, 1e3];
plotMarkers = {'-d', '-o', '-^', '-*', '-+'};

% load measurements
try
    load('LF_4_aver.mat');
catch
    error('Please check the original data file LF_4_aver.mat');
end

% Load density data
try
    load('densityData.mat');
catch
    error('densityData.mat not found. Try to run preprocessRho.m first.');
end
    
%% Go through all the test cases
fig = figure('pos',[0 0 600 400]);
for i = 1:length(gamma)
    % load Opt K according to Nk
    dataFileName = ['invK', num2str(NK), 'rho', num2str(NRho), '_gamma', num2str(gamma(i)), '_maskedBC_longP.mat'];
    
    try
        load(dataFileName);
    catch
        error([dataFileName, ' not found. Try to run testInverseHeat.m with Nk=', num2str(NK), ...
            ' NRho=', num2str(NRho), ' gamma=', num2str(gamma(i)), '.']);
    end
    
    % fetch data at specific year
    K = K_opt{year, dataIndex};
    rho = rho_opt{year, dataIndex};
    t_data = t_data_opt{year, dataIndex};
    nanFlagK = isnan(K(:,2));
    nanFlagRho = isnan(rho(:,2));
    K = K(~nanFlagK, :);
    rho = rho(~nanFlagRho, :);
    
    % plot K
    subplot(2, 1, 1)
    plot(K(:,1), K(:,2), plotMarkers{i}, 'linewidth', 1.5)
    hold on
    
    % plot Rho
    subplot(2, 1, 2)
    plot(rho(:,1), rho(:,2), plotMarkers{i}, 'linewidth', 1.5)
    hold on
    
    legendList{i} = ['$\gamma$=', num2str(gamma(i))];
end


% add labels and legends
subplot(2, 1, 1)
title(['Optimal K in 201', num2str(year+1)],'Interpreter','latex');
xlim([1, 8])
ylim([0, 2.5])
xlabel('z')
ylabel('K')
legend(legendList, 'Location', 'northwest','Interpreter','latex')

% plot measured rho in subplot(2)
subplot(2, 1, 2)
plot(rhoData{year}.z, rhoData{year}.rho, 'linewidth', 1)
title(['Optimal $\rho$ in 201', num2str(year+1)],'Interpreter','latex');
xlim([1, 8])
ylim([300, 900])
xlabel('z')
ylabel('$\rho$','Interpreter','latex')
legendList{i+1} = 'data';
legend(legendList, 'Location', 'northwest','Interpreter','latex')

%% save figures
if saveFig
    print(fig, ['Figures/Optimal_K', num2str(NK), '_rho', num2str(NRho),'_201', num2str(year+1)], '-depsc');
end

%%
figK_Rho = figure('pos',[0 0 400 300]);
dataIndex = 1;
yearInds = [1:4];
gamma = 10;
for i = 1:length(yearInds)
    % load Opt K according to Nk
    dataFileName = ['invK', num2str(NK), 'rho', num2str(NRho), '_gamma', num2str(gamma), '_maskedBC_longP.mat'];
    
    try
        load(dataFileName);
    catch
        error([dataFileName, ' not found. Try to run testInverseHeat.m with Nk=', num2str(NK), ...
            ' NRho=', num2str(NRho), ' gamma=', num2str(gammai), '.']);
    end
    
    % fetch data at specific year
    K = K_opt{yearInds(i), dataIndex};
    rho = rho_opt{yearInds(i), dataIndex};
    t_data = t_data_opt{yearInds(i), dataIndex};
    nanFlagK = isnan(K(:,2));
    nanFlagRho = isnan(rho(:,2));
    K = K(~nanFlagK, :);
    rho = rho(~nanFlagRho, :);
    
    % plot K
    plot(K(:,2), rho(:,2), plotMarkers{i}(2), 'linewidth', 1.5)
    hold on
	legendKRhoList{i} = ['201', num2str(yearInds(i)+1)];
end
xlim([0, 2])
ylim([400, 650])

xlabel('K')
ylabel('$\rho$','Interpreter','latex')
legend(legendKRhoList, 'Location', 'northwest','Interpreter','latex')
if saveK_Rho
    print(figK_Rho, ['Figures/Optimal_K', num2str(NK), '_To_rho', num2str(NRho)], '-depsc');
end
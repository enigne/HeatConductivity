% Script for sensitivity analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
close all
clear

%% Initialize
% Predefined parameters
NK = 8;
NRho = NK;
gamma = 10;
saveOptK = 1;

% load Opt K according to Nk
% optKFileName = ['invK', num2str(Nk), '_maskedBC.mat'];
% optKFileName = ['invK', num2str(NK), '_maskedBC_longP.mat'];
optKFileName = ['invK', num2str(NK), 'rho', num2str(NRho),'_gamma', num2str(gamma), '_maskedBC_longP.mat'];
try
    load(optKFileName);
catch
    error([optKFileName, ' not found. Try to run testInverseHeat.m first.']);
end
yearIndex = [1:4];


%% Test sensitivity
for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
        % get specific data
        zK = K_opt{i,j}(:,1);
        K = K_opt{i,j}(:,2);
        zRho = rho_opt{i,j}(:,1);
        rho = rho_opt{i,j}(:,2);
        t = t_data_opt{i,j};
        
        % remove nan in K and zK
        nanFlag = isnan(K);
        K = K(~nanFlag);
        zK = zK(~nanFlag);
        
        % solve for sensitivity
        [weightedAK{i,j}, weightedB{i,j}, weightedSE{i,j}, weightedAz{i,j}, weightedD{i,j},...
        weightedARho{i,j}, dTdz{i,j},T_data{i,j}, mask{i,j}, t_cell{i,j}, z_cell{i,j}] = ...
            solveSensitivity(yearIndex(i), K, zK, timePeriods{yearIndex(i)}{j}, rho, zRho);
    end
end

%% Plot OptK and error bar
figOptK = figure('pos',[0 0 900 600]);

for i = 1: length(yearIndex)
    n = 1;
    legendList= {};
    subplot(2, 2, i)
    for j = 1: length(timePeriods{yearIndex(i)})
        zK = K_opt{i,j}(:,1);
        K = K_opt{i,j}(:,2);
        nanFlag = isnan(K);
        K = K(~nanFlag);
        zK = zK(~nanFlag);
        errorbar(zK, K, weightedSE{i,j}, 'linewidth', 1.5);
        hold on;
        t_conv = scaleTimeUnit(t_data_opt{i,j},'','');
        daytemp = datestr(t_conv,'yyyy-mm-dd');
        legendList{n} = [daytemp(1,:),' to ', daytemp(2,:)];
        n = n+1;
    end
    title(['Optimal K in 201', num2str(yearIndex(i)+1)]);
    xlim([1, 8])
    ylim([0, 2.5])
    xlabel('z')
    ylabel('K')
    legend(legendList)
end
if saveOptK
    print(figOptK, ['Figures/Optimal_K', num2str(NK), '_error' ], '-depsc');
end

%% Save data
dataFileName = ['sensitivity_K', num2str(NK), '_maskedBC_longP.mat'];
save(dataFileName, 'weightedAK', 'weightedAz', 'weightedB', 'weightedD', 'weightedSE', ...
    'weightedARho', 'K_opt','t_data_opt', 't_cell', 'z_cell', 'timePeriods', 'yearIndex');

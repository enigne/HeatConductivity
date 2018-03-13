% Script for testing inverse solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear
close all

%% Solve for K
yearIndex = [1:4];
dataIndex = 0;
K_opt = {};
NK = 5;
zK = linspace(1, 8, NK);

timePeriods{1} = {[0, 59]./181}; % 2012
timePeriods{2} = {[0, 40]./80}; % 2013
timePeriods{3} = {[0, 78]./360, [160, 360]./360}; % 2014
timePeriods{4} = {[1, 84]./360}; % 2015

% timePeriods{1} = {[0, 20]./181, [20, 40]./181, [40, 59]./181, [60, 73]./181, [76, 90]./181}; % 2012
% timePeriods{2} = {[0, 20]./80, [20, 40]./80, [47, 59]./80}; % 2013
% timePeriods{3} = {[0, 20]./360, [20, 40]./360, [40, 60]./360, [60, 78]./360}; % 2014
% timePeriods{4} = {[1, 20]./360, [20, 40]./360, [40, 60]./360, [60, 84]./360}; % 2015


for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
        for l = 1:length(dataIndex)
            [K_opt_temp, t_data] = solveInverseHeat(yearIndex(i), dataIndex(l), zK, timePeriods{yearIndex(i)}{j});
            K_opt{i, j, l} = K_opt_temp;
            t_data_opt{i, j, l} = [t_data(1), t_data(end)];
        end
    end
end

%% Plot Optimal K
figure('pos',[0 0 900 600])
for i = 1: length(yearIndex)
    % plot
    legendList = {};
    n = 1;
    subplot(2, 2, i)
    for j = 1: length(timePeriods{yearIndex(i)})
        for l = 1:length(dataIndex)         
            plot(K_opt{i,j, l}(:, 1) , K_opt{i,j, l}(:,2), 'linewidth', 1.5);
            hold on;
            t_conv = scaleTimeUnit(t_data_opt{i,j,l},'','');
            daytemp = datestr(t_conv,'yyyy-mm-dd');
            legendList{n} = [daytemp(1,:),' to ', daytemp(2,:)];
            n = n+1;
        end
    end
    xlim([1, 8])
    ylim([0, 2.5])
    xlabel('z')
    ylabel('K')
    legend(legendList)
end

%% Save data
dataFileName = ['invK', num2str(NK), '_maskedBC_longP.mat'];
save(dataFileName, 'K_opt', 't_data_opt', 'timePeriods', 'yearIndex');

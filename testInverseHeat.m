% Script for testing inverse solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear
close all

%% Solve for K
yearIndex = [1:4];
dataIndex = 0;
K_opt = {};
NK = 10;
zK = linspace(1, 8, NK);

timePeriods{1} = {[0, 3.8/12], [4/12, 4.8/12], [5/12, 6.2/12]}; % 2012
timePeriods{2} = {[0, 1/2], [2.35/4, 2.95/4]}; % 2013
timePeriods{3} = {[0, 5/36], [5/36, 9/36], [18/36,1]}; % 2014
timePeriods{4} = {[0, 2/9], [3/9, 1]}; % 2015


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
            plot(zK , K_opt{i,j, l}, 'linewidth', 1.5);
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
dataFileName = ['invK', num2str(NK), '_maskedBC.mat'];
save(dataFileName, 'K_opt', 't_data_opt', 'timePeriods', 'yearIndex');

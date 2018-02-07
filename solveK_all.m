clear
close all

%%
yearIndex = [1];
dataIndex = 0;
K_opt = {};
NK = 5;
zK = linspace(1, 8, NK);
timePeriods{1} = {[0, 1/3], [1/3, 5/12], [5/12, 1]}; % 2012
timePeriods{2} = {[0, 1/2], [1/2, 2.95/4]}; % 2013
timePeriods{3} = {[0, 5/36], [5/36, 13/36], [13/36,1]}; % 2014
timePeriods{4} = {[0, 2/9]}; % 2015


% timePeriods{1} = {[0, 1/3], [0, 1], [1/3, 1]}; % 2012
% timePeriods{2} = {[0, 2/3], [0, 1], [2/3, 1]}; % 2013
% timePeriods{3} = {[0, 5/36], [0, 5/18], [5/36, 1/3], [0, 1/3], [1/3,1], [0,1]}; % 2014
% timePeriods{4} = {[0, 1/6], [0,1],[0, 5/18], [0,5/12], [5/12,1]}; % 2015


for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
        for l = 1:length(dataIndex)
            [K_opt_temp, t_data] = testheat(yearIndex(i), dataIndex(l), zK, timePeriods{yearIndex(i)}{j});
            K_opt{i, j, l} = K_opt_temp;
            t_data_opt{i, j, l} = [t_data(1), t_data(end)];
        end
    end
end

%%
for i = 1: length(yearIndex)
    % plot
    legendList = {};
%     figure

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
    ylim([0, 3.5])
    xlabel('z')
    ylabel('K')
    legend(legendList)
end

%%
    
%%
% load('invK_realRho.mat')
%
% optK = z_1_8_weighted.K_mat(:, yearIndex);
% averK = mean(K_mat,2);
% stdK = std(K_mat, 0, 2);
% figure
% errorbar(zK, averK, stdK, 'linewidth', 1.5);
% hold on
% plot(zK, optK, 'linewidth', 1.5)
% xlim([1, 8])
% ylim([0, 3.5])
% xlabel('z')
% ylabel('K')
% legend({'averaged OptK', 'OptK using averaged data'})

%
% legend({'2012', '2013', '2014', '2015'})
%
% z_0_.K_mat = K_mat;
% z_0_.K_opt = K_opt;

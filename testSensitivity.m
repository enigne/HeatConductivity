% Script for sensitivity analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
close all
clear

%% Initialize
% Predefined parameters
Nk = 15;

% load Opt K according to Nk
optKFileName = ['invK', num2str(Nk), '_maskedBC.mat'];
load(optKFileName);
yearIndex = [1];

% define zK
if ~exist('zK')
    zK = linspace(1, 8, Nk)';
end
%% 
for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
      K = K_opt{i,j};
      t = t_data_opt{i,j};
      
       [weightedB{i,j}, weightedSE{i,j}, weightedAz{i,j}, weightedD{i,j},dTdz{i,j}, T_data{i,j}, mask{i,j}] = solveSensitivity(yearIndex(i), K, zK, timePeriods{yearIndex(i)}{j});
       
    end
end

%%
figure
for i = 1: length(yearIndex)
    n = 1;
    legendList= {};
    subplot(2, 2, i)
    for j = 1: length(timePeriods{yearIndex(i)})
            errorbar(zK, K_opt{i,j}, weightedSE{i,j}, 'linewidth', 1.5);
            hold on;
            t_conv = scaleTimeUnit(t_data_opt{i,j},'','');
            daytemp = datestr(t_conv,'yyyy-mm-dd');
            legendList{n} = [daytemp(1,:),' to ', daytemp(2,:)];
            n = n+1;
    end
    xlim([1, 8])
    ylim([0, 2])
    xlabel('z')
    ylabel('K')
    legend(legendList)
end




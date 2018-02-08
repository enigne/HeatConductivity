close all
clear

%%
load('invK_seasonal.mat');
yearIndex = [1:4];

zK = linspace(1, 8, 5)';
%%

for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
      K = K_opt{i,j};
      t = t_data_opt{i,j};
      
       [weightedB{i,j}, weightedSE{i,j}, weightedAz{i,j}, weightedD{i,j},dTdz{i,j}, T_data{i,j}, mask{i,j}] = testSensitivity(yearIndex(i), K, zK, timePeriods{yearIndex(i)}{j});
       
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




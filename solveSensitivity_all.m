close all
clear

%%
load('invK_seasonal.mat');
yearIndex = 1;

zK = linspace(1, 8, 5)';
%%

for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
      K = K_opt{i,j};
      t = t_data_opt{i,j};
      
       [weightedB{i,j}, weightedSE{i,j}, weightedAz{i,j}, weightedD{i,j},dTdz{i,j}, T_data{i,j}, mask{i,j}] = testSensitivity(yearIndex(i), K, zK, timePeriods{yearIndex(i)}{j});
       
    end
end
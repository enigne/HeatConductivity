clear 
close all

%%
yearIndex = 3;
dataIndex = 0;
K_opt = {};
NK = 5;
zK = linspace(1, 8, NK);
timePeriod = [0, 0.4];

for i = 1: length(yearIndex)
    for j = 1:length(dataIndex)
        K_opt{i,j} = testheat(yearIndex(i), dataIndex(j), zK, timePeriod);
    end
end

%%
K_mat = cell2mat(K_opt) ;

figure
for i = 1: length(dataIndex)
    plot(zK , K_mat(:, i), 'linewidth', 1.5)
    hold on
end
xlim([1, 8])
ylim([0, 3.5])
xlabel('z')
ylabel('K')

%%

load('invK_realRho.mat')

optK = z_1_8_weighted.K_mat(:, yearIndex);
averK = mean(K_mat,2);
stdK = std(K_mat, 0, 2);
figure
errorbar(zK, averK, stdK, 'linewidth', 1.5);
hold on 
plot(zK, optK, 'linewidth', 1.5)
xlim([1, 8])
ylim([0, 3.5])
xlabel('z')
ylabel('K')
legend({'averaged OptK', 'OptK using averaged data'})


% legend({'2012', '2013', '2014', '2015'})
% 
% z_0_.K_mat = K_mat;
% z_0_.K_opt = K_opt;

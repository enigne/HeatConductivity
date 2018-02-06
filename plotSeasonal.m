clear
% close all

%% Load data
load('invK_seasonal.mat')
load('LF_4_aver.mat')

%% Plot K
% figure
subplot(2, 1, 2)
n = 1;

for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
%         K(:, 2*n-1:2*n) = kron(K_opt{i,j}, [1,1]);
%         t(1, 2*n-1:2*n) = t_data_opt{i,j};
%         t(2*n) = t(2*n) - 1;
                K(:, n) = K_opt{i,j};
                t(n) = mean(t_data_opt{i,j});
        n = n + 1;
    end
end
t = scaleTimeUnit(t);
legendList={};
for i = 1 : 5
    legendList{i} = ['K',num2str(i)];
end

plot(t, K)
ylim([0, 2.5])
xlim([0,max(t)]);
xlabel('t (days)')
ylabel('K')
legend(legendList);

%%

% for i = 1: length(yearIndex)
%     data = LF{yearIndex(i)}.T;
%     [t_data, z_data, T_data, ~] = loadData(data);
%     % mesh for measurment
%     t_data = t_data/24/3600 - LF{1}.T.t(1)+1;
%     
%     [X_data, Y_data] = meshgrid(t_data, z_data);
%     
%     subplot(2, 1, 2)
%     surf(X_data, Y_data, T_data)
%     hold on
%     view(2)
%     shading interp;
%     %     colorbar
%     colormap(jet)
%     axis tight
%     caxis([-20, -2]);
%     xlabel('t (days)')
%     ylabel('z')
%     grid off
% end
% 
% xlim([0,max(t)]);
% ylim([1,8]);

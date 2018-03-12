% Script for plot optimal K, error bar and T_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Load data
% Predefined parameters
NK = 5;

% load Opt K according to Nk
% optKFileName = ['invK', num2str(Nk), '_maskedBC.mat'];
optKFileName = ['sensitivity_K', num2str(NK), '_maskedBC_longP.mat'];
load(optKFileName);
% load measurements
load('LF_4_aver.mat')

%% Generate K field from K_opt
z_offset = [0, -0.5, -1.3, -2.2];
n = 1;
z_cord = [-3:0.1:10]';

for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
        Z(:, 2*n-1:2*n) = kron((K_opt{i,j}(:,1) + z_offset(i)), [1,1]);
        K(:, 2*n-1:2*n) = kron(K_opt{i,j}(:,2), [1,1]);
        t(1, 2*n-1:2*n) = t_data_opt{i,j};
        n = n + 1;
    end
end
t = scaleTimeUnit(t);
t_cord = [t(1):1:t(end)];


%% project K values onto t-z domain
[tMesh, zMesh] = meshgrid(t_cord, z_cord);
KMesh = nan*zeros(size(tMesh));

Nt = length(t);

for i = 1:1: Nt-1
    ind = find(((t_cord <t(i+1)) & (t_cord >= t(i))));
    [tTemp, zTemp] = meshgrid(t(i:i+1), mean(Z(:, i:i+1), 2));
    tempK = interp2(tTemp, zTemp, K(:,i:i+1),tMesh(:,ind), zMesh(:,ind), 'linear');
    KMesh(:, ind) = tempK;
end 


%% Plot K
figure('pos',[0 0 900 600])
subplot(2, 1, 1)
surf(tMesh, zMesh, KMesh)
hold on
view(2)
shading interp;
colorbar
colormap(jet)
axis tight
% caxis([-20, -2]);
xlabel('t (days)')
ylabel('z')
grid off


% 
% legendList={};
% for i = 1 : 5
%     legendList{i} = ['K',num2str(i)];
% end
% 
% plot(t, K)
% ylim([0, 2.5])
% xlim([0,max(t)]);
% xlabel('t (days)')
% ylabel('K')
% legend(legendList);

%%
subplot(2, 1, 2);
for i = 1: length(yearIndex)
    data = LF{yearIndex(i)}.T;
    [t_data, z_data, T_data, ~] = loadData(data);
    % mesh for measurment
    t_data = t_data/24/3600 - LF{1}.T.t(1)+1;
    z_data = z_offset(i) + z_data;
    
    [X_data, Y_data] = meshgrid(t_data, z_data);
    
    subplot(2, 1, 2)
    surf(X_data, Y_data, T_data)
    hold on
    view(2)
    shading interp;
        colorbar
    colormap(jet)
    axis tight
    caxis([-20, -2]);
    xlabel('t (days)')
    ylabel('z')
    grid off
end

xlim([0,max(t)]);
ylim([-2,8]);

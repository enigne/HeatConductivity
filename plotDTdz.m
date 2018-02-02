close all
clear
%%
interpOption = 'linear';
yearIndex = 1 :4;

%% Load data
load('LF_4_aver.mat');
figure

for i = 1: length(yearIndex)
    data = LF{yearIndex(i)}.T;
    [t_data, z_data, T_data] = loadData(data);

    % Cut data
    zK = linspace(1, 8, 5)';

    % Cut the data according to the range of K
    [T_data, z_data] = cutData(T_data, z_data, [zK(1),zK(end)]);
    [T_data, t_data, noNanInd] = cutNan(T_data, t_data);

    %
    dZfine = 5;
    [dTdzData, ~, ~, ~] = computeDTDzFromData(z_data, t_data, T_data, dZfine, [], interpOption);
    t_data = (t_data-t_data(1)+1)/24/3600;
    [X_data, Y_data] = meshgrid(t_data, z_data);

    %% Plot
    subplot(2, 2, i)
    surf(X_data, Y_data, dTdzData);
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    xlabel('t (days)');
    ylabel('z(m)');
    title(['temperature gradient 201', num2str(yearIndex(i)+1)])
    caxis([-4, 4])
end

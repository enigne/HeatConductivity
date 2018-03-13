% Script for preprocessing the standard deviations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%%
yearIndex = [1,3];
%% Load data
load('LF_4_aver.mat');

%%
for y = 1: length(yearIndex)
    data = LF{yearIndex(y)}.T;
    dataIndex = [1: 9];
    nInd = length(dataIndex);
    
    % Time vector (row)
    t = data.t';
    % Position vector (column)
    z_a = data.z_a;
    % Temperature
    T_a = data.T_a;
    
    T_S = zeros(size(T_a));
    T_E = zeros(size(T_a));
    
    for i = 1 : nInd
        % Position vector (column)
        z_data = data.z_i{dataIndex(i)};
        % Temperature
        T_data = data.T_i{dataIndex(i)};
        
        [T_data, tcut, noNanInd] = cutNan(T_data, t);
        % Project data to the computational domain
        %     f_data = project2D(T_data, tcut, z_data, t, z_a);
        
        f_data = project2DY(T_data, z_data, z_a);
        
        res = (f_data - T_a(:, noNanInd)) .^ 2;
        
        %     T_E = T_E + f_data;
        T_S(:, noNanInd) = T_S(:, noNanInd) + res;
    end
    
    % T_E = T_E / (nInd);
    T_S = T_S / (nInd - 1);
    
    
    %%
    [X_data, Y_data] = meshgrid(t, z_a);
    figure
    subplot(2,1,1)
    surf(X_data, Y_data, data.T_sd );
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    % caxis([-1e-14, 1e-14])
    caxis([0, 2.5])
    
    subplot(2,1,2)
    surf(X_data, Y_data, T_S.^0.5);
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    caxis([0, 2.5])
    %%
    dataS{yearIndex(y)}.T_S = T_S;
    dataS{yearIndex(y)}.T_E = T_a;
end

%%
fileName = 'summary.mat';
save(fileName, 'dataS');
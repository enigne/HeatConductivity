clear
close all

%%
yearIndex = 3;
%% Load data
load('LF_4_aver.mat');
data = LF{yearIndex}.T;
dataIndex = [1: 9];
nInd = length(dataIndex);

%%
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
    
    % Project data to the computational domain
    f_data = project2D(T_data, t, z_data, t, z_a);
    
    res = (f_data - T_a) .^ 2;
    
    T_E = T_E + f_data;
    T_S = T_S + res;
end

T_E = T_E / nInd;
T_S = T_S / (nInd - 1);


%% 
[X_data, Y_data] = meshgrid(t, z_a);
surf(X_data, Y_data, T_E-T_a);
view(2)
shading interp;
colorbar
colormap(jet)
axis tight

%%
dataS{yearIndex}.T_S = T_S;
dataS{yearIndex}.T_E = T_E;


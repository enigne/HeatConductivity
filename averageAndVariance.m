clear
close all

%% Load data
load('LF_4_aver.mat');
data = LF{1,1}.T;
dataIndex = [1:9];
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
    z_data = data.z_corr{dataIndex(i)};
    % Temperature
    T_data = data.T_corr{dataIndex(i)};
    
    % Project data to the computational domain
    f_data = projectY2D(T_data, t, z_data, t, z_a);
    
    res = (f_data - T_a) .^ 2;
    
    T_E = T_E + f_data;
    T_S = T_S + res;
end

T_E = T_E / nInd;
T_S = T_S / (nInd - 1);
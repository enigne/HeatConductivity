clear 
% close all

%%
yearIndex = 3;
dataIndex = 0;
noise.noise = 0.1;
noise.zInd = 31;
noise.tInd = 3300;

K_ref = testheat(yearIndex, 0);

K_opt = testheat(yearIndex, 0, noise);

%%
% K_mat = cell2mat(K_opt) ./ 24./3600;
% z = linspace(1, 8, 5);
% 

((K_opt-K_ref)/24/3600) /noise.noise*1e4
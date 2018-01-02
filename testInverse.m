clear
close all

%%
load('LF_4_aver.mat');
data = LF{1,1}.T;
dataIndex = [1:9];
nInd = length(dataIndex);

%%
Nk = 5;
zK = linspace(0, 12, Nk)';
Kmat = zeros(Nk, nInd);

for i = 1:nInd
    Kmat(:,i) = inverseK(data, dataIndex(i), zK);
end

%%
figure

for i = 1:nInd
   plot(zK, Kmat(:,i));
   hold on
end

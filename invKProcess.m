clear 
close all
%%
load('invK_1e8_Nz241.mat');

[Nz, Nk] = size(K_opt_ind);
z = linspace(0,12,Nz);
legendList = cell(Nk,1);
figure
for i = 1 : Nk
    plot(z, K_opt_ind(:,i));
    hold on 
    legendList{i} = ['Index = ', num2str((i))];
end
legend(legendList);

%%
Ntest = size(K_opt_ave,2);
legendList = cell(Ntest,1);
figure
for i = 1 : Ntest
    plot(z, K_opt_ave(:,i));
    hold on 
    legendList{i} = ['Nx = ', num2str(N_opt_ave(i))];
end

legend(legendList);

%%

Nnoise = length(noise);
legendList = cell(Nnoise,1);
figure
for i = 1 : Nnoise
    plot(z, K_opt_noiseT(:,i));
    hold on 
    legendList{i} = ['Noise = ', num2str(noise(i))];
end

legend(legendList);
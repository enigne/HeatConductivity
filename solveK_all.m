clear 
close all

%%
yearIndex = 1:4;
K_opt = {};

for i = 1: length(yearIndex)
    K_opt{i} = testheat(yearIndex(i));
end

%%
K_mat = cell2mat(K_opt) ./ 24./3600;
z = linspace(0, 8, 5);

figure
for i = 1: length(yearIndex)
    plot(z , K_mat(:, i))
    hold on
end
xlim([0, 8])
ylim([0, 3.5])
xlabel('z')
ylabel('K')
legend({'2012', '2013', '2014', '2015'})
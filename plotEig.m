% Script for plot optimal K and rho and eigenvalues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2019-03-29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% Load data
% Predefined parameters
gamma = [1e-5, 1e1, 1e3 1e7];
plotMarkers = {'-d', '-o', '-^', '-*', '-+', '-',':'};

% load measurements
try
    load('gammaTest.mat');
catch
    error('gammaTest.mat not found!');
end

%% Go through all the test cases
fig = figure('pos',[0 0 800 600]);
for i = 1:length(gamma)
    % plot K
    subplot(2, 2, 1)
    plot(z, Klines(:,i), plotMarkers{i}, 'linewidth', 1.5)
    hold on

    % plot Rho
    subplot(2, 2, 3)
    plot(z, rholines(:,i), plotMarkers{i}, 'linewidth', 1.5)
    hold on

    % eigenvectors
    subplot(2, 2, 2)
	plot(eigVals(:,i), plotMarkers{i}(2), 'linewidth', 1.5)
    hold on
    
    subplot(2, 2, 4)
	semilogy(abs(eigVals(:,i)), plotMarkers{i}(2), 'linewidth', 1.5)
    hold on
    % legends
    legendList{i} = ['$\gamma$=', num2str(gamma(i),2)];
end


% add labels and legends
subplot(2, 2, 1)
title(['Optimal K in 2012'],'Interpreter','latex');
xlim([1, 8])
ylim([0, 2.5])
xlabel('z')
ylabel('K')
legend(legendList, 'Location', 'best','Interpreter','latex')

subplot(2, 2, 2)
xlim([1,16])
xlabel('v')
title('eigenvectors corresponds to the smallest eigenvalue','Interpreter','latex')   
legend(legendList, 'Location', 'best','Interpreter','latex')

subplot(2, 2, 4)
xlim([1,16])
xlabel('v')
title('eigenvectors corresponds to the smallest eigenvalue','Interpreter','latex')   
legend(legendList, 'Location', 'best','Interpreter','latex')


% plot measured rho in subplot(2)
subplot(2, 2, 3)
plot(rhoDataZ, rhoDataLine, 'linewidth', 1)
title(['Optimal $\rho$ in 2012'],'Interpreter','latex');
xlim([1, 8])
ylim([300, 900])
xlabel('z')
ylabel('$\rho$','Interpreter','latex')
legendList{i+1} = 'data';
legend(legendList, 'Location', 'northwest','Interpreter','latex')

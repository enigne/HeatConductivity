% Script for plot optimal K and rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2019-03-21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% options
year = 1;
dataIndex = 1;
saveFig = 1;
saveK_Rho = 1;

%% Load data
% Predefined parameters
NK = 8;
NRho = NK;
gamma = [1e-5, 1e1, 1e3 1e7];
plotMarkers = {'-d', '-o', '-^', '-*', '-+', '-',':'};

% load measurements
try
    load('LF_4_aver.mat');
catch
    error('Please check the original data file LF_4_aver.mat');
end

% Load density data
try
    load('densityData.mat');
catch
    error('densityData.mat not found. Try to run preprocessRho.m first.');
end

% load Opt K according to Nk
dataFileName = ['sensitivity_K', num2str(NK), '_halfK.mat'];
try
    load(dataFileName);
catch
    error([dataFileName, ' not found. Try to run testSensitivity.m with Nk=', num2str(NK), ...
        ' NRho=', num2str(NRho)'.']);
end

% Load matrix for F
Q11 = weightedB{year, 1};
Q12 = A12{year, 1};
Q21 = Q12';
Q22 = A22{year, 1};

%%
z = [1:NK]';
Klines = zeros(NK,length(gamma));
rholines = zeros(NK,length(gamma));
eigVals= zeros(NK+NRho,length(gamma));

%% Go through all the test cases
fig = figure('pos',[0 0 800 600]);
for i = 1:length(gamma)
    % load Opt K according to Nk
    dataFileName = ['invK', num2str(NK), 'rho', num2str(NRho), '_gamma', num2str(gamma(i)), '_maskedBC_longP.mat'];
    
    try
        load(dataFileName);
    catch
        error([dataFileName, ' not found. Try to run testInverseHeat.m with Nk=', num2str(NK), ...
            ' NRho=', num2str(NRho), ' gamma=', num2str(gamma(i)), '.']);
    end
    
    % fetch data at specific year
    K = K_opt{year, dataIndex};
    rho = rho_opt{year, dataIndex};
    t_data = t_data_opt{year, dataIndex};
    nanFlagK = isnan(K(:,2));
    nanFlagRho = isnan(rho(:,2));
    K = K(~nanFlagK, :);
    rho = rho(~nanFlagRho, :);
    
    % plot K
    subplot(2, 2, 1)
    plot(K(:,1), K(:,2), plotMarkers{i}, 'linewidth', 1.5)
    hold on
    Klines(:,i) = K(:,2);
    % plot Rho
    subplot(2, 2, 3)
    plot(rho(:,1), rho(:,2), plotMarkers{i}, 'linewidth', 1.5)
    hold on
    rholines(:,i) = rho(:,2);
    % plot the eigenvalues of [A', B']*W*[A; B]

    Q = [Q11, Q12;
        Q21, Q22+gamma(i)*eye(NRho)];
    [V,D]= eig(Q);
    [d,ind] = sort(diag(D),'descend');
    V = V(:,ind);
    
    subplot(2, 2, 2)
	plot((V(1:end,end)), plotMarkers{i}, 'linewidth', 1.5)
    hold on
    eigVals(:, i) = V(1:end,end);
    
    subplot(2, 2, 4)
	semilogy(abs(V(1:end,end)), plotMarkers{i}, 'linewidth', 1.5)
    hold on
    % legends
    legendList{i} = ['$\gamma$=', num2str(gamma(i),2)];
end


% add labels and legends
subplot(2, 2, 1)
title(['Optimal K in 201', num2str(year+1)],'Interpreter','latex');
xlim([1, 8])
ylim([0, 2.5])
xlabel('z')
ylabel('K')
legend(legendList, 'Location', 'best','Interpreter','latex')

subplot(2, 2, 2)
xlim([1,16])
xlabel('z')
title('eigenvectors corresponds to the smallest eigenvalue','Interpreter','latex')   
legend(legendList, 'Location', 'best','Interpreter','latex')

subplot(2, 2, 4)
xlim([1,16])
xlabel('z')
title('eigenvectors corresponds to the smallest eigenvalue','Interpreter','latex')   
legend(legendList, 'Location', 'best','Interpreter','latex')


% plot measured rho in subplot(2)
subplot(2, 2, 3)
plot(rhoData{year}.z, rhoData{year}.rho, 'linewidth', 1)
rhoDataLine = rhoData{year}.rho;
rhoDataZ = rhoData{year}.z;
title(['Optimal $\rho$ in 201', num2str(year+1)],'Interpreter','latex');
xlim([1, 8])
ylim([300, 900])
xlabel('z')
ylabel('$\rho$','Interpreter','latex')
legendList{i+1} = 'data';
legend(legendList, 'Location', 'northwest','Interpreter','latex')

%% save figures
if saveFig
    print(fig, ['Figures/Optimal_K', num2str(NK), '_rho', num2str(NRho),'_201', num2str(year+1)], '-depsc');
end

%%
% figK_Rho = figure('pos',[0 0 400 300]);
% dataIndex = 1;
% yearInds = [1:4];
% gamma = 10;
% for i = 1:length(yearInds)
%     % load Opt K according to Nk
%     dataFileName = ['invK', num2str(NK), 'rho', num2str(NRho), '_gamma', num2str(gamma), '_maskedBC_longP.mat'];
%     
%     try
%         load(dataFileName);
%     catch
%         error([dataFileName, ' not found. Try to run testInverseHeat.m with Nk=', num2str(NK), ...
%             ' NRho=', num2str(NRho), ' gamma=', num2str(gammai), '.']);
%     end
%     
%     % fetch data at specific year
%     K = K_opt{yearInds(i), dataIndex};
%     rho = rho_opt{yearInds(i), dataIndex};
%     t_data = t_data_opt{yearInds(i), dataIndex};
%     nanFlagK = isnan(K(:,2));
%     nanFlagRho = isnan(rho(:,2));
%     K = K(~nanFlagK, :);
%     rho = rho(~nanFlagRho, :);
%     
%     % plot K
%     plot(K(:,2), rho(:,2), plotMarkers{i}(2), 'linewidth', 1.5)
%     hold on
% 	legendKRhoList{i} = ['201', num2str(yearInds(i)+1)];
% end
% xlim([0, 2])
% ylim([400, 650])
% 
% xlabel('K')
% ylabel('$\rho$','Interpreter','latex')
% legend(legendKRhoList, 'Location', 'northwest','Interpreter','latex')
% if saveK_Rho
%     print(figK_Rho, ['Figures/Optimal_K', num2str(NK), '_To_rho', num2str(NRho)], '-depsc');
% end
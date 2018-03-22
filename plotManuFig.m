% Script for generating plots in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
%% Load data
% Predefined parameters
NK = 8;
NRho = NK;
gamma = 10;
z_offset = [3.34, 2.22, 0.9, 0];

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
dataFileName = ['sensitivity_K', num2str(NK), '_maskedBC_longP.mat'];
try
    load(dataFileName);
catch
    error([dataFileName, ' not found. Try to run testSensitivity.m with Nk=', num2str(NK), ...
        ' NRho=', num2str(NRho)'.']);
end

%% Plot K with errorbar
figEbar = figure('pos',[0 0 600 300]);
n = 1;
for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
        % title
        if j == 1
            term = 'spring';
        else
            term = 'fall';
        end
        
        subplot(1, 5, n)
        zK = K_opt{i,j}(:,1) + z_offset(i);
        K = K_opt{i,j}(:,2);
        zRho = rho_opt{i,j}(:,1) + z_offset(i);
        rho = rho_opt{i,j}(:,2);
        nanFlag = isnan(K);
        K = K(~nanFlag);
        zK = zK(~nanFlag);
        
        % plot error bar
        line(K, zK,'linewidth',2);
        hold on
        %         errorbar(K, zK, weightedSE{i,j},'horizontal','LineStyle','none','linewidth',2);
        errorbar(K, zK, weightedSE{i,j},'LineStyle','none','linewidth',2);
        
        xlim([0,2.5]);
        ylim([1,12]);
        view([0, -90]);
        xlabel('$K$','Interpreter','latex');
        %         ylabel('$z$','Interpreter','latex');
        
        ax1 = gca;
        ax1.XColor = 'b';
        ax2 = axes('Position',get(ax1,'Position'),...
            'XAxisLocation','top',...
            'YAxisLocation','right',...
            'Color','none',...
            'XColor','r','YColor','k');
        %Eliminating second y-Axis
        set(gca,'ytick',[])
        ax2.XColor = 'r';
        
        line(rho, zRho,'Parent',ax2,'Color','r','linewidth',1.5);
        xlabel('$\rho$','Interpreter','latex');
        xlim([350,750]);
        ylim([1,12]);
        view([0, -90]);
        
        title(['201',num2str(i+1), ' ', term]);
        
        n = n + 1;
    end
end

%%
% matlab2tikz('optK_rho8.tex','height', '\fheight', 'width', '\fwidth' );







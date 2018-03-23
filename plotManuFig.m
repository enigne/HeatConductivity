% Script for generating plots in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
%% options
plotKandRho = 0;
plotT = 0;
%% Load data
% Predefined parameters
NK = 8;
NRho = NK;
gamma = 10;
z_offset = [3.34, 2.22, 0.9, 0];
timePeriods{1} = {[0, 59]./181}; % 2012
timePeriods{2} = {[0, 40]./80}; % 2013
timePeriods{3} = {[0, 78]./360, [160, 360]./360}; % 2014
timePeriods{4} = {[1, 84]./360}; % 2015

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

%% Plot K with errorbar and rho
if plotKandRho
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
    % Export to tex file
    % matlab2tikz('optK_rho8.tex','height', '\fheight', 'width', '\fwidth' );
end
%% Plot T Data and rho data
if plotT
    fig = figure('pos',[0 0 700 250]);
    t_bc = zeros(5,2);
    n = 1;
    gapConst = 40;
    t_offset = 0;
    t_tick = [];
    for i = 1: length(yearIndex)
        for j = 1: length(timePeriods{yearIndex(i)})
            
            data = LF{yearIndex(i)}.T;
            [t_data, z_data, T_data, ~] = loadData(data, 0, timePeriods{i}{j});
            % mesh for measurment
            t_data = t_data/24/3600 - LF{1}.T.t(1)+1;
            t_bc(n, :) = [t_data(1), t_data(end)];
            
            if n > 1
                t_offset = t_offset + t_bc(n, 1)-t_bc(n-1, 2);
                if n == 5
                    t_offset = t_offset - 4.45;
                else
                    t_offset = t_offset - gapConst;
                end
            end
            t_tick(n) = t_offset;
            
            z_data = z_offset(i) + z_data;
            
            [X_data, Y_data] = meshgrid(t_data-t_offset, z_data);
            
            surf(X_data, Y_data, T_data)
            hold on
            n = n + 1;
        end
    end
    view([0,-90])
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    caxis([-20, -2]);
    grid off
    xlabel('t (days)')
    ylabel('$z$', 'Interpreter','latex')
    xlim([0,max(t_data-t_offset)]);
    ylim([1,11.5]);
    title('Temperature measurements')
    
    % Compute ticks
    ticks = floor(reshape((t_bc - repmat(t_tick',1, 2))', 1,  2*(n-1)));
    tickLabel = (floor(reshape(t_bc', 1,  2*(n-1))));
    
    % remove too closed ticks
    hideTMask = find(diff(ticks) < gapConst/2);
    ticks(hideTMask+1) = [];
    tickLabel(hideTMask+1) = [];
    
    xticks(ticks);
    xticklabels(tickLabel);
    
    
    % Save png
    print(fig, ['Figures/tempMeasurements'], '-dpng', '-r600');
end
%%







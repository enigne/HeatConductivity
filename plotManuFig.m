% Script for generating plots in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-04-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
%% options
plotKandRho = 0;
plotT = 0;
plotRho = 0;
plotAk = 0;
plotAz = 0;
plotTRho = 0;
plotFitRhoK = 0;
printEAz = 1;
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
dataFileName = ['sensitivity_K', num2str(NK), '_halfK.mat'];
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
            %                     errorbar(K, zK, weightedSE{i,j},'horizontal','LineStyle','none','linewidth',2);
            %             errorbar(K, zK, weightedSE{i,j},'LineStyle','none','linewidth',2);
            errorbar(K, zK, weightedSE_t_indep{i,j},'LineStyle','none','linewidth',2);
            
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
    matlab2tikz('optK_rho8.tex','height', '\fheight', 'width', '\fwidth' );
end
%% Plot T Data and rho data
if plotT
    fig = figure('pos',[0 0 600 200]);
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
            t_data = t_data - t_offset;
            
            [X_data, Y_data] = meshgrid(t_data, z_data);
            
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
    xlim([0,max(t_data)]);
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
    set(gca,'fontsize',12);
    
    ax = gca;
    ax.XAxis.FontSize = 8;
    
    % Save png
    print(fig, ['Figures/tempMeasurements'], '-dpng', '-r600');
end
%% Plot rho data
if plotRho
    figRho = figure('pos',[0 0 600 300]);
    for i = 1: length(yearIndex)
        subplot(1, length(yearIndex), i)
        % optimal rho
        zRho = rho_opt{i,1}(:,1) + z_offset(i);
        rho = rho_opt{i,1}(:,2);
        % measurements
        zmin = min(zRho);
        zmax = max(zRho);
        
        z_mea = rhoData{i}.z + z_offset(i);
        rho_mea = rhoData{i}.rho;
        optRange = ((z_mea <= zmax) & (z_mea >= zmin));
        % plot
        plot(rho_mea(optRange), z_mea(optRange), 'LineWidth' ,1)
        hold on
        plot(rho, zRho, 'LineWidth' ,1.5)
        xlim([350,750]);
        ylim([1,12]);
        %         view([0, -90]);
        xlabel('$\rho$','Interpreter','latex');
        title(['201',num2str(i+1)])
    end
    
    matlab2tikz('measuredRho.tex','height', '\fheight', 'width', '\fwidth' );
end

%% Plot AK
if plotAk && plotT % need to run plotT to remove the gap
    figAk = figure('pos',[0 0 736 253]);
    for k = 1: NK
        %         subplot(ceil(NK/2), 2, k)
        %             figAk = figure('pos',[0 0 736 253]);
        
        n = 1;
        for i = 1: length(yearIndex)
            for j = 1: length(timePeriods{yearIndex(i)})
                if k <= size(weightedAK{i,j}, 1)
                    t = (t_cell{i,j}-t_cell{1,1}(1))/24/3600+1-t_tick(n);
                    z = z_cell{i,j}+  z_offset(i);
                    Ak = reshape(weightedAK{i,j}(k,:), length(z), length(t));
                    [tmesh, zmesh] = meshgrid(t,z);
                    surf(tmesh, zmesh, Ak)
                    hold on
                    n = n + 1;
                end
            end
        end
        view([0,-90])
        shading interp;
        colormap(jet)
        axis tight
        
        xticks(ticks);
        if k <= NK -2
            xticklabels([]);
        else
            xticklabels(tickLabel);
            xlabel('t (days)')
        end
        
        if mod(k,2)
            ylabel('z')
        else
            yticklabels([]);
            %                     colorbar
        end
        
        grid off
        %         title(['$A_{K',num2str(k),'}$'], 'Interpreter','latex')
        caxis([-5e-4, 5e-4]);
        ylim([1,11.5])
        xlim([ticks(1), ticks(end)])
        set(gca,'FontSize',12);
        %
        %         ax = gca;
        %         ax.XAxis.FontSize = 8;
        hAxes = gca;
        hAxes.XRuler.Axle.LineStyle = 'none';
        axis off
        hold off
        pause(0.5)
        print(figAk, ['Figures/sensitivityT/Sensitivity_AK_Nk', num2str(NK),'K',num2str(k) ], '-dpng', '-r600');
        
    end
    %     colorbar('southoutside')
    %     print(figAk, ['Figures/Sensitivity_AK_Nk', num2str(NK),'_leg' ], '-dpng', '-r600');
    % 	matlab2tikz('Sensitivity_AK_Nk.tex','height', '\fheight', 'width', '\fwidth' );
end

%% Plot Az and Az*R
if plotAz
    figAz = figure('pos',[0 0 120 300]);
    n = 1;
    for i = 1: length(yearIndex)
        for j = 1: length(timePeriods{yearIndex(i)})
            figAz = figure('pos',[0 0 120 300]);
            tempAz = weightedAz{i,j};
            tempD = weightedD{i,j};
            [nK, nz] = size(tempAz);
            
            % get date
            t_conv = scaleTimeUnit(t_data_opt{i,j},'','');
            daytemp = datestr(t_conv,'yyyy-mm-dd');
            
            
            zK = linspace(1, 8, nK);
            z_data = 1:0.1:8;
            
            zK = zK(1:nK);
            z_data = z_data(1:nz)+ 0* z_offset(i);
            
            % plot
            [X, Y] = meshgrid(zK, z_data);
            surf(X, Y, tempD')
            view([0,-90])
            shading interp;
            %             colorbar
            colormap(jet)
            grid off
            xlim([1,8])
            ylim([1,8])
            xlabel('$z_K$', 'Interpreter','latex')
            ylabel('deepth')
            caxis([-5, 5])
            %             title(['$A_zR$ from ', daytemp(1,:),' to ', daytemp(2,:) ], 'Interpreter','latex');
            set(gca,'FontSize',12);
            hAxes = gca;
            hAxes.XRuler.Axle.LineStyle = 'none';
            axis off
            print(figAz, ['Figures/sensitivityZ/Sensitivity_z', num2str(NK),'_',num2str(n)], '-dpng', '-r600');
            
            n = n +1;
        end
    end
    %     colorbar('southoutside')
    %         matlab2tikz('sensitivityZ.tex','height', '\fheight', 'width', '\fwidth' );
    %     print(figAz, ['Figures/sensitivityZ/Sensitivity_z_leg'], '-dpng', '-r600');
    
end

%%
if plotTRho
    %     for i = 1: length(yearIndex)
    %         for j = 1: length(timePeriods{yearIndex(i)})
    i = 1; j = 1;
    figRho = figure('pos',[0 0 600 600]);
    tempARho = weightedARho{i,j}*1e4;
    zK = 1:size(tempARho, 1);
    zRho = 1:size(tempARho, 2);
    [X,Y] = meshgrid(zRho, zK);
    surf(X, Y, tempARho)
    view([0,-90])
    shading interp;
    %     colorbar('southoutside')
    xlabel('$z_K$', 'Interpreter','latex')
    ylabel('deepth')
    colormap(jet)
    grid off
    caxis([-10, 10])
    hAxes = gca;
    hAxes.XRuler.Axle.LineStyle = 'none';
    axis off
    %         end
    %     end
    print(figRho, ['Figures/Sensitivity_Rho'], '-dpng', '-r600');
    %     print(figRho, ['Figures/Sensitivity_Rho_leg'], '-dpng', '-r600');
    
end
%%
if plotFitRhoK
    figK_Rho = figure('pos',[0 0 600 400]);
    plotMarkers = {'-d', '-x', '-^', '-*', '-+'};
    K_vec = [];
    rho_vec = [];
    
    n = 1;
    for i = 1: length(yearIndex)
        for j = 1: length(timePeriods{yearIndex(i)})
            % fetch data at specific year
            K = K_opt{yearIndex(i), j}(1:end-1,:);
            rho = rho_opt{yearIndex(i), j}(1:end-1,:) ./1e3;
            err = weightedSE{i,j}(1:end-1,:);
            
            nanFlagK = isnan(K(:,2));
            nanFlagRho = isnan(rho(:,2));
            K = K(~nanFlagK, :);
            rho = rho(~nanFlagRho, :);
            
            K_vec = [K_vec;K(:,2)];
            rho_vec = [rho_vec;rho(:,2)];
            
            % plot K
            plot(rho(:,2), K(:,2),  plotMarkers{n}(2), 'linewidth', 2)
            hold on
            %             errorbar(K(:,2), rho(:,2), err,'horizontal','LineStyle','none')
            if j == 2
                legendKRhoList{n} = ['201', num2str(yearIndex(i)+1),'-2'];
            else
                legendKRhoList{n} = ['201', num2str(yearIndex(i)+1)];
            end
            n = n +1;
        end
    end
    ylim([0, 2])
    xlim([400, 650]./1e3)
    
    p = polyfit( rho_vec,K_vec, 2);
    x = [min(rho_vec):0.01:max(rho_vec)];
    y = polyval(p, x);
    plot(x, y,'--')
    
    legendKRhoList{n} = ['$k$=', num2str(p(1),'%.2f'),'$\rho^2$+',num2str(p(2),'%.2f'),'$\rho$', num2str(p(3),'%.2f')];
    ylabel('k')
    xlabel('$\rho$','Interpreter','latex')
    legend(legendKRhoList, 'Location', 'best','Interpreter','latex')
    matlab2tikz('K_rho.tex','height', '\fheight', 'width', '\fwidth' );
end

%% Print out the expected value and variance from Az
if printEAz
    epsilon = 1/800;
    for i = 1: length(yearIndex)
        for j = 1: length(timePeriods{yearIndex(i)})
            Az = weightedD{i,j};
            eMn = ones(size(Az,2),1);
            Edk = -epsilon * Az*eMn;
            Vdk = sqrt(diag(Az*Az'));
            disp(num2str(Edk', '%3.3f & '));
            disp(num2str(Vdk', '%3.3f & '));
            disp(' ');
        end
    end
    
    
end
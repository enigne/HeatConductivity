% Script for plot optimal K, error bar and T_data, Az and Az*R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%% options

% deepth sensitivity
plotAz = 0;
saveAz = 0;

% T sensitivity
plotAk = 0;
saveAk = 0;

% K errors
plotOptK = 0;
saveK = 0;

% rho sensitivity
printARho = 1;

%% Load data
% Predefined parameters
NK = 5;

% load Opt K according to Nk
% optKFileName = ['invK', num2str(Nk), '_maskedBC.mat'];
optKFileName = ['sensitivity_K', num2str(NK), '_maskedBC_longP.mat'];
try
    load(optKFileName);
catch
    error([optKFileName, ' not found. Try to run testSensitivity.m first.']);
end
% load measurements
load('LF_4_aver.mat')

%% Generate K field from K_opt
z_offset = [3.34, 2.22, 0.9, 0];
n = 1;
z_cord = [0:0.1:12]';

for i = 1: length(yearIndex)
    for j = 1: length(timePeriods{yearIndex(i)})
        Z(:, 2*n-1:2*n) = kron((K_opt{i,j}(:,1) + z_offset(i)), [1,1]);
        K(:, 2*n-1:2*n) = kron(K_opt{i,j}(:,2), [1,1]);
        t(1, 2*n-1:2*n) = t_data_opt{i,j};
        n = n + 1;
    end
end
t = scaleTimeUnit(t);
t_cord = [t(1):1:t(end)];


%% project K values onto t-z domain
[tMesh, zMesh] = meshgrid(t_cord, z_cord);
% t = repmat(t, NK, 1);
% KMesh = griddata(t, Z, K, tMesh, zMesh, 'linear');

KMesh = nan*zeros(size(tMesh));

Nt = length(t);

for i = 1:2: Nt-1
    ind = find(((t_cord <t(i+1)) & (t_cord >= t(i))));
    [tTemp, zTemp] = meshgrid(t(i:i+1), mean(Z(:, i:i+1), 2));
    tempK = interp2(tTemp, zTemp, K(:,i:i+1),tMesh(:,ind), zMesh(:,ind), 'linear');
    KMesh(:, ind) = tempK;
end

%% Plot Opt K and data
if plotOptK
    fig = figure('pos',[0 0 900 600]);
    subplot(2, 1, 1)
    surf(tMesh, zMesh, KMesh)
    hold on
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    % caxis([-20, -2]);
    xlabel('t (days)')
    ylabel('z')
    grid off
    ylim([1,11.5]);
    title('Optimal K')
    caxis([0, 2]);
    
    subplot(2, 1, 2);
    for i = 1: length(yearIndex)
        data = LF{yearIndex(i)}.T;
        [t_data, z_data, T_data, ~] = loadData(data);
        % mesh for measurment
        t_data = t_data/24/3600 - LF{1}.T.t(1)+1;
        z_data = z_offset(i) + z_data;
        
        [X_data, Y_data] = meshgrid(t_data, z_data);
        
        subplot(2, 1, 2)
        surf(X_data, Y_data, T_data)
        hold on
    end
    view(2)
    shading interp;
    colorbar
    colormap(jet)
    axis tight
    caxis([-20, -2]);
    grid off
    xlabel('t (days)')
    ylabel('z')
    xlim([0,max(t_cord)]);
    ylim([1,11.5]);
    title('Temperature measurements') 
    
	if saveK
        print(fig, ['Figures/Optimal_K', num2str(NK), '_z_t'], '-dpng', '-r600');
	end
end
%% Plot Az and Az*R
if plotAz
    for i = 1: length(yearIndex)
        for j = 1: length(timePeriods{yearIndex(i)})
            tempAz = weightedAz{i,j};
            tempD = weightedD{i,j};
            [nK, nz] = size(tempAz);
            
            % get date
            t_conv = scaleTimeUnit(t_data_opt{i,j},'','');
            daytemp = datestr(t_conv,'yyyy-mm-dd');
            
            
            zK = linspace(1, 8, nK);
            z_data = 1:0.1:8;
            
            zK = zK(1:nK);
            z_data = z_data(1:nz);
            
            % plot
            [X, Y] = meshgrid(zK, z_data);
            fig = figure;
            subplot(2, 1, 1)
            surf(X, Y, tempAz')
            view(2)
            shading interp;
            colorbar
            colormap(jet)
            axis tight
            grid off
            xlabel('$z_K$', 'Interpreter','latex')
            ylabel('deepth')
            caxis([-1, 1])
            title(['$A_z$ from ', daytemp(1,:),' to ', daytemp(2,:) ], 'Interpreter','latex');
            
            subplot(2, 1, 2)
            surf(X, Y, tempD')
            view(2)
            shading interp;
            colorbar
            colormap(jet)
            axis tight
            grid off
            xlabel('$z_K$', 'Interpreter','latex')
            ylabel('deepth')
            caxis([-5, 5])
            title(['$A_zR$ from ', daytemp(1,:),' to ', daytemp(2,:) ], 'Interpreter','latex');
            if saveAz
                print(fig, ['Figures/Sensitivity_z', num2str(NK),'_',daytemp(1,:)], '-dpng', '-r600');
            end
        end
    end
end

%% Plot AK
if plotAk
    for k = 1: NK
        fig = figure('pos',[0 0 900 300]);
        for i = 1: length(yearIndex)
            for j = 1: length(timePeriods{yearIndex(i)})
                if k <= size(weightedAK{i,j}, 1)
                    t = (t_cell{i,j}-t_cell{1,1}(1))/24/3600+1;
                    z = z_cell{i,j};
                    Ak = reshape(weightedAK{i,j}(k,:), length(z), length(t));
                    [tmesh, zmesh] = meshgrid(t,z);
                    surf(tmesh, zmesh, Ak)
                    hold on
                    view(2)
                    shading interp;
                    colorbar
                    colormap(jet)
                    axis tight
                    xlabel('t (days)')
                    ylabel('z')
                    grid off
                    ylim([-2,8]);
                    title(['$A_{K',num2str(k),'}$'], 'Interpreter','latex')
                    caxis([-5e-4, 5e-4]);
                    ylim([1,8])
                end
            end
        end
        if saveAk
            print(fig, ['Figures/Sensitivity_AK_Nk', num2str(NK), '_K', num2str(k) ], '-dpng', '-r600');
        end
    end
end

%% print out ARho
if printARho
    for i = 1: length(yearIndex)
        for j = 1: length(timePeriods{yearIndex(i)})
            tempARho = weightedARho{i,j}*1e4;
            disp(['% 201', num2str(yearIndex(i)+1), '-', num2str(j) ]);
            disp(['\[A_{\rho 201', num2str(yearIndex(i)+1), '}=\begin{bmatrix}']);
            for m = 1:size(tempARho,1)
                outputString = sprintf('%.3f & ', tempARho(m,1:end-1));
                outputString = [outputString, sprintf('%.3f \\\\ ', tempARho(m,end))];
                disp(outputString);
            end
            disp('\end{bmatrix}\]');
        end
    end
end
% Script for preprocessiong density data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

%%
fileList = {'LF12.mat', 'LF13.mat', 'LF14.mat', 'LF15.mat'};
rhoData = {};


for i = 1:length(fileList) 
    load(fileList{i});
    
    % Special treatment for LF12
    if contains(fileList{i}, '12')
        rho = rho_rs;
    end
    
    % Convert the coordinates to the middle points
    rawZ = 0.5 * (rho(:,1) + rho(:,2));
    
    % Padding rho=rho(1) at z= 0
    [~, rawRhoIndex0] = min(rho(:,1));
    rhoPadHead = [0, rho(rawRhoIndex0,3)];
    
    % Padding rho=900 at z= 20
    rhoPadTail = [20, 900];
    
    % Padding
    rawRhoUni = [rhoPadHead; rawZ,rho(:,3); rhoPadTail];
    
    % z coordinate
    z = linspace(0,12,500)';
    
    % Linear interpolation
    rhoInterp = interp1(rawRhoUni(:,1), rawRhoUni(:,2), z, 'linear', 'extrap');
    
    % Plot
    subplot(2,2,i)
    plot(rho_reg(:,1),rho_reg(:,2))
    hold on
    plot(z, rhoInterp)
    xlim([0,12]);
    ylim([300, 900]);
    title(fileList{i});
    xlabel('$z (m)$', 'Interpreter', 'latex');
    ylabel('$\rho (kg/m^3)$', 'Interpreter', 'latex');
   
    % Save data
    rhoData{i}.z = z;
    rhoData{i}.rho = rhoInterp;  
end

legend({'$\rho\in P_0$', '$\rho\in P_1$'}, 'Interpreter', 'latex')

%%
fileName = 'densityData.mat';
save(fileName, 'rhoData');

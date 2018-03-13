function [weightedAK, weightedB, weightedSE, weightedAz, weightedD, dTdz, T_data, mask, t_data, z_data] = solveSensitivity(yearIndex, K_opt, zK, timePeriod)
    %% Initialize
    % Settings
    interpOption = 'linear';

    %% Load data
    load('LF_4_aver.mat');
    % Load density data
    load('densityData.mat');
    % Load average and mean
    load('summary.mat');
    
    % Assign data
    data = LF{yearIndex}.T;
    rho = rhoData{yearIndex};

    % Load Measurements
    [t_data, z_data, T_data, indt] = loadData(data, 0, timePeriod);

    % Physical parameters
    C = 152.5 + 7.122 * (273.15 - 10);
    
    %% Initialization
    % Initial Heat conductivity
    zK = zK(:);
    Nk = length(zK);
    K_opt = K_opt(:);

    % cut the data according to the range of K
    [T_data, z_data, indCutZ] = cutData(T_data, z_data, [zK(1),zK(end)]);
    
    % size of measurments
    Nz = length(z_data);

    % Take weights into account
    if ((yearIndex == 1) || (yearIndex == 3))
        T_S = dataS{yearIndex}.T_S(indCutZ, indt);
    elseif (yearIndex == 2)
        T_S = data.T_sd(indCutZ, indt);
        T_S = T_S.^2;
    else 
        T_S = 0.02*ones(size(T_data));
    end
    
    % Cut the time series with Nan in T_data
    [T_data, t_data, noNanInd] = cutNan(T_data, t_data);

    T_S = T_S(: , noNanInd);

    % mask for T >= -2 put 0
    mask = find( T_data > -2);
    
    %% Compute dTdz
    dZfine = 5;
    [dTdz, matDTDz, ~, ~] = computeDTDz(z_data, t_data, T_data, dZfine, zK, K_opt, rho, C, mask, interpOption);
%     [dTdzData, ~, ~, ~] = computeDTDzFromData(z_data, t_data, T_data, dZfine, [], interpOption);
% 
%     t_data = scaleTimeUnit(t_data);
%     [X_data, Y_data] = meshgrid(t_data, z_data);
%     figure
%     subplot(2,1,1)
%     surf(X_data, Y_data, dTdz);
%     view(2)
%     shading interp;
%     colorbar
%     colormap(jet)
%     axis tight
%     xlabel('t (days)');
%     ylabel('z (m)');
%     title(['Temperature gradient 201', num2str(yearIndex+1)])
%     caxis([-4, 4])
%     subplot(2,1,2)
%     surf(X_data, Y_data, dTdzData);
%     view(2)
%     shading interp;
%     colorbar
%     colormap(jet)
%     axis tight
%     xlabel('t (days)');
%     ylabel('z (m)');
%     title(['Temperature gradient 201', num2str(yearIndex+1)])
%     caxis([-4, 4])
    %% Compute A  
    dZfine = 5;
    dK = 1e-6;
    [A, ~, ~, error] = computeKSensitivity(z_data, t_data, T_data, dZfine, zK, K_opt, dK, rho, C, mask, interpOption);

    %%
%     Nt = (numel(T_data)-numel(mask))/ length(z_data);
    Nt = length(t_data);
    %% compute weight
    w = 1./ ((T_S));
    W = w(:);
    W(mask) = 0;
    weightedA = A' * spvardiag(W);
    weightedB = (weightedA * A);
    weightedAK = weightedB \ weightedA;

    weightedSE = sqrt( diag( inv(weightedB./Nt) ) );

%     %% Plot Weighted AK and A
%     figure
%     [X_data, Y_data] = meshgrid(t_data, z_data);
% 
%     for i = 1 : Nk+1
%         if i > 1
%     %         p_data = reshape(weightedA(i-1,:), size(X_data));
%             p_data = reshape(weightedAK(i-1,:), size(X_data));
%             subTitle = ['K',num2str(i-1)];
%         else
%             p_data = reshape(W.^0.5, size(X_data));
%             subTitle = ['Weights in 201',num2str(yearIndex+1)];
%         end
% 
%         p_data(mask) = 0;
% 
%         subplot(ceil((Nk+1)/2), 2, i)
%         surf(X_data, Y_data, p_data);
%         view(2)
%         shading interp;
%         colorbar
%         colormap(jet)
%         axis tight
%         title(subTitle);
%         if i == 1
%             caxis([0, 15]);
%         else
%             caxis([-1e-3, 1e-3]);
%         end
%         grid off
%     end

    %% Plot Az and Az*R

    zALeg = {};
    for i = 1: Nk
        zALeg{i} = ['K', num2str(i)];
    end
    
    R = tril(ones(Nz));

    weightedAz = weightedAK * matDTDz;
    weightedD = weightedAz * R;
    
%     figure
%     subplot(2,1,1)
%     plot(z_data, weightedAz');
%     xlim([min(z_data), max(z_data)])
%     % legend(zALeg);
%     xlabel('z');
%     ylabel('A_z');
%     subplot(2,1,2)
%     plot(z_data, weightedD');
%     xlim([min(z_data), max(z_data)])
%     xlabel('z');
%     ylabel('A_zR');
%     legend(zALeg);


% Sensitivity analysis with the optimal K and Rho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'yearIndex'     - index of the year to be used, from 2012-2015;
%   'K_opt'         - the optimal K solutions;
%   'zK'            - z-coordinate of the K parameter;
%   'timePeriod'    - only part of t_data are taken into account. [0,1] indicates
%                       the full data piece;
%   'rho_opt'       - the optimal rho solutions;
%   'zRho'          - z-coordinate of the rho parameter;    
% The return values:
%   'weightedAK'    - weighted sensitivity matrix dK = weightedAK*dT;
%   'weightedB'     - weighted matrix A^T*W*A;
%   'weightedSE'    - sensitivity of K;
%   'weightedSE_t_indep'    - sensitivity of K, computed by time independent assumption;
%   'weightedAz'    - weighted sensitivity matrix dK = weightedAz*dz;
%   'weightedD'     - accumulated weighted sensitivity matrix dK = weightedD*de, 
%                       where dz = R*de;
%   'weightedARho' 	- weighted sensitivity matrix dk = weightedARho*dRho;
%   'dTdz'        	- temprature gradient;
%   'T_data'        - measurement after cutoff;
%   'mask'          - mask matrix used in sensitivity analysis;
%   't_data'        - time points used in the computation, in the unit of seconds;
%   'z_data'        - spatial grids used in the computation.
%   'A12'           - A^T*W*ARho
%   'A22'           - ARho^T*W*ARho
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2019-03-28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [weightedAK, weightedB, weightedSE, weightedSE_t_indep, weightedAz, ...
    weightedD, weightedARho, dTdz, T_data, mask, t_data, z_data, A12, A22] = ...
    solveSensitivity(yearIndex, K_opt, zK, timePeriod, rho_opt, zRho)
    if nargin < 5
        loadRhoFlag = 1;
    else 
        loadRhoFlag = 0;
    end
    %% Initialize
    % Settings
    interpOption = 'linear';

    %% Load data
    try
        load('LF_4_aver.mat');
    catch
        error('Please check the original data file LF_4_aver.mat');
    end
    % Assign data
    data = LF{yearIndex}.T;

    % Load average and mean
    try
        load('summary.mat');
    catch
        error('summary.mat not found. Try to run averageAndVariance.m first.');
    end
    
    % Load density data
    if loadRhoFlag
        try
            load('densityData.mat');
        catch
            error('densityData.mat not found. Try to run preprocessRho.m first.');
        end
        % Load rho data
        rho = rhoData{yearIndex};
    else
        rho.rho = rho_opt;
        rho.z = zRho;
    end
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

    %% Compute AK  
    dZfine = 5;
    dK = 1e-6;
    [A, ~, ~, ~] = computeKSensitivity(z_data, t_data, T_data, dZfine, zK, K_opt, dK, rho, C, mask, interpOption);

    %% Compute ARho
    dRho = 1e-3;
    ARho = computeRhoSensitivity(z_data, t_data, T_data, dZfine, zK, K_opt, dRho, rho, C, mask, interpOption);

    %% Compute weight
    Nt = length(t_data);
    Nz = length(z_data);
    
    w = 1./ ((T_S));
    W = w(:);
    W(mask) = 0;
    weightedA = A' * spvardiag(W);
    % A^T*W*A
    weightedB = (weightedA * A);
    % Ar
    weightedAK = weightedB \ weightedA;
    % variance
    weightedSE = sqrt( diag( inv(weightedB./Nt) ) );

    % Cov(dk)=Ar*e*s^2*e^T*Ar^T with the assumption that the errors are
    % time independent
    tempTS = T_S(:);
    tempTS(mask) = nan;
    
    sigmaW = nanmean(reshape(tempTS.^0.5,Nz, Nt), 2);
    sigmaW(isnan(sigmaW)) = 0;
    et = speye(Nz);
    eE = repmat(et, Nt, 1);
    E =eE * spvardiag(sigmaW.^2) * eE';
    testS = weightedAK * E * weightedAK';
    weightedSE_t_indep = sqrt(diag(testS));
    
    weightedARho = -weightedAK * ARho;
    %% Plot Az and Az*R
    zALeg = {};
    for i = 1: Nk
        zALeg{i} = ['K', num2str(i)];
    end
    
    R = tril(ones(Nz));

    weightedAz = weightedAK * matDTDz;
    weightedD = weightedAz * R;

    %% prepare for the eigenvalue check
    A12 = weightedA * ARho;
    A22 = ARho' * ARho;

% Solve 1D Heat equation with variable coefficient K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   't'             - grid point of the time discretization;       
%   'z'             - grid point of the spatial discretization;
%   'K'             - variable coefficient K(z);
%   'heatParam'     - other coefficients: 
%                       'dt'        - time step;
%                       'rho'       - density of the snow;
%                       'C'         - heat capacity;
%                       'TbcUp'     - Dirichlet boundary condition at z=0
%                       'TbcDown'   - Dirichlet boundary condition at z=12
%                       'T0'        - intial condition         
% The return values:
%   'T'             - the solution in time*space;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T]=solveHeat(t, z, K, heatParam)
    %% 
    dt = heatParam.dt;
    rho = heatParam.rho;
    C = heatParam.C;
    T0 = heatParam.T0;
    TbcUp = heatParam.TbcUp;
    TbcDown = heatParam.TbcDown;
    
    %% Get Kp on each z
    Kp = HeatConductivity(heatParam.zK, K, z);
    
    %% Initialization
    T_old = T0;
    nt = length(t);
    dx = abs(z(2)-z(1));
    Nx = length(z);
    T = zeros(Nx,nt);
    
    D1m = Dm(Nx, dx, 0);
    D1p = Dp(Nx, dx, 0);
    Kc = spvardiag(Kp);
    D = 1./rho./C.*(D1p*Kc*D1m);

    %% Time iterations
    for i = 1: nt    
        % Trapzoidal Rule
        T_new = AdamsMoulton(T_old, dt, D, [TbcUp(i),TbcDown(i)]);
        
        % Dirichlet B.C.
        T_new(1) = TbcUp(i);
        T_new(end) = TbcDown(i);
        
        % Save and update
        T(:,i) = T_new;
        T_old = T_new;
    end
end
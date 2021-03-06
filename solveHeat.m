% Solve 1D Heat equation with variable coefficient K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   't'             - grid point of the time discretization;       
%   'z'             - grid point of the spatial discretization;
%   'x'             - variable coefficient x including K (and rho)
%   'heatParam'     - other coefficients: 
%                       'dt'        - time step;
%                       'rho'       - density of the snow in vector form;
%                       'zRho'      - depth of the density;
%                       'C'         - heat capacity;
%                       'TbcUp'     - Dirichlet boundary condition at z=0
%                       'TbcDown'   - Dirichlet boundary condition at z=12
%                       'T0'        - intial condition         
% The return values:
%   'T'             - the solution in time*space;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T]=solveHeat(t, z, x, heatParam)
    %% Set up 
    dt = heatParam.dt;
    C = heatParam.C;
    T0 = heatParam.T0;
    Tbc = heatParam.Tbc;
    
    % split K and rho from x if needed
    Nk = length(heatParam.zK);
    if length(x) == Nk
        rho = heatParam.rho;
        K = x;
    else
        K = x(1:Nk);
        % scale rho back to the unit system for solving the heat equation
        rho.rho = x(Nk+1:end) * heatParam.rhoScale;
        rho.z = heatParam.zRho;
    end
    
    %% Initialization
    T_old = T0;
    nt = length(t);
    dx = abs(z(2)-z(1));
    Nx = length(z);
    
    % Get Kp on each z
    zHalf = z - 0.5*dx;
    zHalf(1) = z(1);
    zHalf(end) = z(end);
    Kp = heatConductivity(heatParam.zK, K, zHalf);
    
    % Get rho on each z
    rhop = density(rho, z);

    T = zeros(Nx,nt);
    T(:,1) = T_old;
    
    D1m = Dm(Nx, dx, 0);
    D1p = Dp(Nx, dx, 0);
    Kc = spvardiag(Kp);
    D = 1./rhop./C.*(D1p*Kc*D1m);

    %% Time iterations
    for i = 2: nt    
        % Trapzoidal Rule
        T_new = AdamsMoulton(T_old, dt, D, Tbc, i);
        
        % Dirichlet B.C.
        T_new(1) = Tbc.Up(i);
        T_new(end) = Tbc.Down(i);
        
        % More Dirichlet for the masked area
        dBCInd = Tbc.mask(:, i);
        T_new(dBCInd) = Tbc.data(dBCInd, i);
        
        % Save and update
        T(:,i) = T_new;
        T_old = T_new;
    end
end
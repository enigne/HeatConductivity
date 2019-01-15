% Set numerical and physical parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'dt'            - time step size;
%   'Nt'            - number of time steps;
%   'dz'            - spatial step size;
%   'Nz'            - number of spatial discretization;
%   'rho'           - density;
%   'C'             - heat capacity;
%   'T0'            - initial condition(temperature);
%   'Tbc'           - boundary condition(temperature);
%   'zK'            - z coordinates for K;
%   'zRho'          - z coordinates for Rho;
%   'rhoScale'      - scaling factor of rho such that it remains at the
%                     same level as K in x, but in a correct unit when use it;
%   'maskConst'     - constant threashold for temperature measurements;   
% The return values:
%   'heatParam'     - structure contains parameters for heat equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-15
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function heatParam = setHeatParam(dt, Nt, dz, Nz, rho, C, T0, Tbc, zK, zRho, rhoScale, maskConst)
    if nargin < 12
        maskConst = -2;
        if nargin < 11
            rhoScale = 330;
            if nargin < 10
                zRho = zK;
            end
        end
    end
    
    %%
    heatParam.dt = dt;
    heatParam.Nt = Nt;
    heatParam.dz = dz;
    heatParam.Nz = Nz;
    heatParam.rho = rho;
    heatParam.C = C;
    heatParam.T0 = T0;
    heatParam.Tbc = Tbc;   
    heatParam.zK = zK;   
    heatParam.zRho = zRho;   
    heatParam.rhoScale = rhoScale;
    heatParam.maskConst = maskConst;
end
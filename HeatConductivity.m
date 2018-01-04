% Function to interpolate heat conductivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'zC'            - given coarse grid of KC;
%   'KC'            - given variable coefficient KC;      
%   'z'             - z-coordinate.     
% The return values:
%   'K'             - a piecewise linear function of z;       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K =HeatConductivity(zC, KC, z)
    K = interp1(zC, KC, z,'linear');
end
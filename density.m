% Function to interpolate the density, if rhoC is a constant, then fill in
% the whole vector by it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'zC'            - given coarse grid of KC;
%   'rhoC'          - given variable coefficient KC;      
%   'z'             - z-coordinate.     
% The return values:
%   'rho'           - a piecewise linear function of z;       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rho =density(zC, rhoC, z)
    if (length(rhoC) == 1)
        rho = ones(size(z)) * rhoC;
    else
        rho = interp1(zC, rhoC, z,'linear');
    end
end
% Function to interpolate the density, if rhoC is a constant, then fill in
% the whole vector by it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'rhoC'          - given variable coefficient rhoC.rho and the
%                     coordinates at rhoC.z.
%   'z'             - z-coordinate.     
% The return values:
%   'rho'           - a piecewise linear function of z;       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rho =density(rhoC, z)
    if (isstruct(rhoC))
        rho = interp1(rhoC.z, rhoC.rho, z, 'linear', 'extrap');
    else
        rho = ones(size(z)) * rhoC;
    end
end
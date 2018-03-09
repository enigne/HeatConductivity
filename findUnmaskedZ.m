% Find all the unmasked z values from z_data based on maskTval for T_data,
% which is the largest z-coordinates for the unknowns in the heat equation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'z_data'	- z-coordinates from data;
%   'T_data'	- temperature measurements,
%   'maskTval'	- mask value for imposing Dirichlet conditions from T_data.
% The return values:
%   'umaskedZ'	- largest unmasked z-coordinates in the current data piece.
%   'NmasksZ'	- number of masks at each z coordinate in the whole data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [umaskedZ, NmasksZ] = findUnmaskedZ(T_data, z_data, maskTval)
    if nargin < 3
        maskTval = -2;
    end
    % determine the dimension of the data piece
    Nt = size(T_data, 2);
    
    % take mask
    T_masked = (T_data > maskTval);
    
    % compute number of masks
    NmasksZ = sum(T_masked, 2);
    umaskedZ = z_data(NmasksZ < Nt);
end
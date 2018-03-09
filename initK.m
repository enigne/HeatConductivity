% Initialize K and zK according to the masks and the user setup, remove 
% user's zK from the second one which falls out of the unmaskedZ range.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'zK'            - z-coordinate of the K parameter;
%   'umaskedZ'      - largest unmasked z-coordinates in the current data piece;
%   'constScaleK0'  - scaling factor for K0.
% The return values:
%   'zKmasked'      - zK after the mask;
%   'K0'            - provide initial guess for K with certain unit and same dimension as zKmasked;
%   'zKMflag'       - indicator for zK with mask.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [zKmasked, K0, zKMflag] = initK(zK, umaskedZ, constScaleK0)
    if nargin<3
        constScaleK0 = 1;
    end
    % find the range of umaskedZ
    uMZmax = max(umaskedZ);
    uMZmin = min(umaskedZ);
    
    % compare zK and umaskedZ
    zK = zK(:);
    Nk = length(zK);
    zKMflag = ones(Nk, 1);
    zKMflag = logical(zKMflag);
    
    for i = 1: Nk-1
        % both larger than uMZmax
        if ((zK(i) >= uMZmax) && (zK(i+1)) > uMZmax)
            zKMflag(i+1) = 0;
        end
        % both smaller than uMZMin
        if ((zK(i) < uMZmin) && (zK(i+1)) <= uMZmin)
            zKMflag(i) = 0;
        end
    end
    zKmasked = zK(zKMflag);
    
    Nk = length(zKmasked);
    K0 = constScaleK0*ones(Nk, 1);
end
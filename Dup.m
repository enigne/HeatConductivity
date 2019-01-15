% D1 Operator of first derivative for N-vectors with space size dx 
% upwind scheme: if q> 0 qx = (q(i)-q(i-1))/dx; if q<0, qx = (q(i+1)-q(i))/dx;
% upFlag = 0 for q > 0; upFlag = 1 for q < 0;
% periodFlag = 0 then B.C. is Direchlet(Default), so D(end,end) === 0;
% periodFlag = 1 then B.C. is periodical, so D(end,1) === 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'N'             - the size of the operator;
%   'dx'            - the space size;
%   'upFlag'        - the upwind flag, upFlag = 0 for q > 0; upFlag = 1 for q < 0;
%   'periodFlag'    - periodFlag = 0 then B.C. is Direchlet(Default), so D(end,end) === 0; 
%                     periodFlag = 1 then B.C. is periodical, so D(end,1) === 1
% The return values:
%   'D'             - The D1 Operator.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = Dup(N, dx, upFlag, periodFlag)
    % Set default value for periodFlag
    switch nargin
        case 3
        periodFlag = 0;
    end
    
    % Get the operator only
    D1p = Dp(N, 1, periodFlag);
    D1m = Dm(N, 1, periodFlag);
        
    % Sign matrix
    Pup = spvardiag(upFlag);
    Pdown = spvardiag(ones(N,1)- upFlag);
 
    % Upwind shceme
    D = Pup*D1p + Pdown*D1m;
    
    % Scalling
    D = D/dx;
end


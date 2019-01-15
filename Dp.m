% D+ Operator of first derivative for N-vectors with space size dx 
% one side finite differences qx = (q(i+1)-q(i))/dx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'N'             - the size of the operator;
%   'dx'            - the space size;
%   'periodFlag'    - periodFlag = 0 then B.C. is Direchlet(Default), so D(end,end) === 0; 
%                     periodFlag = 1 then B.C. is periodical, so D(end,1) === 1
% The return values:
%   'D'             - The D+ Operator.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = Dp(N, dx, periodFlag)
    % Set default value for periodFlag
    switch nargin
        case 2
        periodFlag = 0;
    end
    
    % One side finite differences
    I = speye(N);
    D = -I + circshift(I,[0,1]);
    D(end,1) = 0;

    if periodFlag == 0
        % Direhelet 
         D(end,end) = -1;
         D(end,end-1) = 1;
    else
        % Periodical
        D(end,1) = 1;
    end
   
    D = sparse(D)/dx;
end
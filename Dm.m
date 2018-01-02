% D- Operator of first derivative for N-vectors with space size dx 
% one side finite differences qx = (q(i)-q(i-1))/dx
% periodFlag = 0 then B.C. is Direchlet(Default), so D(end,end) === 0
% periodFlag = 1 then B.C. is periodical, so D(end,1) === 1
%
function D = Dm(N, dx, periodFlag)
    % Set default value for periodFlag
    switch nargin
        case 2
        periodFlag = 0;
    end
    
    % One side finite differences
    I = speye(N);
    D = I - circshift(I,[1,0]);
    D(1,end) = 0;

    if periodFlag == 0
        % Direhelet 
         D(1,1) = -1;
         D(1,2) = 1;
    else
        % Periodical
        D(1,end) = -1;
    end
   
    D = sparse(D)/dx;
end
% D0 Operator of first derivative for N-vectors with space size dx 
% Central differences qx = (q(i+1)-q(i-1))/2dx
% periodFlag = 0 then B.C. is Direchlet(Default), so D(end,end) === 0
% periodFlag = 1 then B.C. is periodical, so D(end,1) === 1
%
function D = Dcd(N, dx, periodFlag)
    % Set default value for periodFlag
    switch nargin
        case 2
        periodFlag = 0;
    end

    % Central differences
    I = 0.5*speye(N);
    D = circshift(I,[0,1])-circshift(I,[1,0]);
    D(1,end) = 0;
    D(end,1) = 0;
    
%     v = 0.5*ones(N-1,1);
%     D = sparse(diag(v,1))+sparse(diag(-v,-1));

    % Boundary Condition
    if periodFlag == 0
        % Direhelet 
        D(1,1) = -1;
        D(1,2) = 1;
        D(end,end) = 1;
        D(end,end-1) = -1;    
    else 
        % Periodical
        D(1,end) = -0.5;
        D(end,1) = 0.5;
    end
    D = sparse(D)/dx;
end


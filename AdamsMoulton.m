% Adams-Moulton second order implicit time stepping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'y_old'         - solution at t^n;       
%   'dt'            - time step;       
%   'D'             - spatial discretization operator;
%   'BC'            - boundary conditions.
% The return values:
%   'y_new'         - new solution of at t^{n+1}.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-03-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y_new = AdamsMoulton(y_old, dt, D, BC, tInd, zInd)
    n = size(D, 1);
    % Default values
    if nargin < 6
        zInd = [1, n];
    end
    
    I = speye(n);
    rhs = (I + dt/2*D) *y_old;
    Q = I - dt/2*D;
        
    % normal lower and upper boundaries
    if nargin > 3
        rhs(1) = BC.Up(tInd);
        rhs(n) = BC.Down(tInd);
        Q(1, :) = 0;
        Q(n, :) = 0;
        Q(1, 1) = 1;
        Q(n, n) = 1;
    end
    
    % fill-in all the mask area with Dirichlet B.C.
    dBCInd = BC.mask(:, tInd);
    rhs(dBCInd) = BC.data(dBCInd,tInd);
    Q(dBCInd, : ) = 0;
    Q(dBCInd, dBCInd) = speye(sum(dBCInd));
    
    y_new = Q \ rhs;
end
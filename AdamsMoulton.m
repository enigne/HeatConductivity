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
% Date: 2018-01-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y_new = AdamsMoulton(y_old, dt, D, BC)
    n = size(D, 1);
    I = speye(n);
    rhs = (I + dt/2*D) *y_old;
    Q = I - dt/2*D;
        
    if nargin > 3
        rhs(1) = BC(1);
        rhs(n) = BC(2);
        Q(1, :) = 0;
        Q(n, :) = 0;
        Q(1, 1) = 1;
        Q(n, n) = 1;
    end
    y_new = Q \ rhs;
end
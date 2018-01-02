% Function of heat conductivity
% K is a piecewise function of z
% with KC distributed on zC

function K =HeatConductivity(zC, KC, z)
    K = interp1(zC, KC, z,'linear');
end
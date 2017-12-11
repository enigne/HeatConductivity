function [ C ] = Cice( T )
%Cice returns specific heat capacity [C] of ice at given temperature [T]
% T - ice temperature in [degC]
% C - specific heat capacity in [J/kg/K]
% C is calculated as a function of T based on expression in Cuffey and Paterson, 2010, page 400

% test: T = [-50:0];

C = 152.5+7.122.*(273.15+T);

end


function [ k ] = keffcalonne( r )
%keffsturm( r ) returns the effective thermal conductivity of snow/firn (k) in [J (s m K)^-1] 
% as a function of density (r) in [kg m^-3] following Calonne et al., 2011
k = 0.024 - 1.23*10^-4.*r + 2.5*10^-6.*r.^2;
end
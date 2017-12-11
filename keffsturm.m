function [ k ] = keffsturm( r )
%keffsturm( r ) returns the effective thermal conductivity of snow/firn (k) in [J (s m K)^-1] 
% as a function of density (r) in [kg m^-3] following Sturm et al., 1997
k = 0.138 - 1.01*10^-3.*r + 3.233*10^-6.*r.^2;
end
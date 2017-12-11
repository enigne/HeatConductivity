function [K] = densification_Ligt11(T, rho, Tavg, b)
%densification_Ligt11 estimate gravitational densification
% OUT: K densification rate, [kg m^-3 year^-1]
% IN: layer variables:
%                     T    = [ 273.15+[-5] 273.15+[-6] ];  % temperature [K]
%                     rho  = [500 600];          % density [kg/m^3]
%     profile parameters:
%                     Tavg = [ 273.15+[-5] 273.15+[-6] ];  % temporal mean subsurface temperature [K]
%                     b    = 540;          % accumulation rate [mm we/year]
%     constants
g       = 9.81;      % gravitational acceleration [m/s^2]
rho_ice = 917;       % density of ice [kg/m^3]
R       = 8.314;     % universal gas constant [J/mol/K]
Ec      = 60000;     % activation energy associated with creep by lattice diffusion [J/mol]
Eg      = 42400;     %                                            grain growth

%% calculation
% T    = 273.15+D{2}.Traw(:,t-1);
% rho  = D{2}.R( :,t );
% Tavg = 273.15+D{2}.Tm;
% b    = 660;


case1 = rho <  550;
case2 = rho >= 550;

C(case1) = 0.0991 - 0.0103*log(b);  % 0.07;
C(case2) = 0.0701 - 0.0086*log(b);  % 0.03;

K = C'.*b.*g.*(rho_ice-rho).*exp(-Ec./R./T+Eg./R./Tavg);
% K = K*3600*24*365.25;

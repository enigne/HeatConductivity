function [R] = Rdensity(NRA)
homedir = '';
% homedir = 'C:\DATA\Cores\';
% NRA = 100;
R.z = [ 0:0.1:25 ]';    % define the vertical model domain
R.NRA = NRA;
% load density data
load([homedir 'LF12.mat'], 'rho_rs_fill'); r{1} = rho_rs_fill;
load([homedir 'LF13.mat'], 'rho'        ); r{2} = rho;
load([homedir 'LF14.mat'], 'rho'        ); r{3} = rho;
load([homedir 'LF15.mat'], 'rho'        ); r{4} = rho;
clear rho_rs_fill rho homedir

for y = 1:4
    rn{y} = r{y};
    rn{y}(:,3) = normrnd(r{y}(:,3), NRA/3); %  figure; subplot(1,2,1); plot(NRA*randn( 1000,1 ), '.'); subplot(1,2,2); plot(normrnd(0, NRA, [1000 1]), '.')
    ind = find( rn{y}(:,3)>900 );
    rn{y}(ind,3) = 900; clear ind

    tmp = NaN(2*size(rn{y},1),2);
    tmp(1:2:length(tmp(:,1)),:) = [rn{y}(:,1)+0.0001 rn{y}(:,3)];
    tmp(2:2:length(tmp(:,1)),:) = [rn{y}(:,2)-0.0001 rn{y}(:,3)];

    % interpolate to a regular grid, extrapolate to have 900 kg/m^3 at 20 m
    [~, up  ] = min( abs( R.z - tmp( 1 ,1) ) );
    [~, dn  ] = min( abs( R.z - tmp(end,1) ) );
    [~, dn20] = min( abs( R.z - 20         ) );
    R.R{y}( 1     : up-1        ) = tmp(1,2);
    R.R{y}(up     : dn          ) = interp1( tmp(:,1)   , tmp(:,2)        , R.z(up:dn  ), 'linear', 'extrap' );
    R.R{y}(dn     : dn20        ) = interp1([R.z(dn) 20], [R.R{y}(dn) 900], R.z(dn:dn20), 'linear', 'extrap' );
    R.R{y}(dn20+1 : size(R.z,1) ) = 900;
    R.R{y} = R.R{y}';

    R.Rs{y} = smooth(R.R{y}, 10, 'moving'); % smoothing with running mean

end; clear y up dn dn20 tmp

% close all;
% figure; hold on; set(gca, 'YDir', 'reverse')
% c = lines(4);
% for y = 1:4
%     plot( R.R {y}  , R.z      , 'Color', c(y,:), 'LineWidth', 2, 'DisplayName', [num2str(2011+y) ' inter- and extra- polated'])
%     plot( R.Rs{y}  , R.z      , 'Color', c(y,:), 'LineWidth', 2, 'DisplayName', [num2str(2011+y) ' inter- and extra- polated, smoothed'])
%     for l = 1:size(r{y},1)
%         plot( [ r{y}(l,3)  r{y}(l,3)], [ r{y}(l,1)  r{y}(l,2)]/100, 'Color', c(y,:), 'DisplayName', [num2str(2011+y) ' raw']);
% %         plot( [rn{y}(l,3) rn{y}(l,3)], [rn{y}(l,1) rn{y}(l,2)]/100, 'Color', c(y,:), 'DisplayName', [num2str(2011+y) ' noizy'], 'LineWidth', 2);
%         if l == 1; legend show; end
%     end; clear l
%     [mean(R.R {y}) mean(R.Rs{y})]
% end;
% clear c y r rn

% save('C:\FILES\_PHD\4\R100.mat')
end

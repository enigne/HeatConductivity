clear; clc;
param{1} = [   0.5     0.2    0.1   0.05 ];    % 3*std in temperature measurements, degC (random noise)
param{2} = [   0.5     0.2    0.1   0.05 ];    % 3*std in temperature measurements, degC (systematic noise)
param{3} = [   1.02    1.05   1.1        ];    %          tortuosity of the cable, m
param{4} = [ 200     100     50    25    ];    %          density measurements, kg m-^3 (random noise)
param{5} = [ 200     100     50    25    ];    %          density measurements, kg m-^3 (systematic noise)

%%
for par  = 3%size(param,2)                  % choose parameter to explore: temperature (1), depth (2), density (3)
for unv  = 2:3%:size(param{par},2)            % uncertainty value
for itn  = 1                             % iterations
clearvars -except param par unv itn
% par = 3; unv = 1; itn = 1;
% NTA = 0; NTAS = 0; NZA = 0; NRA = 0; NRAS = 0;
rng('shuffle');
NTA  = param{par}(unv)*double(par==1);
NTAS = param{par}(unv)*double(par==2);
NZA  = param{par}(unv)*double(par==3);
NRA  = param{par}(unv)*double(par==4);
NRAS = param{par}(unv)*double(par==5);

%% load subsurface temperature data
load('LF_4_aver.mat');
% load('C:\DATA\temperature\LF_4_aver.mat');
for y = 1:4
    if   y == 1; Y = y  ;
    else         Y = y+1;    end
    D.t     {y} = LF{Y}.T.t';
    D.z_a   {y} = LF{Y}.T.z_a;
    D.T_a   {y} = LF{Y}.T.T_a;
    D.T_corr{y} = LF{Y}.T.T_corr;
    D.z_corr{y} = LF{Y}.T.z_corr;
end; clear y Y LF
D.T_a{2}(end-1:end, 25) = 0;

%% track density changes due to gravitational densification

% densification effects taken into account for each year: normally [2 2 1 1]: densification occurs always but in 2015 and 2016 the string was encased in plastic tubes
%     0 - density does not change in time;
%     1 - density increases with time, layers become thinner and advect upwards but thermistors stay in place;
%     2 - density increases with time, layers become thinner and advect upwards together with thermistors.
dr_opt = [ 0 0 0 0 ];

R = Rdensity(NRA);           % load raw density data, add noise (the argument), resample to an even grid and extrapolate
R.R (2) = []; R.Rs(2) = [];  % remove data from 2013

% figure; hold on; plot(R.Rs{1}, -R.z, 'k')
for y = 1:size(R.R,2)   % systematic noise on density profile: sine pattern with amplitude NRAS, period 3 m and random phase offset.
R.Rs{y} = R.Rs{y} + NRAS*sin(2*pi/3*(R.z+normrnd(0, 20)));
R.Rs{y}(R.Rs{y}>910) = 910;
end; clear y
% plot(R.Rs{1}, -R.z, 'k.'); ylim([-12 0])

% average temperature evolution during one complete year
Tmz = D.z_a{2};

n = datevec(D.t{2}( 1 )) + [1 0 0 0 0 0];            % April 2014 - April 2015: the mean is weighted with time intervals between measurements as they are not even in time
n = datenum(n) - D.t{2}(end);
n = round(n*24);
Traw = [D.T_a{2} repmat(D.T_a{2}(:,end), 1, n )];
traw = [D.t{2}   D.t{2}(end)+cumsum( repmat(datenum([0 0 0 1 0 0]), 1, n ) )  ];
Tm14 = wmean( Traw(:,1:end-1), repmat( diff(traw), size(Traw,1), 1), 2 ); clear n Traw traw

n = datevec(D.t{3}( 1 )) + [1 0 0 0 0 0];            % April 2015 - April 2016
n = datenum(n) - D.t{3}(end);
n = round(n*24);
Traw = [D.T_a{3} repmat(D.T_a{3}(:,end), 1, n )];
Traw(end+1:end+2,:) = [Traw(end,:); Traw(end,:)];
Tm15 = mean(Traw,2);    clear n Traw traw
Tm1415 = mean([Tm14 Tm15],2);
% figure; hold on; plot([Tm14 Tm15 Tm1415],-Tmz); title('Mean annual subsurface temperature profile: 2014 (blue), 2015 (red) and mean (yellow)')

[~, indt] = min( abs( D.t{4}(1) - datenum([1 0 0 0 0 0]) - D.t{3} ) ); indt = indt-1;
[~, indz] = min( abs( D.z_a{3}                           - 1.8    ) ); indz = indz+1;
n = ( D.t{3}(1) - D.t{3}(end) + datenum([1 0 0 0 0 0]) ) * 24-1;
ins = D.T_a{3}(indz:end,1:indt);
ins = [ins; zeros(indz-1, size(ins,2))];
D.T_a{4} = [ repmat(ins(:,1),1,round(n)) ins D.T_a{4} ]; clear indt indz ins n
D.t{4} = [ D.t{3}(end)+datenum([0 0 0 1 0 0]):...
                       datenum([0 0 0 1 0 0]):...
           D.t{4}( 1 )-datenum([0 0 0 1 0 0]) D.t{4}];
% figure; hold on; uimagesc(D.t{3}, D.z_a{3}    , D.T_a{3});
%                  uimagesc(D.t{4}, D.z_a{4}-1.8, D.T_a{4});
% set(gca, 'YDir', 'reverse'); datetickzoom('x', 'yy-mmm-dd'); axis tight;

for y = 1:4

    % add random and systematic noise to temperature measurements and associated depths
    for is = 1:size( D.T_corr{y},2 )
                        % temperature data + random noise + systematic noise
        D.T_corr{y}{is} = normrnd(D.T_corr{y}{is}, NTA/3) + normrnd(zeros(size(D.T_corr{y}{is},1),1), NTAS/3);
%         D.z_corr{y}{is} = D.z_corr{y}{is} - 0*[0:0.1:0.1*( size(D.z_corr{y}{is},1)-1 )]'; % depth data + systematic noise

        % random negative noise cumulative with depth: 0 for the uppermost sensor, triangular PDF for the others
        if NZA>0
%              noise = [ 0; random(makedist('Triangular', 'a', -NZA, 'b', 0, 'c', 0), size(D.z_corr{y}{is}) - [1 0]) ];
             NZAt = 1-1/NZA;
             noise = [0; -diff(D.z_corr{y}{is}) * NZAt];
%              noise = [0; -diff(D.z_corr{y}{is}) .* abs( normrnd(zeros(size(D.z_corr{y}{is},1)-1,1), NZA/3/100) )];
        else
             noise = zeros(size(D.z_corr{y}{is}));
        end
%         if y >= 3; noise(noise>0.12) = 0.12;                                    % not larger than 12 cm for 2015 data
%         else       noise(noise>0.24) = 0.24; end                                %                 24 cm for 2012 and 2014 data
        D.z_corr{y}{is} =  D.z_corr{y}{is} + cumsum(noise); clear noise
        D.z_corr{y}{is} = [D.z_corr{y}{is} nan(size(D.T_corr{y}{is}) - [0 1])]; % add nans to allow for time evolution
        D.z_corr{y}{is}(1,2:end) = D.z_corr{y}{is}(1,1);                        % set the uppermost sensor depth constant in time

    end; clear is

    D.R{y} = nan( size(D.z_a{y},1), size(D.t{y},2));
    if     y == 1;  D.R{y}(:, 1 ) = interp1(R.z           , R.Rs{y  }       , D.z_a{y}, 'spline' );  % strings installed in April 2012
    elseif y == 2;  D.R{y}(:, 1 ) = interp1(R.z           , R.Rs{y  }       , D.z_a{y}, 'spline' );  %                            2014
                    D.R{y}(:,end) = interp1(R.z-0.9       , R.Rs{y+1}       , D.z_a{y}, 'spline' );  % (use data from the core drilled in April 2015 as that one has the effects refreezing during summer 2014)
    elseif y == 3;  D.R{y}(:, 1 ) = interp1(R.z           , R.Rs{y  }       , D.z_a{y}, 'spline' );  %                            2015
    elseif y == 4;  D.R{y}(:, 1 ) = interp1(D.z_a{y-1}+1.8, D.R {y-1}(:,end), D.z_a{y}, 'spline' );  %                            2016
    end

    if  dr_opt(y)==0    % density does not change in time

        for is = 1:size( D.T_corr{y},2 )
             D.z_corr{y}{is} = repmat( D.z_corr{y}{is}(:,1), 1, size(D.t{y},2));
        end; clear is

        if     y == 1; D.R{y} = repmat(D.R{y  }(:,1), 1, size(D.t{y},2));
        elseif y == 2; [~, ind] = min( abs(D.t{y} - datenum([2014 9 1 0 0 0])) );
                       D.R{y}(:, 1   :ind) = repmat(D.R{y}(:,1  ), 1, ind);
                       D.R{y}(:,ind+1:end) = repmat(D.R{y}(:,end), 1, size(D.t{y},2) - ind); clear ind
        elseif y == 3; D.R{y} = repmat(D.R{y  }(:,1), 1, size(D.t{y},2));
        elseif y == 4; D.R{y} = repmat(D.R{y-1}(:,1), 1, size(D.t{y},2));
        end

    else  % density increases with time, layers become thinner and advect upwards together with thermistors

        % gravitational densification
        Tm = interp1(Tmz, Tm1415, D.z_a{y}, 'linear', 'extrap');  % interpolate mean annual temperature data to the depth grid of temperature data
        if     any(y == [1 3 4]) % strings installed in 2012, 2015 and 2016: forward densification
               for it = 2:size(D.R{y},2)
                   D.R{y}( :,it ) = D.R{y}( :,it-1 ) + 1*densification_Ligt11( 273.15+D.T_a{y}(:,it-1), D.R{y}(:,it-1), 273.15+Tm, 660 )*( D.t{y}(it) - D.t{y}(it-1) )/365.25;
               end; clear it
        elseif y == 2            % string installed in 2014: forward densification before 2014 September 1 and backwards after
               [~, ind] = min( abs(D.t{y} - datenum([2014 9 1 0 0 0])) );
               for it = 2:ind
                   D.R{y}( :,it ) = D.R{y}( :,it-1 ) + 1*densification_Ligt11( 273.15+D.T_a{y}(:,it-1), D.R{y}(:,it-1), 273.15+Tm, 660 )*( D.t{y}(it) - D.t{y}(it-1) )/365.25;
               end; clear it
               for it = size(D.R{y},2)-1:-1:ind+1
                   D.R{y}( :,it ) = D.R{y}( :,it+1 ) - 1*densification_Ligt11( 273.15+D.T_a{y}(:,it+1), D.R{y}(:,it+1), 273.15+Tm, 660 )*( D.t{y}(it+1) - D.t{y}(it) )/365.25;
               end; clear it ind
        end

        % add 30 kg m^-3 to density to account for the effect of melt water refreezing in summer 2015
        if y == 3
           [~, indt] = min( abs( D.t{y} - datenum([2015 9 1 0 0 0]) ) );
           D.R{y}(:,indt:end) = D.R{y}(:,indt:end) + 30;
        end

        % add 30 kg m^-3 to density to account for the effect of melt water refreezing in summer 2016
        if y == 4
            D.R{4} = D.R{4} + 30;   
        end

        % resample density data to a constant grid and track changes in depth of individual layers (sensors) driven by settling
        ztmp = [ D.z_a{y} nan( size(D.z_a{y},1), size(D.t{y},2)-1 ) ];         % updated depth grid for resampling density data
        ztmp(1,:) = ztmp(1,1);
        for it = 2:size(D.t{y},2)
            th  = diff(ztmp( :, it-1 ));                                       % thicknesses of layers before densification
            thn = th .* D.R{y}( 2:end, it-1 ) ./ D.R{y}( 2:end, it);   %                               after

            if or( [y == 2 && abs(D.t{y}(it)-datenum([2014 9 1 3 0 0]))<0.01], [y == 3 && abs(D.t{y}(it)-datenum([2015 9 1 0 0 0]))<0.01])
            thn = th;             % leave layer thickness unchanged during sudden artificially introduced density changes
%                     D.z_corr{y}{is}(iz, it) = D.z_corr{y}{is}(iz, it-1);
            end
            ztmp( 2:end, it ) = ztmp( 1, it ) + cumsum(thn);

            for is = 1:size(D.z_corr{y},2)
                if dr_opt(y)==1 % prevent thermistors from advecting upwards
                   D.z_corr{y}{is}(2:end, it) = D.z_corr{y}{is}(2:end, it-1);    continue;
                end

                tmp1 = D.z_corr{y}{is}(1 , 1   );                   % depth of the topmost sensor
                [~, ind1] = min( abs( tmp1 - ztmp(:,it-1) ) );      % index of the layer closest to the topmost sensor
            for iz = 2:size(D.z_corr{y}{is},1)
                tmp2 = D.z_corr{y}{is}(iz, it-1);                   % depth of the sensor in question
                [~, ind2] = min( abs( tmp2 - ztmp(:,it-1) ) );      % index of the layer closest to the layer's bottom
                thinning = ( sum( th (ind1:ind2-1) ) - ...          % relative thinning: decrease in summed thicknesses divided by the whole thickness
                             sum( thn(ind1:ind2-1) )   ) / ( ztmp(ind2, it-1) - ...
                                                             ztmp(ind1, it-1) );
                D.z_corr{y}{is}(iz, it) = D.z_corr{y}{is}(iz, it-1) - thinning*(tmp2-tmp1);
            end
            end; clear is iz tmp* ind* thinning

        end; clear it th thn

        if     any(y ==[1 3 4])
               for it = 2:size(D.t{y},2)                     % resample density to a depth grid that is constant in time
                    D.R{y}(:,it) = interp1( ztmp( : , it ), D.R{y}( : ,it), D.z_a{y}, 'linear', 'extrap' );
               end; clear it ztmp
        elseif y == 2            % string installed in 2014: upwards advection of density layers before 2014 September 1 and downwatds after
               [~, ind] = min( abs(D.t{y} - datenum([2014 9 1 0 0 0])) );
               for it = 2:ind
                    D.R{y}(:,it) = interp1( ztmp( : , it ), D.R{y}( : ,it), D.z_a{y}, 'linear', 'extrap' );
               end; clear it
               ztmp( :, end ) = ztmp( :, 1 );
               for it = size(D.R{y},2)-1:-1:ind+1
                   th  = diff(ztmp( :, it+1 ));                                       % thicknesses of layers before densification
                   thn = th .* D.R{y}( 2:end, it+1 ) ./ D.R{y}( 2:end, it);   %                               after
                   ztmp( 2:end, it ) = ztmp( 1, it ) + cumsum(thn);
               end
               for it = size(D.R{y},2)-1:-1:ind+1
                    D.R{y}(:,it) = interp1( ztmp( : , it ), D.R{y}( : ,it), D.z_a{y}, 'linear', 'extrap' );
               end; clear it ind ztmp
        end
    end

    % remove data before 2016 Aug 16, when subsurface temperature measurements are not available
    if y == 4
        [~, indt] = min( abs( D.t{4} - datenum([2016 8 16 16 0 0]) ) );
        D.t     {4}   (:,1 : indt-1) = [];
        D.R     {4}   (:,1 : indt-1) = [];
        D.z_corr{4}{1}(:,1 : indt-1) = [];
    end

    disp(y)

end; clear y Tm* R dr_opt

% close all; figure; hold on; box on; set(gca, 'YDir', 'reverse');
% pcolor(D.t{1}, D.z_a{1}+1.8+0.9+2.4, D.R{1});
% pcolor(D.t{2}, D.z_a{2}+1.8+0.9    , D.R{2});
% pcolor(D.t{3}, D.z_a{3}+1.8        , D.R{3});
% pcolor(D.t{4}, D.z_a{4}            , D.R{4}); shading flat; datetickzoom('x', 'yy-mmm-dd'); colormap(jet(20)); colorbar;

%% interpolate and average
for y  = 1:4
    % interpolation
    for is = 1:size(D.z_corr{y},2)
        D.z_i{y}{is} = [ ceil( D.z_corr{y}{is}( 1 )*10 )/10 : 0.1 : ...
                        floor( D.z_corr{y}{is}(end)*10 )/10]';
        for it = 1:size(D.t{y},2)
            x = D.z_corr{y}{is}(:,it);
            T = D.T_corr{y}{is}(:,it);
            ind = find (isnan(T));
            T    (ind) = [];
            x    (ind) = [];

            if   numel(T) >= 2
                 Ti = interp1(x, T, D.z_i{y}{is}, 'pchip');
                 Ti(Ti==5) = NaN;
            else
                 disp('Not enough data points for interpolation!')
                 Ti = nan(size( D.z_i{y}{is} ));
            end
            D.T_i{y}{is}(:,it) = Ti;
        end; clear it x T ind Ti

    end; clear is

    % averaging
    for is = 1:size(D.z_corr{y},2)
        dn(is) = D.z_i{y}{is}(end);         % vertical grid 
    end
    D.z_r{y} = [1:0.1:max(dn)]'; if y == 4; D.z_r{y} = [D.z_i{y}{is}:0.1:max(dn)]'; end
    clear is dn

    tmp = nan(size( D.z_r{y},1), size(D.T_corr{y},2), size(D.t{y},2 )); % matrix to be filled by data from 9 individual t-strings
    for is = 1:size(D.z_corr{y},2)
        [~, b] = min( abs( 1        - D.z_i{y}{is}      ) );   % 
        [~, f] = min( abs( D.z_r{y} - D.z_i{y}{is}(end) ) );   % index of the lower element             
        tmp(1:f,is,:) = D.T_i{y}{is}(b:end,:);
    end; clear is b f

    D.T_r{y} = squeeze( nanmean(tmp,  2) );
%     D.T_s{y} = squeeze( std    (tmp,1,2) );
    clear tmp

    [~, ind1] = min( abs( D.z_a{y} - D.z_r{y}( 1 ) ) );
    [~, ind2] = min( abs( D.z_a{y} - D.z_r{y}(end) ) );
    D.R_r{y} = D.R{y}(ind1:ind2,:); clear ind*

end; clear y

% close all;
% figure('units','normalized','outerposition',[0 0 1 1]); hold on; grid on;
% y = 2;
% c =  jet(size(D.T_corr{y},2));
% % datevec(D.t{y});
% for it = 2309:6:size(D.t{y},2)
%     for is = 1:size(D.T_corr{y},2)
%         plot(D.T_corr{y}{is}(:,it), D.z_corr{y}{is}(:,1 ), 'o' , 'Color', c(is,:), 'DisplayName', ['raw ' num2str(is)]         )
%         plot(D.T_corr{y}{is}(:,it), D.z_corr{y}{is}(:,it), '+' , 'Color', c(is,:), 'DisplayName', ['densified ' num2str(is)]   )
% %         plot(D.T_i   {y}{is}(:,it), D.z_i   {y}{is}      , '.-', 'Color', c(is,:), 'DisplayName', ['interpolated ' num2str(is)])
%     end;
%     plot(D.T_r{y}(:,it), D.z_r{y}, '.-k', 'LineWidth', 2, 'DisplayName', 'averaged' )
%     legend show
% %     plot(-[0.1 0.1],[0 13],'k')
%     title(['Temperature distribution on ' datestr(D.t{y}(it), 'yyyy-mmm-dd HH') ', time step ' num2str(it)])
%     set(gca, 'YDir', 'reverse')
%     xlim([-10 1]); ylim([0 15]) % 12.3
%     legend show
%     pause%(0.1)
%     cla
% end; clear it is y

% close all;
% figure('units','normalized','outerposition',[0 0 1 1]); hold on; grid on;
% y = 2;
% c = jet( 15 );
% for s = 1:9
% % s = 2;
% n = size(D.T_corr{y}{s},1);
% for iz = 1:n
% plot( D.t{y}, D.T_corr{y}{s}(iz,:), '.-', 'Color', c(iz,:))
% end; clear is n
% end; clear s y
% datetickzoom('x', 'yy-mmm-dd')

D.z = D.z_r;
D.T = D.T_r;
D.R = D.R_r;
D = rmfield(D, {'z_a', 'T_a', 'T_corr', 'z_corr', 'z_i', 'T_i', 'z_r', 'T_r', 'R_r'});
save('D.mat', 'D')
% save('D_sens.mat', 'D')

%% define the time domains for the 6 datasets
clear D      L;       load('D.mat')
% clear D_sens L_sens;  load('D_sens.mat');
yind = [1 2 2 3 3 4];

tmp = [ datenum([2012  4 26  0 0 0]) datenum([2012 7 02 12 0 0]);... % spring 2012
		datenum([2014  4 27  0 0 0]) datenum([2014 7 04 23 0 0]);... % spring 2014
		datenum([2014 12  1  3 0 0]) datenum([2015 4 11 14 0 0]);... % fall   2014: [2014  9 25  3 0 0]
		datenum([2015  4 15  9 0 0]) datenum([2015 7 09  9 0 0]);... % spring 2015
		datenum([2015 10 09  1 0 0]) datenum([2016 4 12  0 0 0]);... % fall   2015
		datenum([2016 10 26 11 0 0]) datenum([2017 4 23 13 0 0]) ];  % fall   2016

for i = 1:6
	[~, ind(i,1)] = min( abs( D.t{yind(i)} - tmp(i,1) ) );
	[~, ind(i,2)] = min( abs( D.t{yind(i)} - tmp(i,2) ) );
	tind(i,:) = ind(i,:);
end; clear tmp i ind

% depth corrections for accumulation: 2.4 m is from stake measurements at S11,
%                                     0.9 m is found by visual fitting of the temperature profiles in April 2015 measured by the old and new setup (see the plotting script below),
%                                     1.8 m is measurements of the distance from SR50 to the surface in April 2016
% off = [ 1.8+0.9+2.4    1.8+0.9    1.8+0.9    1.8    1.8    0 ];
off = [ 0              0          0.9          0      0      0 ];

for y = 1:6
    L.t  {y} = D.t{ yind(y) }(   tind(y,1):tind(y,2) );  L.t_r{y} = L.t{y};           % pick up time        data
    L.z  {y} = D.z{ yind(y) } + off(y);                  L.z_r{y} = L.z{y};           %         depth
    L.T  {y} = D.T{ yind(y) }(:, tind(y,1):tind(y,2) );  L.T_r{y} = L.T{y};           %         temperature
    L.R  {y} = D.R{ yind(y) }(:, tind(y,1):tind(y,2) );  L.R_r{y} = L.R{y};           %         density

    if or(NTA>0, NTAS>0)
        tmp = load('O.mat');
        [~, t1] = min( abs(L.t{y} -  tmp.L.t{y}( 1 )) );
        [~, t2] = min( abs(L.t{y} -  tmp.L.t{y}(end)) );
        [~, z1] = min( abs(L.z{y} -  tmp.L.z{y}( 1 )) );
        [~, z2] = min( abs(L.z{y} -  tmp.L.z{y}(end)) );
        L.t{y} = L.t{y}(  :  ,t1:t2);
        L.z{y} = L.z{y}(z1:z2,  :  );
        L.T{y} = L.T{y}(z1:z2,t1:t2);
        L.R{y} = L.R{y}(z1:z2,t1:t2);

        L.T{ y }( isnan(tmp.L.T{y})  ) = NaN;  % remove temperature datapoints that are set to NaN in the reference dataset
    else
        L.T{ y }( L.T{ y } >= -2 ) = NaN;  % remove temperature datapoints with T>-1 degC

        ind = find( sum(isnan(L.T{y}),2) == size(L.T{y},2) ); % eliminate depth levels with only NaNs
        L.z{y}(ind, :) = [];
        L.T{y}(ind, :) = [];
        L.R{y}(ind, :) = []; clear ind

        ind = find( sum( isfinite(L.T{y}) ) < 3 );   % eliminate time steps with <3 values
        L.t{y}(:, ind) = [];
        L.T{y}(:, ind) = [];
        L.R{y}(:, ind) = []; clear ind
    end

end; clear y off tmp

for y = 1:6
    tmp = isfinite(L.T{y});                               % eliminate datapoints sticking out
           up{y}(size(L.t{y},2)  ) = min( find( tmp(:,size(L.t{y},2)) ) );
           dn{y}(size(L.t{y},2)  ) = max( find( tmp(:,size(L.t{y},2)) ) );
    for it =  size(L.t{y},2)-1:-1:1
           up{y}(       it       ) = min( find( tmp(:,      it      ) ) );
           dn{y}(       it       ) = max( find( tmp(:,      it      ) ) );
        if up{y}(       it       ) > up{y}(it+1)
           up{y}(       it       ) = up{y}(it+1);
        end
        if dn{y}(       it       ) > dn{y}(it+1)
           dn{y}(       it       ) = dn{y}(it+1);
        end
        L.T{y}( 1:up{y}(it)-1    ,it) = NaN;
        L.T{y}(   dn{y}(it)+1:end,it) = NaN;
    end; clear it tmp

    L.T_m{y} = L.T{y};

    for it = 2:size(L.t{y}, 2)-1
        L.T_m{y}( up{y}(it+1)+1:dn{y}(it-1)-1, it  ) = NaN;
    end; clear it
        L.T_m{y}( up{y}(end )+1:dn{y}(end )-1, end ) = NaN;
end; clear y

% close all; figure; hold on;
% plot( L.T{3}(:,end), L.z{3}, '.', 'DisplayName', num2str([datevec(L.t{3}(end))]) )
% plot( L.T{4}(:, 1 ), L.z{4}, '.', 'DisplayName', num2str([datevec(L.t{4}( 1 ))]) )
% legend show; set(gca, 'YDir', 'reverse'); axis tight
% z1 = L.z{3}(isfinite(L.T{3}(:,end))     ); t1 = L.T{3}(isfinite(L.T{3}(:,end)), end);
% z2 = L.z{4}(isfinite(L.T{4}(:, 1 ))     ); t2 = L.T{4}(isfinite(L.T{4}(:, 1 )),  1 );
% [ ~, i1] = min( abs(z2 - z1( 1 )) );
% [ ~, i2] = min( abs(z2 - z1(end)) );
% sqrt( mean( (t2(i1:i2) - t1).^2 ) )
% clear z1 t1 z2 t2 i1 i2 ans

clear up dn
save('O.mat', 'L')
% save('O_sens.mat', 'L')

%%
% % clear; clc;
% load('O.mat')
% % load('O_sens.mat')
% close all;
% figure('units','normalized','outerposition',[0 0 1 1]);
% ax1 = subplot(2,1,1); hold on;
% ax2 = subplot(2,1,2); hold on
% 
% H{1} = [          1    2.9  4.9   6.9  8.9]';
% H{2} = [     1    3.3  5.3  7.3   9.8     ]';
% H{3} = [1    2.2  4.2  6.2  8.2  10.4     ]';
% H{4} = [1    2.2  4.2  6.2  8.2  10.4     ]';
% H{5} = [1    2.2  4.2  6.2                ]';
% H{6} = [2.3  4    6.6                     ]';
% 
% H{2} = [     1    3.3  5.3  7.3 ]';
% H{3} = [     2.2  4.2  6.2      ]';
% H{4} = [1    2.2  4.2  6.2  8.2 ]';
% 
% 
% off = [ 1.8+0.9+2.4    1.8+0.9    1.8    1.8    1.8    0 ];
% % off = [ 0              0          0      0      0      0 ];
% for y = 1:6
%     pcolor(ax1, L.t_r{y}, L.z_r{y}+off(y), L.T_r{y} )
%     pcolor(ax1, L.t  {y}, L.z  {y}+off(y), L.T  {y} );
% %     pcolor(ax1, L.t  {y}, L.z  {y}+off(y), L.T_m{y} )
%     colormap(jet(20))
%     plot  (ax1, L.t_r{y}, frfr(L.T_r{y}, L.z_r{y}+off(y), -0.05), 'k')
%     plot  (ax1, L.t_r{y}, frfr(L.T_r{y}, L.z_r{y}+off(y), -1), 'k')
%     plot  (ax1, L.t_r{y}, frfr(L.T_r{y}, L.z_r{y}+off(y), -2), 'k', 'LineWidth', 2)
%     
%     pcolor(ax2, L.t_r{y}, L.z_r{y}+off(y), L.R_r{y} )
%     pcolor(ax2, L.t  {y}, L.z  {y}+off(y), L.R  {y} )
%     
%     plot(ax1, repmat(L.t{y}(1), size(H{y},1), 1), H{y}+off(y), 'k>', 'MarkerFaceColor', lines(1))
% end; clear y
% 
% set(ax1, 'YDir', 'reverse', 'XTick', [] )
% set(ax2, 'YDir', 'reverse', 'XTick', datenum([2013 4 1 0 0 0; 2014 4 1 0 0 0; 2015 4 1 0 0 0; 2016 4 1 0 0 0; 2017 4 1 0 0 0]));
% datetickzoom(ax2, 'x', 'yy-mmm-dd', 'keepticks');
% caxis   (ax1, [-14 0]); caxis   (ax2, [0 900]);
% shading (ax1, 'flat'  ); shading (ax2, 'flat'  );
% axis    (ax1, 'tight' ); axis    (ax2, 'tight' );
% colorbar(ax1          ); colorbar(ax2          );
% linkaxes([ax1,ax2],'x')
% clear ax1 ax2
 
%% find optimal K distribution
% clear;
load('O.mat')
% load('O_sens.mat')

H{1} = [          1    2.9  4.9   6.9  8.9]';
H{2} = [     1    3.3  5.3  7.3   9.8     ]';
H{3} = [1    2.2  4.2  6.2  8.2  10.4     ]';
H{4} = [1    2.2  4.2  6.2  8.2  10.4     ]';
H{5} = [1    2.2  4.2  6.2                ]';
H{6} = [2.3  4    6.6                     ]';

H{2} = [     1    3.3  5.3  7.3 ]';
H{3} = [     2.2  4.2  6.2      ]';
H{4} = [1    2.2  4.2  6.2  8.2 ]';

% H{1} = [0 100]; H{2} = H{1}; H{3} = H{1}; H{4} = H{1}; H{5} = H{1}; H{6} = H{1};


for y  = [2 3 4]

        H{y}(1:end-1,2) = H{y}(2:end,1); H{y}(end,:) = [];
        H{y}( 1 ,1) = L.z{y}(1);
        H{y}(end,2) = L.z{y}(end);
        tau  = L.t  {y};
        Z    = L.z  {y};
        Texp = L.T  {y};
        T    = L.T_m{y};
        R    = L.R  {y};
        h    = H    {y};

    if y == 9%4
        H{y-1}(1:end-1,2) = H{y-1}(2:end,1); H{y-1}(end,:) = [];
        H{y-1}( 1 ,1) = L.z{y-1}(1);   
        H{y-1}(end,2) = L.z{y-1}(end);
        tau3  = L.t  {y-1};
        Z3    = L.z  {y-1};
        Texp3 = L.T  {y-1};
        T3    = L.T_m{y-1};
        R3    = L.R  {y-1};
        h3    = H    {y-1};
             save('tmp.mat', 'tau', 'Z', 'Texp', 'T', 'R', 'h', 'tau3', 'Z3', 'Texp3', 'T3', 'R3', 'h3');
    else;    save('tmp.mat', 'tau', 'Z', 'Texp', 'T', 'R', 'h');
    end

%     Kg = repmat(1, size(h,1), 1 );
    [Kopt{y}, qv(y)] = fmincon(@compTempirical, repmat(1, size(h,1), 1 ), [], [], [], [], repmat(0, size(h,1), 1), repmat(2.5, size(h,1), 1));

    disp(['year: ' num2str(y) '; Qf: ' num2str( qv(y) ) '; Kopt: ' num2str( Kopt{y}' ) ', Z: ' num2str(reshape(H{y}',1, size(H{y},1)*size(H{y},2)))])

%     out{y} = [ NTA NTAS NZA NRA NRAS itn qv(y) Kopt{y}' reshape(H{y}',1, size(H{y},1)*size(H{y},2))];

%     disp(out{y})

end; clear y

% load('Kopt_Jul4.mat', 'A');
% % for y = [2 4]; A{y}(1    , :                    ) = out{y}; end; clear y; save('Kopt_Jul4.mat', 'A')
% for y = [2 3 4]; A{y}(end+1,1:size(out{y},2)      ) = out{y};
%                  A{y}(end  ,  size(out{y},2)+1:end) = NaN;    end; clear y; save('Kopt_Jul4.mat', 'A')

end; clear itn
end; clear unv
end; clear par


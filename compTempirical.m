function [ Q ] = compTempirical( Kg )
%compTempirical: quantify agreement between measured and simulated temperature
%evolution provided an array of thermal conductivities Kg (in [J sec^-1 m^-1 K^-1])

    load('tmp.mat');
    
        n = size(h,1);
        Ktmp = nan( size(Z) );              % descrete layers
        for iz = 1:max(1, n-1)     
            [~, ind1] = min( abs( Z - h(iz,1) ) );
            [~, ind2] = min( abs( Z - h(iz,2) ) );  
            Ktmp(ind1:ind2) = Kg(iz);
        end
        if n>1; Ktmp(ind2:end ) = Kg(iz+1); end; clear iz ind*
        [ Tm ] = condfill( tau, Z, T, R, Ktmp );
    
    if exist('h3', 'var')
        n = size(h3,1);
        Ktmp3 = nan( size(Z3) );
        for iz = 1:max(1, n-1)     
            [~, ind1] = min( abs( Z3 - h3(iz,1) ) );
            [~, ind2] = min( abs( Z3 - h3(iz,2) ) );
            Ktmp3(ind1:ind2) = Kg(iz);
        end
        if n>1; Ktmp3(ind2:end ) = Kg(iz+1); end; clear iz ind*
        [ Tm3 ] = condfill( tau3, Z3, T3, R3, Ktmp3 );
    end
    
    %% errors in temperature simulation with respect to the measurements
%     e = abs(Tm - Texp);
    e = (Tm - Texp).^2;
    e(~isnan(T)) = NaN;     % remove IC and BC
    Q = nanmean(e(:));
    
    
   if exist('h3', 'var')
      e3 = abs(Tm3 - Texp3);
      e3(~isnan(T3)) = NaN;
      Q = mean( [nanmean(e3(:)) Q]); 
   end

end
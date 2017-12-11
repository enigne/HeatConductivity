function [ Tm ] = condfill( tau, Z, T, R, K )
%CONDFILL  use FTCS scheme to solve thermal conductivity equation for filling a matrix containing IC and BC
% tau  time vector                                                , size N
% Z    depth vector                                               , size M
% T    temperature matrix to be filled in should contain IC and BC, size M*N
% R    density matrix                                             , size M*N
% K    thermal conductivity vector                                , size M

% tau = D{y}.t;
% Z   = D{y}.z;
% T   = D{y}.Tmod{l};
% R   = D{y}.R;
% K   = repmat(K, size(D{y}.t));

sta = min(find(nansum(abs(T))>0))+1;

for t = sta:size( tau, 2)
    
%     disp(['t = ' num2str(t) ' of ' num2str(size( tau, 2))] )

    ind = find( isfinite(T(:,t)) );
    [~, dind] = max(diff(ind));
    up = ind(dind)+1;
    dn = ind(dind+1)-1; clear ind dind
    
    dz  =    Z ( 2 ) -  Z ( 1 );
    dt  = ( tau( t ) - tau(t-1) ) * 24 * 3600;
    Ktmp = K( up-1:dn+1      );
    Rtmp = R( up-1:dn+1, t-1 );
    Ttmp = T( up-1:dn+1, t-1 );
    Ctmp = 152.5 + 7.122.*(273.15 + Ttmp);
    CFL =  max(    dt   .*  Ktmp./ ...
                   dz^2 ./ Rtmp ./ Ctmp );
%     CFL
%     pause
    
    n_dt = ceil(CFL/0.45); clear Ktmp Rtmp Ttmp Ctmp CFL
    
    if  n_dt == 1
        for z = up:dn
             T11 = T( z-1, t-1 );
             T21 = T( z  , t-1 );
             T31 = T( z+1, t-1 );
             R21 = R( z  , t-1 );
             R22 = R( z  , t   );
             K1  = K( z-1      );
             K2  = K( z        );
             K3  = K( z+1      );
             C = 152.5+7.122*(273.15+T21);    % formulation for specific heat capacity of ice from Cuffey and Paterson, 2010
             
             T(z,t) = T21 + dt / ( R21 + R22) / C / dz^2 * ( (K3+K2)*(T31-T21) - (K2+K1)*(T21-T11));
             
         end; clear z T11 T21 T31 R21 R22 K1 K2 K3 C
    else
        tmp = NaN( size(T(:,t),1), n_dt+1  );
        tmp( up-1 , : ) = interp1([tau(t-1) tau(t)], [T(up-1, t-1) T(up-1, t)], [tau(t-1):dt/24/3600/n_dt:tau(t)]);
        tmp( dn+1 , : ) = interp1([tau(t-1) tau(t)], [T(dn+1, t-1) T(dn+1, t)], [tau(t-1):dt/24/3600/n_dt:tau(t)]);
        tmp( up:dn, 1 ) = T( up : dn, t-1);
        for  tn = 2:size(tmp,2)
        for  z  = up:dn
             T11 = tmp( z-1, tn-1 );
             T21 = tmp( z  , tn-1 );
             T31 = tmp( z+1, tn-1 );
             R21 =   R( z  ,  t-1 );
             R22 =   R( z  ,  t   );
             K1  =   K( z-1       );
             K2  =   K( z         );
             K3  =   K( z+1       );
             C = 152.5+7.122*(273.15+T21);    % formulation for specific heat capacity of ice from Cuffey and Paterson, 2010
             
             tmp(z,tn) = T21 + dt/n_dt / ( R21 + R22) / C / dz^2 * ( (K3+K2)*(T31-T21) - (K2+K1)*(T21-T11));
        end; clear z T11 T21 T31 R21 R22 K1 K2 K3 C
        end;
        T(up:dn, t) = tmp( up:dn, end);
        
        clear tmp tn
    end; clear dt CFL n_dt
end; clear t

Tm = T;
end
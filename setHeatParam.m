function heatParam = setHeatParam(dt, Nt, dz, Nz, rho, C, T0, Tbc, zK)
    heatParam.dt = dt;
    heatParam.Nt = Nt;
    heatParam.dz = dz;
    heatParam.Nz = Nz;
    heatParam.rho = rho;
    heatParam.C = C;
    heatParam.T0 = T0;
    heatParam.Tbc = Tbc;   
    heatParam.zK = zK;   
end
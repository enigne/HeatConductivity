function heatParam = setHeatParam(dt, Nt, dz, Nz, rho, C, T0, TbcUp, TbcDown, zK)
    heatParam.dt = dt;
    heatParam.Nt = Nt;
    heatParam.dz = dz;
    heatParam.Nz = Nz;
    heatParam.rho = rho;
    heatParam.C = C;
    heatParam.T0 = T0;
    heatParam.TbcUp = TbcUp;   
    heatParam.TbcDown = TbcDown;   
    heatParam.zK = zK;   
end
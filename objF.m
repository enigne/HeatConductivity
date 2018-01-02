function F = objF(K, HeatSolver, heatParam, z, t, fdata)
    T_sol = HeatSolver(t, z, K, heatParam);
 
    err = T_sol - fdata;   
    mask = (fdata<-2);
    
    compareErr = err(mask);
    F = compareErr(:);
end
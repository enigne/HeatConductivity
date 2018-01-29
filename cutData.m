function [T_data, z_data, indCutZ] = cutData(T_data, z_data, z)
    zminInd = find(z_data >= z(1), 1, 'first');
    zmaxInd = find(z_data <= z(2), 1, 'last');
    
    indCutZ = zminInd:zmaxInd;
    z_data = z_data(indCutZ);
    T_data = T_data(indCutZ, :);
end
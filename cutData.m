function [T_data, z_data] = cutData(T_data, z_data, z)
    zminInd = find(z_data >= z(1), 1, 'first');
    zmaxInd = find(z_data <= z(2), 1, 'last');
    z_data = z_data(zminInd:zmaxInd);
    T_data = T_data(zminInd:zmaxInd, :);
end
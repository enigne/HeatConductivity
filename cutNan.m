function [T_data, t_data, noNan] = cutNan(T_data, t_data)
    flagNan = sum(isnan(T_data), 1);
    noNan = (flagNan == 0);
    T_data = T_data(:, noNan);
    t_data = t_data(noNan);
end
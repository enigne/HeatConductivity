function [T_data, cord_data, indCut] = cutData(T_data, cord_data, cutRange, timeFlag)
    if nargin < 4
        timeFlag = 0;
    end
    
    minInd = find(cord_data >= cutRange(1), 1, 'first');
    maxInd = find(cord_data <= cutRange(2), 1, 'last');
    
    indCut = minInd:maxInd;
    cord_data = cord_data(indCut);
    if (~timeFlag) 
        % cut in space 
        T_data = T_data(indCut, :);
    else
        % cut in time 
        T_data = T_data(:, indCut);
    end
    
end
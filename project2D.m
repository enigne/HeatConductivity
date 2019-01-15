% Project data with coordinates (xdata,ydata) to a 
% new grid (X,Y) with linear interpolation on Y coordinates only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Pdata = project2D(data, xdata, ydata, X, Y)
%     nX = size(data, 2);
%     nY = length(Y);
%     
%     Pdata = zeros(nY, nX);
% 
%     for i = 1: nX
%         Pdata(:,i) = interp1(ydata, data(:,i), Y, 'spline', 'extrap');
%     end

    Pdata = interp2(xdata, ydata, data, X, Y);

end
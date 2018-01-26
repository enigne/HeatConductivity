% Project data with coordinates (xdata,ydata) to a 
% new grid (X,Y) with linear interpolation on Y coordinates only

function Pdata = project2DY(data, ydata, Y)
    nX = size(data, 2);
    nY = length(Y);
    
    Pdata = zeros(nY, nX);

    for i = 1: nX
        Pdata(:,i) = interp1(ydata, data(:,i), Y, 'spline', 'extrap');
    end
end
% Project data with coordinates (xdata,ydata) to a 
% new grid (X,Y) with linear interpolation on Y coordinates only

function Pdata = projectY2D(data, xdata, ydata, X, Y)
    nX = length(X);
    nY = length(Y);
    
    Pdata = zeros(nY, nX);

    for i = 1: nX
        Pdata(:,i) = interp1(ydata, data(:,i), Y, 'spline', 'extrap');
    end

end
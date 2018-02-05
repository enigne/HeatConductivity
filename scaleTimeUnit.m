function t = scaleTimeUnit(t_old, oldUnit, newUnit)
    % Check the input variables
    if nargin < 3
        newUnit = 'days';
        if nargin < 2
            oldUnit = 'seconds';
        end
    end

    % convert anyunit to seconds 
    switch oldUnit
        case 'seconds'
            t_conv = t_old;
        case 'days'
            t_conv = t_old * 24*3600;
        otherwise
            t_conv = datenum(t_old);
    end

    % convert seconds to newUnit
    switch newUnit
        case 'days'
            t = (t_conv-t_conv(1))/24/3600 +1;
        case 'seconds'
            t = t_conv;
        otherwise
            t = datevec(t_conv);
    end
    
end
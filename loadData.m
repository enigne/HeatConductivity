% Load the mearsurements from given data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The input variables:
%   'data'      - the whole data set as a structure, with:
%               	'data.t'        - the time vector of the measurements;
%                   'data.z'        - the depth of the measurements;
%                   'data.T'        - the temperature measurements;
%                   'data.T_corr'   - the temperature after correction;
%                   'data.z_corr'   - the depth after correction;
%                   'data.T_i'      - the interpolated temperature;
%                   'data.z_i'      - the depth for the interpolation;
%                   'data.z_a'      - the averaged depth;
%                   'data.T_a'      - the averaged interpolated temperature;
%                   'data.T_sd'     - the standard deviation of the
%                                     temperature over the 9 holes;
%   'dataIndex'	- index of the holes to be used, 0 means to use the average
%                 date;
%   'timePeriod'- part of t_data are taken into account. [0,1] indicates
%                 the full data piece;
%   'timeScaling'- The scaling factor of time unit, currently using (days).
% The return values:
%   't_data'	- time series;
%   'z_data'	- z-coordinates;
%   'T_data'	- temperature measurements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-02-01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t_data, z_data, T_data, indt] = loadData(data, dataIndex, timePeriod, timeScaling)
    % optional input
    if nargin < 4
        timeScaling = 24*3600;
        if nargin < 3
            timePeriod = [0, 1];
            if nargin < 2
                dataIndex = 0;
            end
        end
    end

    % Load Measurements
    if (dataIndex > 0) && (dataIndex < length(data.T_i))
        % Take individual measurements
        % Position vector (column)
        z_data = data.z_i{dataIndex};
        % Temperature
        T_data = data.T_i{dataIndex};
    else
        % Take the averaged value
        % Position vector (column)
        z_data = data.z_a;
        % Temperature
        T_data = data.T_a;        
    end
    
    % Time vector (row)
    t_data = data.t'  .* timeScaling;
    indt = 1:length(t_data);
        
    % Cut data in time series
    if (timePeriod(1) > 0) || (timePeriod(2) < 1)
        cutRange = interp1([0, 1],[min(t_data), max(t_data)], timePeriod);
        [T_data, t_data, indt] = cutData(T_data, t_data, cutRange, 1);
    end
end
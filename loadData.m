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
%   'timeScaling'- The scaling factor of time unit, currently using (days).
% The return values:
%   't_data'	- time series;
%   'z_data'	- z-coordinates;
%   'T_data'	- temperature measurements.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Cheng Gong
% Date: 2018-01-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t_data, z_data, T_data] = loadData(data, dataIndex, timeScaling)
    % optional input
    if nargin < 3
        timeScaling = 1;
        if nargin < 2
            dataIndex = 0;
        end
    end

    % Load Measurements
    if (dataIndex > 0) && (dataIndex < length(data.T_i))
        % Take individual measurements
        % Time vector (row)
        t_data = data.t'  .* timeScaling;
        % Position vector (column)
        z_data = data.z_i{dataIndex};
        % Temperature
        T_data = data.T_i{dataIndex};
    else
        % Take the averaged value
        % Time vector (row)
        t_data = data.t'  .* timeScaling;
        % Position vector (column)
        z_data = data.z_a;
        % Temperature
        T_data = data.T_a;        
    end
end
% remove_large_jumps
% This function sets values to NaN where the difference between adjacent points
% in the input `data` exceeds the specified `threshold`.
%
% INPUTS:
%   data      - 1D array of data points (e.g., position data over time).
%   threshold - Threshold for detecting large jumps (e.g., 180 degrees).
%
% OUTPUT:
%   data_processed - Modified data where large jumps are replaced with NaN.
%
function data_processed = remove_large_jumps(data, threshold)
% Calculate the difference between adjacent points
diff_data = abs(diff(data));

% Initialize the output data
data_processed = data;

% Set points to NaN where the difference exceeds the threshold
data_processed([false, diff_data > threshold]) = NaN;
end

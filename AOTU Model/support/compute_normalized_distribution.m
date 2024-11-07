% compute_normalized_distribution
%
% This function computes the normalized histogram of combined input data 
% across multiple trial runs and time points for both right and left inputs,
% using all possible input values from runSettings.DNa02input. The histogram 
% is generated across bins with a specified size and normalized to a maximum of 1.
%
% INPUTS:
%   input_history     - 3D array representing the input history for multiple trial runs.
%                       Dimensions:
%                       - x: Number of trial runs
%                       - y: Time points
%                       - z: Input history to the right (1) vs. left (2)
%   runSettings       - Struct containing the parameters for the run, including:
%                       - DNa02input: Array of all possible input values.
%                       - binSize (optional): Size of each bin for the histogram (default is 0.1).
%
% OUTPUTS:
%   normalized_histogram - Array containing the normalized histogram of the 
%                          combined distribution of inputs, with a maximum value of 1.
%   bin_centers          - Array containing the center values of each bin for plotting.
%
% CREATED: 11/05/2024 - MC

function [normalized_histogram,bin_centers] = compute_normalized_distribution(input_history, runSettings)

% Set bin size
binSize = 0.25;

% Define bin edges based on DNa02input and binSize
minInput = min(runSettings.DNa02input);
maxInput = max(runSettings.DNa02input);
binEdges = minInput:binSize:maxInput;
% Calculate bin centers
bin_centers = binEdges(1:end-1) + binSize / 2;

% Compute histogram counts for the combined data
histogram_counts = histcounts(input_history, binEdges);

% Normalize histogram counts to a maximum of 1
normalized_histogram = histogram_counts / max(histogram_counts);

end
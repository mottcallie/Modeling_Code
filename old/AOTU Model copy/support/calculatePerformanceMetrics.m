% calculatePerformanceMetrics
% Function to compute various performance metrics, including circular standard
% deviation, variance, probability near setpoint, mean Integral of Squared Error (ISE),
% mean Integral of Absolute Error (IAE), and direction changes.
%
% INPUTS:
%   visobj_history - Trajectory of the visual object over time.
%   rotvel_history - Rotational velocity over time.
%   timebase       - Timebase for the simulation.
%   plotSettings   - Structure containing plot settings.
%
% OUTPUTS:
%   metric         - Structure containing performance metrics such as:
%                   - std: Circular standard deviation.
%                   - var: Circular variance.
%                   - prob: Probability of being within +/- 5 degrees of setpoint.
%                   - ISE: Mean Integral of Squared Error.
%                   - IAE: Mean Integral of Absolute Error.
%
function metric = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings)
    % Calculate circular standard deviation
    metric.std = median(circ_std(deg2rad(visobj_history(:,:)), [], [], 2));

    % Calculate circular variance
    metric.var = median(1 - circ_var(deg2rad(visobj_history(:,:)), [], [], 2));

    % Calculate the probability of being within +/- 5 degrees of setpoint
    [h, e] = histcounts(visobj_history(:,:), 'BinWidth', plotSettings.hdBinSize, 'Normalization', 'probability');
    e = e(1:end-1) + plotSettings.hdBinSize / 2; % Center bins
    metric.prob = sum(h(abs(e) <= 5));

    % Calculate mean ISE (Integral of Squared Error) across the trajectory
    squared_error = (visobj_history(:,:) - 0).^2; % Setpoint is 0
    metric.ISE = mean(trapz(squared_error, 2)); % Mean ISE

    % Calculate mean IAE (Integral of Absolute Error) across the trajectory
    absolute_error = abs(visobj_history(:,:) - 0); % Setpoint is 0
    metric.IAE = mean(trapz(absolute_error, 2)); % Mean IAE
end

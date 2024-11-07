% analyzeObjectAndVelocityCrossings
% Function to analyze zero-crossings in object position and rotational velocity,
% calculate the time delay between when the object crosses zero and when the animal
% changes directions (indicated by rotational velocity zero-crossing), and bin the 
% crossing delays by both rotational velocity and instantaneous object velocity.
%
% INPUTS:
%   visobj_history       - Matrix of object position history (numRuns x nTime).
%   rotvel_history       - Matrix of rotational velocity history (numRuns x nTime).
%   timebase             - Time axis for the simulation, in seconds.
%   nTest                - Number of trials to analyze, limited by available data.
%
% OUTPUTS:
%   all_rotvel_crosstimes - Times for rotational velocity to cross zero after each object zero-crossing.
%   vel_bins              - Bin edges for velocity at crossing
%   rotvel_binned_avgs    - Binned averages of crossing times within each rotational velocity bin.
%   visobj_binned_avgs    - Binned averages of crossing times within each object velocity bin.
%
% The function calculates the instantaneous velocity of the object position (`visobj_history`)
% and bins the time delays of rotational velocity zero-crossings by both the rotational
% velocity and object velocity at the time of each object position zero-crossing.
%
% CREATED: 11/01/2024 - MC
%
function [all_rotvel_crosstimes, vel_bins, rotvel_binned_avgs, visobj_binned_avgs] = ...
    analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest)

%% Ensure nTest does not exceed available runs
maxRuns = min(size(visobj_history, 1), size(rotvel_history, 1));
if nTest > maxRuns
    nTest = maxRuns;
    warning('nTest exceeds the number of available runs. Using maximum dimension: %d', maxRuns);
end

% Set bin edges and analysis parameters
vel_bins = 0:60:440;               % Rotational velocity bins
window_idx = 20;                      % Size of analysis window around each zero-crossing
time_in_ms = (timebase(1:window_idx*2 + 1) - timebase(window_idx+1)) * 1000;  % Time array centered at zero-crossing

%% Initialize output arrays
all_rotvel_crosstimes = [];
all_rotvel_at_cross = [];
all_visobj_crosstimes = [];
all_visobj_at_cross = [];

% Loop through each trial and analyze zero-crossings
for runIdx = 1:nTest
    this_obj_history = visobj_history(runIdx, :);
    this_rotvel_history = rotvel_history(runIdx, :);

    % Calculate instantaneous velocity of visobj
    visobj_velocity = diff(this_obj_history) ./ diff(timebase);

    % Detect zero-crossings in object position
    sign_changes = diff(sign(this_obj_history));
    zero_cross_indices = find(sign_changes ~= 0);

    for i = 1:length(zero_cross_indices)
        idx = zero_cross_indices(i);

        % Ensure the analysis window is within bounds
        if idx - window_idx < 1 || idx + window_idx > length(this_obj_history)
            continue;
        end

        % Extract windows around the zero-crossing in object position
        rotvel_window = this_rotvel_history(idx - window_idx : idx + window_idx);

        % Identify the first rotational velocity zero-crossing after the object zero-crossing
        rotvel_sign_changes = diff(sign(rotvel_window(window_idx+1:end)));
        first_rotvel_cross_idx = find(rotvel_sign_changes ~= 0, 1, 'first');

        % Calculate the delay only if a rotational velocity zero-crossing exists after the object crossing
        if ~isempty(first_rotvel_cross_idx)
            time_delay = (first_rotvel_cross_idx + 1) * mean(diff(timebase)) * 1000; % Convert to milliseconds
            rotvel_at_cross = abs(this_rotvel_history(idx));               % Rotational velocity at crossing
            all_rotvel_crosstimes = [all_rotvel_crosstimes, time_delay];
            all_rotvel_at_cross = [all_rotvel_at_cross, rotvel_at_cross];

            % Fetch the instantaneous velocity of visobj at the crossing
            if idx < length(visobj_velocity)
                visobj_at_cross = abs(visobj_velocity(idx));
                all_visobj_crosstimes = [all_visobj_crosstimes, time_delay];
                all_visobj_at_cross = [all_visobj_at_cross, visobj_at_cross];
            end
        end
    end
end

% Bin cross times by rotational velocity at the object crossing
rotvel_binned_avgs = zeros(1, length(vel_bins)-1);
for binIdx = 1:length(vel_bins)-1
    in_bin = all_rotvel_at_cross >= vel_bins(binIdx) & all_rotvel_at_cross < vel_bins(binIdx+1);
    if any(in_bin)
        rotvel_binned_avgs(binIdx) = mean(all_rotvel_crosstimes(in_bin), 'omitnan');
    else
        rotvel_binned_avgs(binIdx) = NaN;
    end
end

% Bin cross times by visobj velocity at the object crossing
visobj_binned_avgs = zeros(1, length(vel_bins)-1);
for binIdx = 1:length(vel_bins)-1
    in_bin = all_visobj_at_cross >= vel_bins(binIdx) & all_visobj_at_cross < vel_bins(binIdx+1);
    if any(in_bin)
        visobj_binned_avgs(binIdx) = mean(all_visobj_crosstimes(in_bin), 'omitnan');
    else
        visobj_binned_avgs(binIdx) = NaN;
    end
end

%% Optional Plotting of Binned Cross Times
optPlot = 0;
if optPlot
    figure;
    set(gcf, 'Position', [100 100 600 300]);

    % Plot binned cross times by rotational velocity bins
    subplot(1, 2, 1);
    plot(vel_bins(1:end-1), rotvel_binned_avgs);
    xlabel('Rotational Velocity at Crossing (deg/s)');
    ylabel('Avg. Cross Time (ms)');
    title('Binned Cross Times by Rotational Velocity');
    grid on;

    % Plot binned cross times by visobj velocity bins
    subplot(1, 2, 2);
    plot(vel_bins(1:end-1), visobj_binned_avgs);
    xlabel('VisObj Velocity at Crossing (deg/s)');
    ylabel('Avg. Cross Time (ms)');
    title('Binned Cross Times by VisObj Velocity');
    grid on;
end
end
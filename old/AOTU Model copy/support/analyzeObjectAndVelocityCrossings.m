% analyzeObjectAndVelocityCrossings
% Function to analyze zero-crossings in object position and rotational velocity,
% and bin cross times by object position peaks and rotational velocity peaks.
%
% INPUTS:
%   visobj_history    - Object position history (numRuns x nTime).
%   rotvel_history    - Rotational velocity history (numRuns x nTime).
%   timebase          - Time axis for the simulation.
%   nTest             - Number of trials to run against.
%
% OUTPUTS:
%   all_rotvel_crosstimes - Times it takes for rotational velocity to cross zero after object crossing.
%   obj_bins              - Object position bins used for binning cross times.
%   rotvel_bins           - Rotational velocity bins used for binning cross times.
%   obj_binned_avgs       - Binned averages of cross times for object position bins.
%   rotvel_binned_avgs    - Binned averages of cross times for rotational velocity bins.
%
function [all_rotvel_crosstimes, obj_bins, rotvel_bins, obj_binned_avgs, rotvel_binned_avgs] = analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest)
    %% Check nTest against the maximum dimension of the input arrays
    maxRuns = min(size(visobj_history, 1), size(rotvel_history, 1));
    if nTest > maxRuns
        nTest = maxRuns;  % Use the maximum available number of runs
        warning('nTest exceeds the number of available runs. Using maximum dimension: %d', maxRuns);
    end

    %% Initialize output arrays
   
    % Define the object and rotational velocity bins
    obj_bins = -150:20:150;  % Object position bins
    rotvel_bins = -15:0.5:15;  % Rotational velocity bins

    % Preallocate for object and rotvel windows
    window_idx = 15; % Window size in indices
    time_in_ms = (timebase(1:window_idx*2 + 1) - timebase(window_idx+1)) * 1000;  % Convert to ms around zero-crossing

    %% Analyze direction changes
    % Initialize arrays
    all_rotvel_crosstimes = [];
    all_objpeaks = [];
    all_rotvelpeaks = [];

    % Loop through each run
    for runIdx = 1:nTest
        this_obj_history = visobj_history(runIdx, :);
        this_rotvel_history = rotvel_history(runIdx, :);

        % Find zero-crossings in object position (where it crosses zero)
        sign_changes = diff(sign(this_obj_history));
        zero_cross_indices = find(sign_changes ~= 0);

        for i = 1:length(zero_cross_indices)
            if i > length(zero_cross_indices), continue; end  % Avoid out-of-bounds error
            idx = zero_cross_indices(i);

            % Fetch data in the window around the zero-crossing
            if idx - window_idx < 1 || idx + window_idx > length(this_obj_history)
                continue; % Skip if window is out of bounds
            end
            obj_window = this_obj_history(idx - window_idx : idx + window_idx);
            rotvel_window = this_rotvel_history(idx - window_idx : idx + window_idx);

            % Find the first zero-crossing in rotational velocity after the object zero-crossing
            rotvel_sign_changes = diff(sign(rotvel_window(window_idx+1:end)));
            first_rotvel_cross_idx = find(rotvel_sign_changes ~= 0, 1, 'first');

            if ~isempty(first_rotvel_cross_idx)
                % Find the peak closest to the zero-crossing in object position
                [obj_peaks, obj_peak_locs] = max(abs(obj_window(1:window_idx)));  % Find peaks in object position
                if ~isempty(obj_peaks)
                    [~, closest_obj_peak_idx] = min(abs(obj_peak_locs - window_idx));  % Closest peak to zero-crossing
                    obj_peaks = obj_peaks(closest_obj_peak_idx);  % Only store the closest peak
                end

                % Find the peak closest to the zero-crossing in rotational velocity
                [rotvel_peaks, rotvel_peak_locs] = findpeaks(abs(rotvel_window(1:window_idx)));  % Find peaks in rotational velocity
                if ~isempty(rotvel_peaks)
                    [~, closest_rotvel_peak_idx] = min(abs(rotvel_peak_locs - window_idx));  % Closest peak to zero-crossing
                    rotvel_peaks = rotvel_peaks(closest_rotvel_peak_idx);  % Only store the closest peak
                end

                % Only store if both object and rotational velocity peaks are found and all variables are valid
                if ~isempty(obj_peaks) && ~isempty(rotvel_peaks) && ~isempty(first_rotvel_cross_idx)
                    % Time it takes for rotational velocity to cross zero after object crossing
                    all_rotvel_crosstimes = [all_rotvel_crosstimes, time_in_ms(window_idx+1 + first_rotvel_cross_idx)];

                    % Store the peaks (combine right and left crossings by taking the absolute value)
                    all_objpeaks = [all_objpeaks, obj_peaks];
                    all_rotvelpeaks = [all_rotvelpeaks, rotvel_peaks];
                end
            end
        end
    end

    %% Bin cross times by object and rotational velocity peaks
    obj_binned_avgs = zeros(1, length(obj_bins)-1);
    rotvel_binned_avgs = zeros(1, length(rotvel_bins)-1);

    % Bin the object peaks and calculate binned averages for cross times
    for binIdx = 1:length(obj_bins)-1
        in_bin = all_objpeaks >= obj_bins(binIdx) & all_objpeaks < obj_bins(binIdx+1);
        if any(in_bin)
            obj_binned_avgs(binIdx) = mean(all_rotvel_crosstimes(in_bin), 'omitnan');
        else
            obj_binned_avgs(binIdx) = NaN;  % Set to NaN if no values are in the bin
        end
    end

    % Bin the rotational velocity peaks and calculate binned averages for cross times
    for binIdx = 1:length(rotvel_bins)-1
        in_bin = all_rotvelpeaks >= rotvel_bins(binIdx) & all_rotvelpeaks < rotvel_bins(binIdx+1);
        if any(in_bin)
            rotvel_binned_avgs(binIdx) = mean(all_rotvel_crosstimes(in_bin), 'omitnan');
        else
            rotvel_binned_avgs(binIdx) = NaN;  % Set to NaN if no values are in the bin
        end
    end

    %% Skip bins where object max < 10 and rotational velocity max < 5
    %obj_binned_avgs(obj_bins(1:end-1) < 10) = NaN;
    %rotvel_binned_avgs(rotvel_bins(1:end-1) < 5) = NaN;

    %% Optional Plotting
    optPlot=0;
    if optPlot
        figure;
        set(gcf, 'Position', [100 100 1200 600]);

        % Plot binned cross times by object bins
        subplot(1, 2, 1);
        plot(obj_bins(1:end-1), obj_binned_avgs);
        xlabel('Object Max (deg)');
        ylabel('Avg. Cross Time (ms)');
        title('Binned Cross Times by Object Max');
        grid on;

        % Plot binned cross times by rotational velocity bins
        subplot(1, 2, 2);
        plot(rotvel_bins(1:end-1), rotvel_binned_avgs);
        xlabel('Rotvel Max (deg/s)');
        ylabel('Avg. Cross Time (ms)');
        title('Binned Cross Times by Rotvel Max');
        grid on;
    end
end

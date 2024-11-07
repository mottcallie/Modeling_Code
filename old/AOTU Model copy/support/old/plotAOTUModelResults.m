function plotAOTUModelResults(timebase, visobj_history, input_history, predicted_RF, runSettings, plotSettings)
    % plotAOTUModelResults
    % Function to plot the outputs of the AOTU steering model.
    %
    % INPUTS:
    % timebase       - Time axis for the simulation.
    % visobj_history - Trajectory of the visual object over time.
    % input_history  - Summed input from all visual neurons.
    % predicted_RF   - Receptive field tuning curves (RF) for AOTU neurons.
    % runSettings    - Structure containing model parameters.
    % plotSettings   - Structure containing plot settings.
    %
    % OUTPUT:
    % A figure with subplots showing RF tuning curves, example runs, and histograms.
    %
    % 10/05/2024 - MC - Created the function to simplify plotting AOTU model outputs.

    %% Pre-process visual object history to avoid steep jumps at -180 to 180 crossings
    % Detect instances where the visual object crosses between -180 and 180 degrees
    for i = 1:size(visobj_history, 1) % Loop through each simulation run
        for t = 2:size(visobj_history, 2) % Loop through each time step
            % Check if the object crosses from -180 to 180 or vice versa
            if abs(visobj_history(i, t) - visobj_history(i, t-1)) > 180
                visobj_history(i, t) = NaN; % Set to NaN to prevent steep jumps in the plot
            end
        end
    end

    %% Plot 1: Visual receptive field (RF) tuning curves
    nexttile; hold on
    % Loop through and plot each RF curve for AOTU neurons (Right and Left)
    for i = 1:3
        % Plot right tuning curves
        plot(runSettings.visObjPosition, predicted_RF{:, i}, "Color", plotSettings.rfColors{i});
        % Plot left tuning curves (flipped)
        plot(runSettings.visObjPosition, flip(predicted_RF{:, i}), "Color", plotSettings.rfColors{i});
    end
    box off; % Remove box border
    xlim([-180 180]); % Set x-axis limits
    xticks(plotSettings.rfTicks); % Set x-ticks at defined angles
    ylabel('Firing Rate (Hz)'); % Label y-axis
    xlabel('Target Position (deg)'); % Label x-axis

    %% Plot 2: Example simulation runs (visual object trajectory over time)
    nexttile([1 3]); % Create a larger tile for example runs
    plot(timebase(1, 2:end), visobj_history(1:plotSettings.nEx, 2:end), '-'); % Plot trajectories of the first nEx runs
    axis tight; ylim([-180 180]); % Set y-axis limits
    yticks([-180 0 180]); % Set y-ticks at -180, 0, and 180 degrees
    xticks([]); % Remove x-ticks
    yline(0); % Add a horizontal line at 0 degrees
    box off; % Remove box border
    xlabel('Time (s)'); % Label x-axis
    ylabel('Target Position (deg)'); % Label y-axis

    %% Plot 3: Polar histogram of target positions
    nexttile; % Create the next tile for polar histogram
    polarhistogram(deg2rad(visobj_history(:, 3:end)), plotSettings.nPolarBins, ... 
        'FaceColor', '#77AC30', 'FaceAlpha', plotSettings.opacity, 'Normalization', 'probability');

    %% Plot 4: Standard histogram of target positions
    nexttile; % Create the next tile for standard histogram
    histogram(visobj_history(:, 3:end), 'BinWidth', plotSettings.hdBinSize, ...
        'FaceColor', '#77AC30', 'FaceAlpha', plotSettings.opacity, 'Normalization', 'probability');
    xlim([-180 180]); % Set x-axis limits
    xticks(plotSettings.rfTicks); % Set x-ticks to match RF tuning plot
    ylim(plotSettings.hHistRange); % Set y-axis limits
    xline(0); % Add a vertical line at 0 degrees
    xlabel('Target Position (deg)'); % Label x-axis

end

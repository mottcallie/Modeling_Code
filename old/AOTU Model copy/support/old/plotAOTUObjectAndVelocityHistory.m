function plotAOTUObjectAndVelocityHistory(timebase, visobj_history, rotvel_history, runSettings)
    % plotAOTUObjectAndVelocityHistory
    % Function to plot the visual object history and rotational velocity
    % for the first 5 trials, each on the same axis.
    %
    % INPUTS:
    %   timebase        - Time axis for the simulation.
    %   visobj_history  - Trajectory of the visual object over time (numRuns x nTime).
    %   rotvel_history  - Rotational velocity output over time (numRuns x nTime).
    %   runSettings     - Structure containing model parameters (for axis settings, etc.).
    %
    % OUTPUT:
    %   A figure with multiple rows:
    %     - Each row corresponds to a trial (first 5 trials), where both
    %       the visual object history and rotational velocity are plotted on the same axis.
    %

    % Number of trials to plot
    nPlot = min(3, size(visobj_history, 1));  % Plot the first 5 trials or fewer if less data

    % Set up the figure and layout
    tiledlayout(nPlot, 1, 'TileSpacing', 'compact');  % Create a tiled layout with nPlot rows

    % Loop through the first nPlot trials
    for i = 1:nPlot
        % Create a new tile for each trial
        nexttile;
        hold on;

        % Plot visual object history
        yyaxis left;  % Use the left y-axis for object position
        plot(timebase, visobj_history(i, :), '-b', 'DisplayName', 'Object Position');
        ylim([-180 180]);  % Object position typically ranges between -180 and 180 degrees
        ylabel('Object Position (deg)');
        yline(0, '--', 'Setpoint');  % Add a horizontal line at 0 degrees (setpoint)

        % Plot rotational velocity history
        yyaxis right;  % Use the right y-axis for rotational velocity
        plot(timebase, rotvel_history(i, :), '-r', 'DisplayName', 'Rotational Velocity');
        ylabel('Rotational Velocity (deg/s)');

        % Set axis limits and labels
        xlim([min(timebase) max(timebase)]);
        xlabel('Time (s)');
        title(['Trial ' num2str(i)]);

        % Optionally, adjust y-limits for rotational velocity based on runSettings
        if isfield(runSettings, 'rotVelLimits')
            yyaxis right;  % Apply limits to the rotational velocity axis
            ylim(runSettings.rotVelLimits);
        end
        box off;
    end
end

% generateSummaryCrossing.m
%
% This function generates summary plots comparing the average cross times for rotational
% velocity bins and object velocity bins across different steering gain (k) values. 
% Each gain value is displayed in a separate tile, and differences between the two conditions 
% are plotted for both metrics, with gridlines and axis limits specified.
%
% INPUTS:
%   nK                  - Number of steering gain values (k) in the experiment.
%   kValues             - Vector of k values used in the experiment.
%   rotvel_binned_avg1  - Matrix of rotational velocity binned averages for Condition 1.
%   rotvel_binned_avg2  - Matrix of rotational velocity binned averages for Condition 2.
%   objvel_binned_avg1  - Matrix of object velocity binned averages for Condition 1.
%   ovjvel_binned_avg2  - Matrix of object velocity binned averages for Condition 2.
%   vel_bins            - Vector defining bin edges for rotational velocity and object velocity.
%   comparisonLabel     - Cell array with labels for the two conditions being compared.
%   folder              - Folder directories
%
% OUTPUTS:
%   Generates a figure with two rows of plots showing cross-time differences for each condition 
%   across k values, with the first row for rotational velocity and the second row for object velocity.
%
% CREATED: 10/30/2024 - MC
% UPDATED: 11/02/2024 - MC added visobj bins
%
function generateSummaryCrossing(nK, kValues, rotvel_binned_avg1, rotvel_binned_avg2, ...
    objvel_binned_avg_1, objvel_binned_avg_2, vel_bins, comparisonLabel, folder)

    % Generate summary plots
    figure;
    set(gcf, 'Position', [100 100 1500 800]);  % Adjust figure size to accommodate two rows
    tiledlayout(2, nK, 'TileSpacing', 'compact');  % Two rows: one for each velocity type

    % Row 1: Binned Cross Times vs Rotational Velocity
    for kIdx = 1:nK
        nexttile; hold on;
        % Plot cross times for both conditions separately
        plot(vel_bins(1:end-1), rotvel_binned_avg1(kIdx, :), '-o', 'DisplayName', comparisonLabel{1});
        plot(vel_bins(1:end-1), rotvel_binned_avg2(kIdx, :), '-o', 'DisplayName', comparisonLabel{2});

        xlabel('Turn at Cross (deg/s)');
        ylabel('Avg. Cross Time');
        title(['k=' num2str(kValues(kIdx))]);

        if kIdx == 1
            legend('show', 'Location', 'northeast');
        end

        set(gca, 'XScale', 'log')  % Set x-axis to log scale
        set(gca, 'XTick', vel_bins)
        ylim([100 250]);  % Consistent y-axis limits
        xlim([0 400]);    % Consistent x-axis limits
        grid on;
    end

    % Row 2: Binned Cross Times vs Object Velocity
    for kIdx = 1:nK
        nexttile; hold on;
        % Plot cross times for both conditions separately for object velocity
        plot(vel_bins(1:end-1), objvel_binned_avg_1(kIdx, :), '-o', 'DisplayName', comparisonLabel{1});
        plot(vel_bins(1:end-1), objvel_binned_avg_2(kIdx, :), '-o', 'DisplayName', comparisonLabel{2});

        xlabel('Object Velocity at Cross (deg/s)');
        ylabel('Avg. Cross Time');
        title(['k=' num2str(kValues(kIdx))]);

        if kIdx == 1
            legend('show', 'Location', 'northeast');
        end

        set(gca, 'XScale', 'log')  % Set x-axis to log scale
        set(gca, 'XTick', vel_bins)
        ylim([100 250]);  % Same y-axis limits as Row 1
        xlim([0 400]);    % Same x-axis limits as Row 1
        grid on;
    end

    % Overall title for the figure
    sgtitle(['Time to Direction Change Binned by Rotational and Object Velocity at Crossing (' comparisonLabel{1} ' vs ' comparisonLabel{2} ')']);

    % Save figure as PNG and SVG in specified folders
    saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_DirChange.png']));
    set(gcf, 'renderer', 'Painters');  % Use 'Painters' renderer to maintain vector quality
    saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_DirChange.svg']));
end

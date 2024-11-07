% generateSummaryPlots2
%
% This function generates summary plots for performance metrics and binned averages across multiple scaling factors 
% and steering gain (`k`) values after running the AOTU steering model. It compares various metrics (probability, 
% variance, ISE, IAE) and visualizes differences in binned averages between different scaling factors.
%
% INPUTS:
%   kValues          - Array of steering gain (`k`) values used in the simulation.
%   metrics_prob     - Matrix of probability metrics (rows: `k` values, columns: scaling factors).
%   metrics_var      - Matrix of circular variance metrics (rows: `k` values, columns: scaling factors).
%   metrics_ISE      - Matrix of Integral of Squared Error (ISE) metrics (rows: `k` values, columns: scaling factors).
%   metrics_IAE      - Matrix of Integral of Absolute Error (IAE) metrics (rows: `k` values, columns: scaling factors).
%   obj_binned_avgs  - 3D array of binned averages for object position (dimensions: `k`, bins, scaling factors).
%   rotvel_binned_avgs - 3D array of binned averages for rotational velocity (dimensions: `k`, bins, scaling factors).
%   evt              - 3D array of binned angular velocity vs object position (dimensions: `k`, bins, scaling factors).
%   obj_bins         - Array of object position bins used for binned averages.
%   rotvel_bins      - Array of rotational velocity bins used for binned averages.
%   posBins          - Array of position bins used for angular velocity vs object position.
%   nK               - Number of steering gain (`k`) values.
%   comparisonLabel  - Cell array of labels for the different scaling factors (e.g., `{'ScaleFactor=0', 'ScaleFactor=1', ...}`).
%   folder           - Struct containing folder paths where the plots will be saved.
%   plotSuffix       - String that will be appended to plot titles and filenames (e.g., 'Feedforward', 'OpenLoop').
%
% OUTPUTS:
%   The function generates and saves the following plots:
%   - Summary of performance metrics (probability, variance, ISE, IAE) across scaling factors and steering gains.
%   - Binned average cross times compared to the baseline (e.g., ScaleFactor=0) for object position and rotational velocity.
%   - Binned angular velocity vs object position for all scaling factors.
%
% DESCRIPTION:
%   1. The function first generates summary plots for probability, circular variance, ISE, and IAE metrics across all `k` values and scaling factors.
%   2. Then, it generates two plots comparing the binned cross times for object position and rotational velocity against the baseline (e.g., ScaleFactor=0) for all `k` values.
%   3. Finally, it generates plots showing the relationship between binned angular velocity and object position for each `k` value, with lines representing each scaling factor.
%   4. The suffix specified in `plotSuffix` is used to customize the plot titles and filenames (e.g., 'Feedforward', 'OpenLoop').
%
% EXAMPLE USAGE:
%   generateSummaryPlots2(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, ...
%       obj_binned_avgs, rotvel_binned_avgs, evt, obj_bins, rotvel_bins, posBins, nK, comparisonLabel, folder, 'Feedforward');
%
% DEPENDENCIES:
%   This function requires pre-calculated metrics and binned averages from the simulation output.
%
% NOTES:
%   - The binned averages for object position and rotational velocity are compared to the baseline (e.g., ScaleFactor=0).
%   - Ensure that all input data arrays (e.g., metrics, binned averages) are correctly populated before calling this function.
%
function generateSummaryPlots2(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, ...
    obj_binned_avgs, rotvel_binned_avgs, evt, obj_bins, rotvel_bins, posBins, nK, comparisonLabel, folder, plotSuffix)

nFactors = size(metrics_prob, 2);  % Dynamically detect the number of manipulation factors (second dimension)
baselineFactor = 1;  % Assuming the first factor (index 1) is the baseline
condColors = lines(nFactors);  % Generate colors for each feedforward level

%% Generate summary plots
% Skip if open loop
if ~contains(plotSuffix,'OL')
    figure;
    set(gcf, 'Position', [100 100 400 800]);  % Set figure size
    tiledlayout(4, 1, 'TileSpacing', 'compact');  % One row for each metric

    % Plot 1: Probability near setpoint across k values for all levels
    nexttile;
    hold on;
    for factorIdx = 1:nFactors
        plot(kValues, metrics_prob(:, factorIdx), '-o', 'DisplayName', comparisonLabel{factorIdx}, 'Color', condColors(factorIdx, :));
    end
    xlabel('Steering Gain (k)');
    ylabel('Probability Near Setpoint');
    legend('show','FontSize',6);
    title('Probability Near Setpoint');
    axis padded;
    ylim([0 0.5]);
    grid on;

    % Plot 2: Circular variance across k values for all levels
    nexttile;
    hold on;
    for factorIdx = 1:nFactors
        plot(kValues, metrics_var(:, factorIdx), '-o', 'DisplayName', comparisonLabel{factorIdx}, 'Color', condColors(factorIdx, :));
    end
    xlabel('Steering Gain (k)');
    ylabel('1- Circular Variance');
    title('Circular Variance');
    axis padded;
    ylim([0 1]);
    grid on;

    % Plot 3: Integral of Squared Error (ISE) across k values for all levels
    nexttile;
    hold on;
    for factorIdx = 1:nFactors
        plot(kValues, metrics_ISE(:, factorIdx), '-o', 'DisplayName', comparisonLabel{factorIdx}, 'Color', condColors(factorIdx, :));
    end
    xlabel('Steering Gain (k)');
    ylabel('ISE');
    title('Squared Error');
    axis padded;
    grid on;

    % Plot 4: Integral of Absolute Error (IAE) across k values for all levels
    nexttile;
    hold on;
    for factorIdx = 1:nFactors
        plot(kValues, metrics_IAE(:, factorIdx), '-o', 'DisplayName', comparisonLabel{factorIdx}, 'Color', condColors(factorIdx, :));
    end
    xlabel('Steering Gain (k)');
    ylabel('IAE');
    title('Absolute Error');
    axis padded;
    grid on;

    % Save the summary plot using join for filename construction
    sgtitle('Summary of Model Performance Across Metrics');
    saveas(gcf, fullfile(folder.final, join(['SummaryMetrics_', plotSuffix, '.png'], '')));
    set(gcf,'renderer','Painters');
    saveas(gcf, fullfile(folder.vectors, join(['SummaryMetrics_', plotSuffix, '.svg'], '')));
end

%% Generate two summary plots for binned averages vs. object position and rotational velocity
figure;
set(gcf, 'Position', [100 100 1500 600]);  % Set figure size
tiledlayout(2, nK, 'TileSpacing', 'compact');  % One row for each k value

% Compare all runs to baseline for object position binned averages
for kIdx = 1:nK
    nexttile;
    hold on;
    for factorIdx = 2:nFactors  % Start from 2, since we compare to the baseline (index 1)
        this_diff = obj_binned_avgs(kIdx, :, baselineFactor) - obj_binned_avgs(kIdx, :, factorIdx);
        plot(obj_bins(1:end-1), this_diff, '-o', 'DisplayName', [comparisonLabel{baselineFactor} ' - ' comparisonLabel{factorIdx}], ...
            'Color', condColors(factorIdx, :));
    end
    xlabel('Object Max (deg)');
    ylabel('Avg. Cross Time Diff');
    title(['k=' num2str(kValues(kIdx))]);
    if kIdx == 1
        legend('show','Location','southeast');
    end
    ylim([-100 100]);
    xlim([0 150]);
    grid on;
end

% Compare all runs to baseline for rotational velocity binned averages
for kIdx = 1:nK
    nexttile;
    hold on;
    for factorIdx = 2:nFactors
        this_diff = rotvel_binned_avgs(kIdx, :, baselineFactor) - rotvel_binned_avgs(kIdx, :, factorIdx);
        plot(rotvel_bins(1:end-1), this_diff, '-o', 'DisplayName', [comparisonLabel{baselineFactor} ' - ' comparisonLabel{factorIdx}], ...
            'Color', condColors(factorIdx, :));
    end
    xlabel('Rotational Velocity Max (deg/s)');
    ylabel('Avg. Cross Time Diff');
    title(['k=' num2str(kValues(kIdx))]);
    if kIdx == 1
        legend('show','Location','northeast');
    end
    ylim([-100 100]);
    xlim([0 10]);
    grid on;
end
sgtitle(['Binned Cross Times vs ' comparisonLabel{baselineFactor}]);
saveas(gcf, fullfile(folder.final, join(['Binned_DirChange_', plotSuffix, '.png'], '')));
set(gcf,'renderer','Painters');
saveas(gcf, fullfile(folder.vectors, join(['Binned_DirChange_', plotSuffix, '.svg'], '')));

%% Generate Summary Plot for Angular Velocity vs Object Position
figure;
set(gcf, 'Position', [100 100 1500 500]);  % Set figure size
tiledlayout(1, nK, 'TileSpacing', 'compact');  % One row for each k value

% Loop over k values for plotting
for kIdx = 1:nK
    nexttile;
    hold on;
    
    % Plot the relationship for all levels (lines only)
    for factorIdx = 1:nFactors
        plot(posBins, evt(kIdx, :, factorIdx), '-', 'DisplayName', comparisonLabel{factorIdx}, 'Color', condColors(factorIdx, :));
    end

    % Customize the plot
    xlabel('Position (deg)');
    ylabel('Avg. Angular Velocity (deg/s)');
    title(['k=' num2str(kValues(kIdx))]);
    if kIdx == 1
        legend('show', 'Location', 'northwest');
    end
    grid on;
    ylim([-10 10]);
    xlim([-50 50]);
end

% Save the summary plot
sgtitle(['Binned Angular Velocity vs Object Position for Different ' plotSuffix ' Levels Across k Values']);
saveas(gcf, fullfile(folder.final, join(['Binned_EVT_', plotSuffix, '.png'], '')));
set(gcf,'renderer','Painters');
saveas(gcf, fullfile(folder.vectors, join(['Binned_EVT_', plotSuffix, '.svg'], '')));

end

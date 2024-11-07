% generateSummaryMetrics.m
%
% This function generates summary plots for various metrics across steering gain values (k values),
% comparing two conditions side by side. Each metric is displayed in a separate tile for easier
% comparison, with axes and labels appropriately set for clear visualization.
%
% INPUTS:
%   kValues           - Vector of steering gain values (k) used in the experiment.
%   metrics_prob      - Matrix containing probability near setpoint for each condition (columns).
%   metrics_var       - Matrix containing circular variance values for each condition (columns).
%   metrics_ISE       - Matrix containing Integral of Squared Error (ISE) values for each condition (columns).
%   metrics_IAE       - Matrix containing Integral of Absolute Error (IAE) values for each condition (columns).
%   comparisonLabel   - Cell array of labels for the two conditions being compared (e.g., {'Condition1', 'Condition2'}).
%   folder            - Folder directories
%
% OUTPUTS:
%   Generates a figure with four summary plots showing metrics for each condition across k values.
%
% CREATED: 10/30/2024 - MC
%

function generateSummaryMetrics(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, comparisonLabel, folder)

    % Generate summary plots
    figure;
    set(gcf, 'Position', [100 100 400 800]);  % Set figure size
    tiledlayout(4, 1, 'TileSpacing', 'compact');  % One row for each metric

    % Plot 1: Probability near setpoint across k values
    nexttile;
    plot(kValues, metrics_prob(:, 1), '-o', 'DisplayName', comparisonLabel{1});
    hold on;
    plot(kValues, metrics_prob(:, 2), '-o', 'DisplayName', comparisonLabel{2});
    xlabel('Steering Gain (k)');
    ylabel('Probability Near Setpoint');
    legend('show');
    title(['Probability Near Setpoint (' comparisonLabel{1} ' vs ' comparisonLabel{2} ')']);
    axis padded
    ylim([0 0.5])
    grid on;

    % Plot 2: Circular variance across k values
    nexttile;
    plot(kValues, metrics_var(:, 1), '-o', 'DisplayName', comparisonLabel{1});
    hold on;
    plot(kValues, metrics_var(:, 2), '-o', 'DisplayName', comparisonLabel{2});
    xlabel('Steering Gain (k)');
    ylabel('1 - Circular Variance');
    title(['Circular Variance (' comparisonLabel{1} ' vs ' comparisonLabel{2} ')']);
    axis padded
    ylim([0 1])
    grid on;

    % Plot 3: Integral of Squared Error (ISE) across k values
    nexttile;
    plot(kValues, metrics_ISE(:, 1), '-o', 'DisplayName', comparisonLabel{1});
    hold on;
    plot(kValues, metrics_ISE(:, 2), '-o', 'DisplayName', comparisonLabel{2});
    xlabel('Steering Gain (k)');
    ylabel('ISE');
    title(['Squared Error (' comparisonLabel{1} ' vs ' comparisonLabel{2} ')']);
    axis padded
    grid on;

    % Plot 4: Integral of Absolute Error (IAE) across k values
    nexttile;
    plot(kValues, metrics_IAE(:, 1), '-o', 'DisplayName', comparisonLabel{1});
    hold on;
    plot(kValues, metrics_IAE(:, 2), '-o', 'DisplayName', comparisonLabel{2});
    xlabel('Steering Gain (k)');
    ylabel('IAE');
    title(['Absolute Error (' comparisonLabel{1} ' vs ' comparisonLabel{2} ')']);
    axis padded
    grid on;

    sgtitle(['Summary of Model Performance Across Metrics (' comparisonLabel{1} ' vs ' comparisonLabel{2} ')']);
    % Save summary metrics plot as PNG for final presentation and SVG for high-quality vector graphics
    saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_SummaryMetrics.png']));
    set(gcf, 'renderer', 'Painters');  % Set renderer to 'Painters' for better vector graphic rendering
    saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_SummaryMetrics.svg']));
end

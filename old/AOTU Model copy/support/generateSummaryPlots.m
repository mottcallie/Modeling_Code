% generateSummaryPlots
% This function generates summary plots to compare performance metrics 
% across different steering gains (k values) for two conditions. It includes 
% visualizations for probability near setpoint, circular variance, 
% integral of squared error (ISE), and integral of absolute error (IAE).
%
% INPUTS:
%   kValues - Array of steering gain values to analyze.
%   metrics_prob - Array of probability metrics for different conditions.
%   metrics_var - Array of circular variance metrics for different conditions.
%   metrics_ISE - Array of integral of squared error metrics for different conditions.
%   metrics_IAE - Array of integral of absolute error metrics for different conditions.
%   obj_binned_avg1 - Binned average object positions for condition 1.
%   obj_binned_avg2 - Binned average object positions for condition 2.
%   rotvel_binned_avg1 - Binned average rotational velocities for condition 1.
%   rotvel_binned_avg2 - Binned average rotational velocities for condition 2.
%   evt_1 - Angular velocity data for condition 1.
%   evt_2 - Angular velocity data for condition 2.
%   obj_bins - Binning edges for object positions.
%   rotvel_bins - Binning edges for rotational velocities.
%   posBins - Binning edges for positions.
%   nK - Number of k values.
%   comparisonLabel - Cell array of labels for the two conditions being compared.
%   folder - Structure containing paths for saving the summary plots.
%
% OUTPUTS:
%   Generates and saves various summary plots comparing metrics across the two conditions.
%
function generateSummaryPlots(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, ...
    obj_binned_avg1, obj_binned_avg2, rotvel_binned_avg1, rotvel_binned_avg2, ...
    evt_1, evt_2, obj_bins, rotvel_bins, posBins, nK, comparisonLabel, folder)

%% Generate summary plots
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
ylabel('1- Circular Variance');
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

% Save the summary plot comparing the metrics across the comparison
sgtitle(['Summary of Model Performance Across Metrics (' comparisonLabel{1} ' vs ' comparisonLabel{2} ')']);
saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_SummaryMetrics.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_SummaryMetrics.svg']));

%% Generate two summary plots for binned averages vs. object position and rotational velocity
figure;
set(gcf, 'Position', [100 100 1500 400]);  % Set figure size
tiledlayout(2, nK, 'TileSpacing', 'compact');  % One row for each k value
for kIdx = 1:nK
    nexttile; hold on;
    % Plot binned averages for the two conditions
    this_difference = obj_binned_avg1(kIdx, :) - obj_binned_avg2(kIdx, :);
    plot(obj_bins(1:end-1), this_difference, '-o', 'DisplayName', [comparisonLabel{1} '-' comparisonLabel{2}]);
    xlabel('Object Max (deg)');
    ylabel('Avg. Cross Time');
    title(['k=' num2str(kValues(kIdx))]);
    if kIdx == 1
        legend('show','Location','southeast');
    end
    ylim([-30 30])
    xlim([0 150])
    grid on;
end
% Summary Plot 2: Binned Cross Times vs Rotational Velocity
for kIdx = 1:nK
    nexttile; hold on;
    % Plot binned averages for the two conditions
    thisDiff = rotvel_binned_avg1(kIdx, :) - rotvel_binned_avg2(kIdx, :);
    plot(rotvel_bins(1:end-1), thisDiff, '-o', 'DisplayName', [comparisonLabel{1} '-' comparisonLabel{2}]);
    xlabel('Rotational Velocity Max (deg/s)');
    ylabel('Avg. Cross Time');
    title(['k=' num2str(kValues(kIdx))]);
    if kIdx == 1
        legend('show','Location','northeast');
    end
    ylim([-30 30])
    xlim([0 10])
    grid on;
end
sgtitle('Binned Cross Times');
saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_DirChange.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_DirChange.svg']));

%% Generate Summary Plot for Angular Velocity vs Object Position
figure;
set(gcf, 'Position', [100 100 1500 500]);  % Set figure size
tiledlayout(1, nK, 'TileSpacing', 'compact');  % One row for each k value

% Loop over k values for plotting
for kIdx = 1:nK
    nexttile; hold on;
    
    % Plot the relationship for both conditions
    plot(posBins, evt_1(kIdx, :), 'DisplayName', comparisonLabel{1});
    plot(posBins, evt_2(kIdx, :), 'DisplayName', comparisonLabel{2});

    % Customize the plot
    xlabel('Position (deg)');
    ylabel('Avg. Angular Velocity (deg/s)');
    title(['k=' num2str(kValues(kIdx))]);
    if kIdx == 1
        legend('show', 'Location', 'northwest');
    end
    grid on;
    ylim([-10 10])
    xlim([-100 100])
end

% Save the summary plot
sgtitle(['Binned Angular Velocity vs Object Position for ' comparisonLabel{1} ' vs ' comparisonLabel{2} ' Across k Values']);
saveas(gcf, fullfile(folder.final, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_EVT.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, [comparisonLabel{1} 'v' comparisonLabel{2} '_Binned_EVT.svg']));

end

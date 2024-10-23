function generateSummaryPlot(kValues, prob, var, ISE, IAE, dirChangeTime, thisShift, thisAlpha, folder)
    % generateSummaryPlot
    % Generates a summary plot comparing strength 1 vs strength 0 across k values.
    % kValues: Array of steering gain (k) values
    % prob, var, ISE, IAE, dirChangeTime: Performance metrics
    % thisShift: Current shift value being simulated
    % thisAlpha: Current alpha value being simulated
    % folder: Directory where the summary plot will be saved

    % Generate summary plot
    disp(['    Generating summary plot for shift = ', num2str(thisShift), ', alpha = ', num2str(thisAlpha)]);
    
    figure;
    set(gcf, 'Position', [100, 100, 400, 900]);
    tiledlayout(5, 1, 'TileSpacing', 'compact');

    % Calculate the maximum values for ISE and IAE, then pad them by 10%
    maxISE = max(ISE(:));
    maxIAE = max(IAE(:));
    ISE_ylim = [0, maxISE * 1.1]; % Padded ylim for ISE
    IAE_ylim = [0, maxIAE * 1.1]; % Padded ylim for IAE

    % Plot 1: Probability near setpoint across k values, comparing strength 1 vs 0
    nexttile;
    plot(kValues, squeeze(prob(:, 1)), '-o', 'DisplayName', 'Strength 1');
    hold on;
    plot(kValues, squeeze(prob(:, 2)), '-o', 'DisplayName', 'Strength 0');
    axis padded; ylim([0 1]);
    xlabel('Steering Gain (k)');
    ylabel('Prob(0)');
    legend show;
    title('Probability Near Setpoint (Strength 1 vs 0)');

    % Plot 2: Circular variance across k values, comparing strength 1 vs 0
    nexttile;
    plot(kValues, squeeze(var(:, 1)), '-o', 'DisplayName', 'Strength 1');
    hold on;
    plot(kValues, squeeze(var(:, 2)), '-o', 'DisplayName', 'Strength 0');
    axis padded; ylim([0 1]);
    xlabel('Steering Gain (k)');
    ylabel('1-Circular Var');
    legend show;
    title('Circular Variance (Strength 1 vs 0)');

    % Plot 3: Integral of Squared Error (ISE) across k values, comparing strength 1 vs 0
    nexttile;
    plot(kValues, squeeze(ISE(:, 1)), '-o', 'DisplayName', 'Strength 1');
    hold on;
    plot(kValues, squeeze(ISE(:, 2)), '-o', 'DisplayName', 'Strength 0');
    axis padded; ylim(ISE_ylim); % Apply padded ylim for ISE
    xlabel('Steering Gain (k)');
    ylabel('ISE');
    title('Integral of Squared Error (ISE) (Strength 1 vs 0)');

    % Plot 4: Integral of Absolute Error (IAE) across k values, comparing strength 1 vs 0
    nexttile;
    plot(kValues, squeeze(IAE(:, 1)), '-o', 'DisplayName', 'Strength 1');
    hold on;
    plot(kValues, squeeze(IAE(:, 2)), '-o', 'DisplayName', 'Strength 0');
    axis padded; ylim(IAE_ylim); % Apply padded ylim for IAE
    xlabel('Steering Gain (k)');
    ylabel('IAE');
    title('Integral of Absolute Error (IAE) (Strength 1 vs 0)');

    % Plot 5: Average direction change time across k values, comparing strength 1 vs 0
    nexttile;
    plot(kValues, squeeze(dirChangeTime(:, 1)), '-o', 'DisplayName', 'Strength 1');
    hold on;
    plot(kValues, squeeze(dirChangeTime(:, 2)), '-o', 'DisplayName', 'Strength 0');
    axis padded; ylim([0 max(dirChangeTime(:)) * 1.1]); % Padded ylim for direction change times
    xlabel('Steering Gain (k)');
    ylabel('Avg. Direction Change Time (ms)');
    title('Direction Change Time (Strength 1 vs 0)');

    % Save the summary plot for the current alpha/shift combination
    sgtitle(['Summary for Shift = ' num2str(thisShift) ', Alpha = ' num2str(thisAlpha) ' (Strength 1 vs 0)']);
    saveas(gcf, fullfile(folder.summary, ['Alpha_' num2str(thisAlpha) '_Shift_' num2str(thisShift) '_K_Comparison_Strength1_vs_0.png']));
end
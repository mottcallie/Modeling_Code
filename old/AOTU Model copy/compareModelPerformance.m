% compareModelPerformance
% This function evaluates the performance of a steering model by running 
% simulations across various steering gains (k values) and strengths. It 
% generates visualizations to compare metrics such as probability near 
% setpoint, circular variance, integral of squared error (ISE), and 
% integral of absolute error (IAE) for two strengths.
%
% INPUTS:
%   alpha - Alpha value representing a parameter in the model.
%   predicted_RF - Predicted receptive fields to use for the model.
%
% OUTPUTS:
%   Generates and saves various plots comparing the model's performance 
%   across different metrics and strengths.
%
% Created:  [Date] - MC
%
function compareModelPerformance(alpha,predicted_RF)
%% Initialize

% Refresh settings
[folder, plotSettings, runSettings] = modelSettings();
runSettings.numRuns = 100; %reduced runs
close all;

% Define the range to test
kValues = 8:2:20;  % Steering gain (k) values to iterate over
nK = length(kValues);  % Number of k values
strengthValues = [1, 0];
noiseLevel = 2.5;
startPos = 0;
simDuration = 30;
thisSynapse = "inhibitory";  % Fixed synapse type

% Convert alpha value to a plot-compatible string
alphaStr = sprintf('a%s', strrep(num2str(alpha, '%.2f'), '.', ''));  % E.g., 0.5 becomes a05

%% Plot RFs and ELU for reference
figure; set(gcf, 'Position', [100 100 900 400]);  % Set figure size
tiledlayout(1,2, 'TileSpacing', 'compact'); 
nexttile
plotRFsTogether(predicted_RF,runSettings)
nexttile
plot_adjELU(alpha, runSettings)
sgtitle(alphaStr)
% Save the plot
saveas(gcf, fullfile(folder.summary, ['Test_Model_Settings_' alphaStr '.png']));

%% Run model

% Preallocate to store performance metrics and binned averages
metrics_prob = zeros(nK, 2);  % Two columns for strength 1 and strength 0
metrics_var = zeros(nK, 2);
metrics_ISE = zeros(nK, 2);
metrics_IAE = zeros(nK, 2);
metrics_dirChangeTime = zeros(nK, 2);

% Preallocate to store binned averages
obj_binned_avg_strength1 = [];
obj_binned_avg_strength0 = [];
rotvel_binned_avg_strength1 = [];
rotvel_binned_avg_strength0 = [];
evt_strength1 = [];
evt_strength0 = [];

% Select up to X evenly spaced k values for plotting
maxPlotK = min(6, nK);  % Limit to X `k` runs for plotting
selected_k_indices = round(linspace(1, nK, maxPlotK));  % Evenly spaced indices from `kValues`
strengthColors = {"#0072BD";"#D95319"};
% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);  % Set figure size
tiledlayout(maxPlotK, 2, 'TileSpacing', 'compact'); 

% Loop over k values (steering gain) for analysis
for kIdx = 1:nK
    thisK = kValues(kIdx);
    % Output progress
    disp(['Running simulations for k = ', num2str(thisK)]);

    % Loop over strength values (1 and 0)
    for strengthIdx = 1:2
        thisStrength = strengthValues(strengthIdx);

        % Adjust tuning for this strength
        thisTuning = predicted_RF;
        thisTuning.AOTU019 = thisTuning.AOTU019 .* thisStrength;  % Adjust tuning by strength

        % Run the AOTU steering model with adjusted tuning
        [timebase, visobj_history, input_history, rotvel_history] = ...
            aotu_steering_model(thisTuning, noiseLevel, startPos, thisSynapse, simDuration, thisK, alpha, runSettings);

        % Calculate performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        % Store metrics for this strength and k value
        metrics_prob(kIdx, strengthIdx) = metrics_results.prob;
        metrics_var(kIdx, strengthIdx) = metrics_results.var;
        metrics_ISE(kIdx, strengthIdx) = metrics_results.ISE;
        metrics_IAE(kIdx, strengthIdx) = metrics_results.IAE;

        % Measure performance at direction changes and store binned averages
        nTest = 100;
        [rotvel_cross_times, obj_bins, rotvel_bins, obj_binned_avgs, rotvel_binned_avgs] = analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest);
        % Store direction changes
        metrics_dirChangeTime(kIdx, strengthIdx) = mean(rotvel_cross_times, 'omitnan');

        % analyze object position v velocity
        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings); 

        % Store the binned averages for strength 1 and strength 0 separately
        if strengthIdx == 1
            obj_binned_avg_strength1(kIdx, :) = obj_binned_avgs;
            rotvel_binned_avg_strength1(kIdx, :) = rotvel_binned_avgs;

            evt_strength1(kIdx, :) = posvang;
        else
            obj_binned_avg_strength0(kIdx, :) = obj_binned_avgs;
            rotvel_binned_avg_strength0(kIdx, :) = rotvel_binned_avgs;

            evt_strength0(kIdx, :) = posvang;
        end

        % Plot example run
        if any(kIdx==selected_k_indices)
            nexttile;
            plot(timebase, visobj_history(1, :),'Color',strengthColors{strengthIdx});
            xlabel('Time (s)');
            ylabel('Pos');
            title(['k = ' num2str(thisK)]);
            axis tight;
            if kIdx==1
                legend([num2str(strengthValues(strengthIdx)) 'X'])
            end
            ylim([-180 180]);
        end
    end
end

% Save the example plot comparing strength 1 and 0 across selected k values
sgtitle([(alphaStr) ' Example Runs for Strength 1 vs Strength 0']);
saveas(gcf, fullfile(folder.summary, ['Test_ExampleRuns_' alphaStr '.png']));

%% Generate summary plots
% Summary Plot Comparing Across Different Metrics (Strength 0 vs Strength 1)
figure;
set(gcf, 'Position', [100 100 400 800]);  % Set figure size
tiledlayout(4, 1, 'TileSpacing', 'compact');  % One row for each metric

% Plot 1: Probability near setpoint across k values
nexttile;
plot(kValues, metrics_prob(:, 1), '-o', 'DisplayName', 'Strength 1');
hold on;
plot(kValues, metrics_prob(:, 2), '-o', 'DisplayName', 'Strength 0');
xlabel('Steering Gain (k)');
ylabel('Probability Near Setpoint');
legend('show');
title('Probability Near Setpoint (Strength 1 vs Strength 0)');
axis padded
ylim([0 0.5])
grid on;

% Plot 2: Circular variance across k values
nexttile;
plot(kValues, metrics_var(:, 1), '-o', 'DisplayName', 'Strength 1');
hold on;
plot(kValues, metrics_var(:, 2), '-o', 'DisplayName', 'Strength 0');
xlabel('Steering Gain (k)');
ylabel('1- Circular Variance');
title('Circular Variance (Strength 1 vs Strength 0)');
axis padded
ylim([0 1])
grid on;

% Plot 3: Integral of Squared Error (ISE) across k values
nexttile;
plot(kValues, metrics_ISE(:, 1), '-o', 'DisplayName', 'Strength 1');
hold on;
plot(kValues, metrics_ISE(:, 2), '-o', 'DisplayName', 'Strength 0');
xlabel('Steering Gain (k)');
ylabel('ISE');
title('Squared Error (Strength 1 vs Strength 0)');
axis padded
grid on;

% Plot 4: Integral of Absolute Error (IAE) across k values
nexttile;
plot(kValues, metrics_IAE(:, 1), '-o', 'DisplayName', 'Strength 1');
hold on;
plot(kValues, metrics_IAE(:, 2), '-o', 'DisplayName', 'Strength 0');
xlabel('Steering Gain (k)');
ylabel('IAE');
title('Absolute Error (Strength 1 vs Strength 0)');
axis padded
grid on;

% Save the summary plot comparing strength 0 and 1 across metrics
sgtitle({[alphaStr 'Summary of Model Performance'],'Across Metrics (Strength 1 vs Strength 0)'});
saveas(gcf, fullfile(folder.summary, ['Test_SummaryMetrics_' alphaStr '.png']));

% Generate two summary plots for binned averages vs. object position and rotational velocity
% Summary Plot 1: Binned Cross Times vs Object Position
figure;
set(gcf, 'Position', [100 100 1500 400]);  % Set figure size
tiledlayout(2, nK, 'TileSpacing', 'compact');  % One row for each k value
for kIdx = 1:nK
    nexttile; hold on;
    % Plot binned averages for strength 1 and 0
    this_difference = obj_binned_avg_strength0(kIdx, :) - obj_binned_avg_strength1(kIdx, :);
    plot(obj_bins(1:end-1), this_difference, '-o', 'DisplayName', '0-1');
    xlabel('Object Max (deg)');
    ylabel('Avg. Cross Time');
    title(['k=' num2str(kValues(kIdx))]);
    if kIdx==1
        legend('show','Location','southeast');
    end
    ylim([00 30])
    xlim([0 150])
    grid on;
end
% Summary Plot 1: Binned Cross Times vs Object Position
for kIdx = 1:nK
    nexttile; hold on;
    % Plot binned averages for strength 1 and 0
    thisDiff = rotvel_binned_avg_strength0(kIdx, :) - rotvel_binned_avg_strength1(kIdx, :);
    plot(rotvel_bins(1:end-1), thisDiff, '-o', 'DisplayName', '0-1');
    xlabel('Rotational Velocity Max (deg/s)');
    ylabel('Avg. Cross Time');
    title(['k=' num2str(kValues(kIdx))]);
    if kIdx==1
        legend('show','Location','southeast');
    end
    ylim([00 30])
    xlim([0 25])
    grid on;
end
sgtitle([alphaStr 'Binned Cross Times']);
saveas(gcf, fullfile(folder.summary, ['Test_Binned_DirChange_' alphaStr '.png']));


% Generate Summary Plot for Strength 0 vs Strength 1
figure;
set(gcf, 'Position', [100 100 1500 500]);  % Set figure size
tiledlayout(1, nK, 'TileSpacing', 'compact');  % One row for each k value

% Loop over k values for plotting
for kIdx = 1:nK
    nexttile; hold on;
    
    % Plot the relationship for strength 1 and strength 0
    plot(posBins, evt_strength1(kIdx, :), 'DisplayName', 'Strength 1');
    plot(posBins, evt_strength0(kIdx, :), 'DisplayName', 'Strength 0');

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
sgtitle([alphaStr 'Binned Angular Velocity vs Object Position for Strength 1 vs Strength 0 Across k Values']);
saveas(gcf, fullfile(folder.summary, ['Test_EVT_' alphaStr '.png']));
end
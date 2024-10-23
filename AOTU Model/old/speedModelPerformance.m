function speedModelPerformance(alpha,predicted_RF)
%% Initialize

% Refresh settings
[folder, plotSettings, runSettings] = modelSettings();
close all;

% Define the range to test
kValues = 8:2:20;  % Steering gain (k) values to iterate over
nK = length(kValues);  % Number of k values
noiseLevel = 2.5;
startPos = 0;
simDuration = 30;
thisSynapse = "inhibitory";

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
saveas(gcf, fullfile(folder.final, ['FastvSlow_Model_Settings_' alphaStr '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['FastvSlow_Model_Settings_' alphaStr '.svg']));

%% Run model

% Preallocate to store performance metrics and binned averages
metrics_prob = zeros(nK, 2);  % Two columns for fast and slow
metrics_var = zeros(nK, 2);
metrics_ISE = zeros(nK, 2);
metrics_IAE = zeros(nK, 2);
metrics_dirChangeTime = zeros(nK, 2);

% Preallocate to store binned averages
obj_binned_avg_fast = [];
obj_binned_avg_slow = [];
rotvel_binned_avg_fast = [];
rotvel_binned_avg_slow = [];
evt_fast = [];
evt_slow = [];

% Define visuomotor delays (fast and slow)
visuomotorDelays = [80, 100];  % Fast = 80 ms, Slow = 100 ms
speedColors = {"#0072BD"; "#D95319"};  % Colors for fast and slow

% Select up to X evenly spaced k values for plotting
maxPlotK = min(6, nK);  % Limit to X `k` runs for plotting
selected_k_indices = round(linspace(1, nK, maxPlotK));  % Evenly spaced indices from `kValues`

% Initialize figure for example plots
figure; set(gcf, 'Position', [100 100 1500 900]);  % Set figure size
tiledlayout(maxPlotK, 2, 'TileSpacing', 'compact');

% Loop over k values (steering gain) for analysis
for kIdx = 1:nK
    thisK = kValues(kIdx);
    % Output progress
    disp(['Running simulations for k = ', num2str(thisK)]);

    % Loop over visuomotor delays (fast and slow)
    for speedIdx = 1:2
        thisDelay = visuomotorDelays(speedIdx);
        runSettings.AOTU019_delay = thisDelay;

        % Run the AOTU steering model with adjusted delay
        [timebase, visobj_history, input_history, rotvel_history] = ...
            aotu_steering_model(predicted_RF, noiseLevel, startPos, thisSynapse, simDuration, thisK, alpha, runSettings);

        % Calculate performance metrics
        metrics_results = calculatePerformanceMetrics(visobj_history, rotvel_history, timebase, plotSettings);
        % Store metrics for this speed (fast/slow) and k value
        metrics_prob(kIdx, speedIdx) = metrics_results.prob;
        metrics_var(kIdx, speedIdx) = metrics_results.var;
        metrics_ISE(kIdx, speedIdx) = metrics_results.ISE;
        metrics_IAE(kIdx, speedIdx) = metrics_results.IAE;

        % Measure performance at direction changes and store binned averages
        nTest = 100;
        [rotvel_cross_times, obj_bins, rotvel_bins, obj_binned_avgs, rotvel_binned_avgs] = analyzeObjectAndVelocityCrossings(visobj_history, rotvel_history, timebase, nTest);
        metrics_dirChangeTime(kIdx, speedIdx) = mean(rotvel_cross_times, 'omitnan');

        % Analyze object position v velocity
        [posvang, posBins] = analyzeErrorVsTurn(visobj_history, rotvel_history, runSettings); 

        % Store the binned averages for fast and slow delays separately
        if thisDelay == 80
            obj_binned_avg_fast(kIdx, :) = obj_binned_avgs;
            rotvel_binned_avg_fast(kIdx, :) = rotvel_binned_avgs;
            evt_fast(kIdx, :) = posvang;
        else
            obj_binned_avg_slow(kIdx, :) = obj_binned_avgs;
            rotvel_binned_avg_slow(kIdx, :) = rotvel_binned_avgs;
            evt_slow(kIdx, :) = posvang;
        end

        % Plot example run
        if any(kIdx == selected_k_indices)
            nexttile;
            plot(timebase, visobj_history(1, :), 'Color', speedColors{speedIdx});
            xlabel('Time (s)');
            ylabel('Pos');
            title(['k = ' num2str(thisK)]);
            axis tight;
            if kIdx == 1
                legend(['Speed = ', num2str(thisDelay), ' ms']);
            end
            ylim([-180 180]);
        end
    end
end

% Save the example plot comparing fast and slow delays across selected k values
sgtitle([alphaStr ' Example Runs for Fast vs Slow Delays']);
saveas(gcf, fullfile(folder.final, ['FastvSlow_ExampleRuns_' alphaStr '.png']));
set(gcf, 'renderer', 'Painters');
saveas(gcf, fullfile(folder.vectors, ['FastvSlow_ExampleRuns_' alphaStr '.svg']));

%% Generate summary
% Call the generateSummaryPlots function after the model run loop
comparisonLabel = {'Fast', 'Slow'};

% Generate the summary plots
generateSummaryPlots(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, ...
    obj_binned_avg_fast, obj_binned_avg_slow, rotvel_binned_avg_fast, rotvel_binned_avg_slow, ...
    evt_fast, evt_slow, obj_bins, rotvel_bins, posBins, nK, alphaStr, comparisonLabel, folder);

end
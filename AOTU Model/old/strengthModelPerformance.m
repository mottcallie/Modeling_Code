function strengthModelPerformance(alpha,predicted_RF)
%% Initialize

% Refresh settings
[folder, plotSettings, runSettings] = modelSettings();
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
saveas(gcf, fullfile(folder.final, ['Strength1vStrength0_Model_Settings_' alphaStr '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['Strength1vStrength0_Model_Settings_' alphaStr '.svg']));

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
saveas(gcf, fullfile(folder.final, ['Strength1vStrength0_ExampleRuns_' alphaStr '.png']));
set(gcf,'renderer','Painters')
saveas(gcf, fullfile(folder.vectors, ['Strength1vStrength0_ExampleRuns_' alphaStr '.svg']));

%% Generate summary plots
% Call the generateSummaryPlots function after the model run loop
comparisonLabel = {'Strength1', 'Strength0'};

% Generate the summary plots
generateSummaryPlots(kValues, metrics_prob, metrics_var, metrics_ISE, metrics_IAE, ...
    obj_binned_avg_strength1, obj_binned_avg_strength0, rotvel_binned_avg_strength1, rotvel_binned_avg_strength0, ...
    evt_strength1, evt_strength0, obj_bins, rotvel_bins, posBins, nK, alphaStr, comparisonLabel, folder);

end